import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
from redolarium.config import CONFIG

def scan_gene_promoter(sequence, gene_start, gene_end, strand, sigma_dict, taxonomy_context, window=600):
    seq_len = len(sequence)
    promoters = []
    
    # Define upstream region based on gene strand
    if strand == 1 or strand == "+":
        up_start = max(0, gene_start - window)
        up_end = gene_start
        up_seq_obj = Seq(sequence[up_start:up_end])
        up_seq = str(up_seq_obj).upper()
        offset_base = up_start
    else:
        up_start = gene_end
        up_end = min(seq_len, gene_end + window)
        up_seq_obj = Seq(sequence[up_start:up_end]).reverse_complement()
        up_seq = str(up_seq_obj).upper()
        offset_base = up_end

    # Scan for Shine-Dalgarno in the 30bp immediately preceding the start codon
    sd_motif = "Not detected"
    sd_pattern = r"(A?GGAGG|GGAG|GAGG)"
    sd_search_region = up_seq[-30:] if len(up_seq) >= 30 else up_seq
    sd_match = re.search(sd_pattern, sd_search_region)
    if sd_match:
        sd_motif = sd_match.group(0)
    
    # Threshold for PSSM scoring (log odds)
    threshold = 3.0
    
    for sigma_name, constraints in sigma_dict.items():
        inst_35 = constraints.get("instances_35", [])
        inst_10 = constraints.get("instances_10", [])
        s_min = constraints.get("spacer_min", 15)
        s_max = constraints.get("spacer_max", 19)
        induction = constraints.get("induction", "Unknown")
        
        if not inst_35 or not inst_10:
            continue
            
        try:
            m35 = motifs.create([Seq(s) for s in inst_35])
            pssm35 = m35.counts.normalize(pseudocounts=0.5).log_odds()
            len35 = len(inst_35[0])
            
            m10 = motifs.create([Seq(s) for s in inst_10])
            pssm10 = m10.counts.normalize(pseudocounts=0.5).log_odds()
            len10 = len(inst_10[0])
        except Exception:
            continue # Bio.motifs failed to build
            
        up_seq_len = len(up_seq)
        if up_seq_len < (len35 + s_max + len10):
            continue
            
        # Calculate scores along the sequence
        # Only compute up to where a full motif fits
        scores35 = [pssm35.calculate(Seq(up_seq[i:i+len35])) for i in range(up_seq_len - len35 + 1)]
        scores10 = [pssm10.calculate(Seq(up_seq[i:i+len10])) for i in range(up_seq_len - len10 + 1)]
        
        for i in range(len(scores35)):
            if scores35[i] < -2.0: # Skip very poor -35 boxes to save time
                continue
                
            for spacer in range(s_min, s_max + 1):
                j = i + len35 + spacer
                if j < len(scores10):
                    total_score = scores35[i] + scores10[j]
                    if total_score > threshold:
                        box35 = up_seq[i:i+len35]
                        box10 = up_seq[j:j+len10]
                        
                        if strand == 1 or strand == "+":
                            loc_offset = offset_base + i
                        else:
                            loc_offset = offset_base - (i + len35 + spacer + len10)
                            
                        quality = "Strong" if sd_motif != "Not detected" else "Weak"
                        
                        promoters.append({
                            "Sigma_Factor": f"{sigma_name} [{taxonomy_context}]",
                            "Minus35_Seq": box35,
                            "Minus10_Seq": box10,
                            "Spacer_Length": spacer,
                            "Upstream_Position": loc_offset,
                            "Shine_Dalgarno": sd_motif,
                            "Quality_Class": quality,
                            "Induction_Condition": induction,
                            "Score": round(float(total_score), 2)
                        })
                        
    # Sort to keep highest scores
    promoters = sorted(promoters, key=lambda x: x["Score"], reverse=True)
    return promoters

def run_promoter_prediction(query_gb, target_bgc, region_genes, out_dir, logger):
    """
    Genome-Aware Promoter Identification using PSSM logic to detect degenerate BGC promoters.
    Also predicts BGC Operon Architecture globally.
    """
    logger.info("Stage 10: Genome-Aware Promoter Identification and PSSM Extraction...")
    
    bgc_id = target_bgc["BGC_ID"]
    
    try:
        record = SeqIO.read(query_gb, "genbank")
    except Exception:
        try:
            record = SeqIO.read(query_gb, "fasta")
        except Exception as e:
            logger.warning(f"Failed to read genome for promoter prediction: {e}")
            return []
            
    sequence = str(record.seq)
    
    taxonomy_list = record.annotations.get("taxonomy", [])
    taxonomy_str = " ".join(taxonomy_list).lower()
    
    all_motifs = CONFIG.get("sigma_motifs", {})
    selected_phylum = "Universal"
    
    if "firmicutes" in taxonomy_str or "bacillota" in taxonomy_str:
        selected_phylum = "Firmicutes"
    elif "proteobacteria" in taxonomy_str or "pseudomonadota" in taxonomy_str:
        selected_phylum = "Pseudomonadota"
    elif "actinobacteria" in taxonomy_str or "actinomycetota" in taxonomy_str:
        selected_phylum = "Actinomycetota"
    elif "bacteroidetes" in taxonomy_str or "bacteroidota" in taxonomy_str:
        selected_phylum = "Bacteroidota"
    
    sigma_dict = all_motifs.get(selected_phylum, all_motifs.get("Universal", {}))
    logger.info(f"Extracted organism taxonomy mapped to: {selected_phylum}. Applying {selected_phylum} motif PSSMs.")
    
    promoter_records = []
    
    for gene in region_genes:
        locus_tag = gene.get("Locus_Tag", "Unknown")
        gene_symbol = gene.get("Gene", "-")
        start = gene.get("Start", 0)
        end = gene.get("End", 0)
        strand = gene.get("Strand", "+")
        role = gene.get("Role", "Other")
        
        if start == 0 and end == 0:
            continue
            
        search_window = CONFIG.get("promoter_search_upstream", 600)
        gene_promoters = scan_gene_promoter(sequence, start, end, strand, sigma_dict, selected_phylum, window=search_window)
        
        # Only keep top 2 non-redundant promoters per gene to avoid noise
        if len(gene_promoters) > 2:
            gene_promoters = gene_promoters[:2]
            
        for p in gene_promoters:
            record_dict = {
                "Locus_Tag": locus_tag,
                "Gene_Symbol": gene_symbol,
                "Role": role
            }
            record_dict.update(p)
            promoter_records.append(record_dict)
            
    # Calculate Operon Architecture Globally
    architecture = "Unknown (No promoters found)"
    if promoter_records:
        core_genes = [g for g in region_genes if g.get("Role") == "Core Biosynthetic"]
        if core_genes:
            # Sort core genes by genomic position
            core_genes_sorted = sorted(core_genes, key=lambda x: min(x["Start"], x["End"]))
            first_core_locus = core_genes_sorted[0]["Locus_Tag"]
            
            promoters_on_core = [p for p in promoter_records if p["Role"] == "Core Biosynthetic"]
            core_promoter_loci = set([p["Locus_Tag"] for p in promoters_on_core])
            
            if len(core_promoter_loci) == 1 and first_core_locus in core_promoter_loci:
                architecture = "Single Operon Driven (Promoter at Core Start)"
            elif len(core_promoter_loci) > 1:
                architecture = "Independently Regulated Modules (Internal Promoters)"
            elif len(core_promoter_loci) == 0:
                architecture = "Putative Single Operon (Promoter on Flanking Regulator)"
                
    # Append Operon Architecture to all records
    for r in promoter_records:
        r["Operon_Architecture"] = architecture

    os.makedirs(os.path.join(out_dir, "bgc_motifs"), exist_ok=True)
    csv_out = os.path.join(out_dir, "bgc_motifs", f"{bgc_id}_promoter_motifs.csv")
    
    if not promoter_records:
        logger.info(f"No promoter motifs detected in BGC {bgc_id} (threshold > 3.0).")
        df_prom = pd.DataFrame(columns=[
            "Locus_Tag", "Gene_Symbol", "Sigma_Factor", "Minus35_Seq", "Minus10_Seq", 
            "Spacer_Length", "Upstream_Position", "Score", "Shine_Dalgarno", "Quality_Class", 
            "Induction_Condition", "Operon_Architecture"
        ])
    else:
        logger.info(f"Detected {len(promoter_records)} putative regulatory elements in BGC {bgc_id}. Architecture: {architecture}")
        df_prom = pd.DataFrame(promoter_records)
        # Drop internal processing column
        if "Role" in df_prom.columns:
            df_prom = df_prom.drop(columns=["Role"])
        
    df_prom.to_csv(csv_out, index=False)
    logger.info(f"Saved targeted BGC promoter motifs data to: {csv_out}")
    
    return promoter_records
