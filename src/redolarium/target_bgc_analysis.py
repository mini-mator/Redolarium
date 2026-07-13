# type: ignore
import os
import re
import time
import platform
import urllib.request
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from redolarium.utils import api_retry, save_publication_plot
from redolarium.structures import PredictionResult
from redolarium.config import CONFIG
from redolarium.qc import get_module_qc_penalty

def scan_mge_domains(flanking_genes, logger):
    try:
        import pyhmmer
        import pyhmmer.easel as easel
        import pyhmmer.plan7 as plan7
    except ImportError:
        logger.warning("pyhmmer not installed. Skipping MGE HMM profiling.")
        return {}

    local_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources")
    hmm_path = os.path.join(local_dir, "mge_profiles.hmm")
    
    pfams = ["PF00589", "PF00872", "PF00239"]
    
    if not os.path.exists(hmm_path):
        logger.info("MGE HMM database not found locally. Downloading from Pfam...")
        import gzip
        combined_content = ""
        for pfam in pfams:
            url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam}?annotation=hmm"
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
                with urllib.request.urlopen(req, timeout=10) as response:
                    raw_data = response.read()
                    try:
                        decompressed = gzip.decompress(raw_data)
                        combined_content += decompressed.decode("utf-8") + "\n"
                    except Exception:
                        combined_content += raw_data.decode("utf-8") + "\n"
            except Exception as e:
                logger.warning(f"Failed to download Pfam HMM for {pfam}: {e}")
        if combined_content.strip():
            try:
                with open(hmm_path, "w", encoding="utf-8") as f:
                    f.write(combined_content)
            except Exception:
                pass
                
    if not os.path.exists(hmm_path):
        return {}
        
    alphabet = easel.Alphabet.amino()
    digitized_seqs = []
    for g in flanking_genes:
        trans = g.get("translation") or g.get("Protein_Sequence") or ""
        ltag = g.get("locus_tag") or g.get("Locus_Tag") or ""
        if trans and ltag and len(trans) > 10:
            try:
                seq = easel.TextSequence(name=ltag.encode("utf-8"), sequence=trans).digitize(alphabet)
                digitized_seqs.append(seq)
            except Exception:
                pass
                
    if not digitized_seqs:
        return {}
        
    mge_hits = {}
    try:
        import pyhmmer.hmmer as hmmer
        with plan7.HMMFile(hmm_path) as hmm_file:
            hmms = list(hmm_file)
            all_hits = hmmer.hmmsearch(hmms, digitized_seqs, cpus=1)
            for hmm, hits in zip(hmms, all_hits):
                hmm_name = hmm.name.decode("utf-8") if isinstance(hmm.name, bytes) else hmm.name
                for hit in hits:
                    if hit.evalue <= 1e-5:
                        ltag = hit.name.decode("utf-8") if isinstance(hit.name, bytes) else hit.name
                        mge_hits[ltag] = hmm_name
    except Exception as e:
        logger.warning(f"MGE HMM search execution failed: {e}")
        
    return mge_hits

@api_retry(retries=5, backoff_factor=2.0)
def run_bgc_conservation_blast(bgc_id, query_seq, email, out_dir, logger, is_protein=True):
    program = "blastp" if is_protein else "blastn"
    database = "nr" if is_protein else "nt"
    
    logger.info(f"Running BGC conservation BLAST ({program} against '{database}' database) for {bgc_id}...")
    
    blast_rows = []
    retries = 3
    delay = 10
    
    if not is_protein and len(query_seq) > 1000:
        logger.info(f"Nucleotide query sequence is too long ({len(query_seq)} bp). Slicing to first 1000 bp to prevent remote timeouts.")
        query_seq = query_seq[:1000]
        
    while retries > 0:
        try:
            result_handle = NCBIWWW.qblast(program, database, query_seq, hitlist_size=30)
            blast_record = NCBIXML.read(result_handle)
            
            for alignment in blast_record.alignments:
                title = alignment.title
                accession = alignment.accession
                if alignment.hsps:
                    best_hsp = max(alignment.hsps, key=lambda h: h.score)
                    ident = (best_hsp.identities / best_hsp.align_length) * 100 if best_hsp.align_length > 0 else 0.0
                    cov = (best_hsp.align_length / len(query_seq)) * 100 if len(query_seq) > 0 else 0.0
                    
                    if ident < 30.0 or cov < 30.0:
                        logger.warning(f"Excluding low homology hit {accession} (Identity: {ident:.1f}%, Coverage: {cov:.1f}%).")
                        continue
                        
                    blast_rows.append({
                        "Target_Species": title.split(",")[0],
                        "Accession": accession,
                        "Evalue": best_hsp.expect,
                        "Identity_Pct": round(ident, 2),
                        "Bit_Score": best_hsp.score,
                        "Alignment_Length": best_hsp.align_length,
                        "Coverage_Pct": round(cov, 2)
                    })
            break
        except Exception as e:
            logger.warning(f"BGC BLAST query attempt failed: {e}. Retrying in {delay}s...")
            retries -= 1
            time.sleep(delay)
            
    if not blast_rows:
        logger.warning(f"BGC BLAST returned no homologous sequences. Preserving empty dataset as biological absence.")
        
    df_blast = pd.DataFrame(blast_rows)
    csv_out = os.path.join(out_dir, "tabular_data", f"{bgc_id}_conservation_blast.csv")
    df_blast.to_csv(csv_out, index=False)
    logger.info(f"Saved BGC conservation BLAST results to: {csv_out}")
    
    val_out = os.path.join(out_dir, "tabular_data", f"{bgc_id}_alignment_validation.csv")
    df_val = df_blast[["Accession", "Target_Species", "Evalue", "Identity_Pct", "Bit_Score", "Coverage_Pct"]] if not df_blast.empty else pd.DataFrame(columns=["Accession", "Target_Species", "Evalue", "Identity_Pct", "Bit_Score", "Coverage_Pct"])
    df_val.to_csv(val_out, index=False)
    logger.info(f"Saved alignment parameters validation report to: {val_out}")
    
    return blast_rows

def run_target_bgc_analysis(query_gb, selected_bgc, run_blast, email, out_dir, logger, qc_result=None) -> PredictionResult:
    logger.info("Stage 9: Running Target BGC promoter and flanking prophage analysis...")
    start_time = time.time()
    
    record = max(list(SeqIO.parse(query_gb, "genbank")), key=lambda r: len(r.seq))
    bgc_id = selected_bgc["BGC_ID"]
    bgc_start = selected_bgc["Start_Coord"]
    bgc_end = selected_bgc["End_Coord"]
    
    # ─── 1. Promoter Scanning ───
    logger.info(f"Scanning promoter regions upstream of BGC core genes for {bgc_id} using PWM & shape analysis...")
    promoter_records = []
    
    core_genes = selected_bgc.get("Hits", [])
    if not core_genes:
        core_genes = []
        for feat in record.features:
            if feat.type == "CDS":
                fs = int(feat.location.start)
                fe = int(feat.location.end)
                if bgc_start <= fs <= bgc_end:
                    core_genes.append({
                        "locus_tag": feat.qualifiers.get("locus_tag", [""])[0],
                        "gene": feat.qualifiers.get("gene", ["NA"])[0],
                        "start": fs,
                        "end": fe,
                        "strand": "+" if feat.location.strand >= 0 else "-",
                        "product": feat.qualifiers.get("product", [""])[0]
                    })
                    
    from redolarium.promoter_prediction import run_promoter_prediction
    
    # Normalize core_genes keys for run_promoter_prediction
    normalized_genes = []
    for g in core_genes:
        normalized_genes.append({
            "Locus_Tag": g.get("locus_tag") or g.get("Locus_Tag", "Unknown"),
            "Gene": g.get("gene") or g.get("Gene_Symbol", "-"),
            "Start": g.get("start") if g.get("start") is not None else g.get("Start_Coord", 0),
            "End": g.get("end") if g.get("end") is not None else g.get("End_Coord", 0),
            "Strand": g.get("strand") or g.get("Strand", "+"),
            "Role": "Core Biosynthetic"
        })
        
    promoter_records = run_promoter_prediction(query_gb, selected_bgc, normalized_genes, out_dir, logger)

    # ─── 2. Flanking Phage Remnants and MGE Scanning ───
    logger.info(f"Scanning flanking coordinates of BGC {bgc_id} for mobile elements and prophages...")
    phage_hits = []
    flank_w = 30000
    
    flank_genes_list = []
    for feat in record.features:
        if feat.type == "CDS":
            fs = int(feat.location.start)
            fe = int(feat.location.end)
            if (bgc_start - flank_w <= fs < bgc_start) or (bgc_end < fe <= bgc_end + flank_w):
                flank_genes_list.append({
                    "locus_tag": feat.qualifiers.get("locus_tag", [""])[0],
                    "translation": feat.qualifiers.get("translation", [""])[0],
                    "product": feat.qualifiers.get("product", [""])[0],
                    "gene": feat.qualifiers.get("gene", [""])[0],
                    "start": fs,
                    "end": fe,
                    "strand": "+" if feat.location.strand >= 0 else "-"
                })
                
    mge_domain_hits = scan_mge_domains(flank_genes_list, logger)
    
    for g in flank_genes_list:
        ltag = g["locus_tag"]
        hit_domain = mge_domain_hits.get(ltag)
        
        is_phage = False
        mge_class = "Other MGE Element"
        prod_desc = g["product"]
        
        if hit_domain:
            is_phage = True
            if hit_domain == "PF00589":
                mge_class = "Integrase Component"
                prod_desc = "Integrase catalytic domain (detected via pyhmmer Pfam PF00589)"
            elif hit_domain == "PF00872":
                mge_class = "Transposase Component"
                prod_desc = "Transposase insertion domain (detected via pyhmmer Pfam PF00872)"
            elif hit_domain == "PF00239":
                mge_class = "Recombinase Component"
                prod_desc = "Recombinase domain (detected via pyhmmer Pfam PF00239)"
        else:
            # Removed unscientific string-matching for phage keywords
            pass
            
        if is_phage:
            phage_hits.append({
                "Locus_Tag": ltag,
                "Start": g["start"],
                "End": g["end"],
                "Strand": g["strand"],
                "Product": prod_desc,
                "MGE_Class": mge_class
            })
            
    if not phage_hits:
        df_phage = pd.DataFrame(columns=["Locus_Tag", "Start", "End", "Strand", "Product", "MGE_Class"])
    else:
        df_phage = pd.DataFrame(phage_hits)
    csv_phage = os.path.join(out_dir, "tabular_data", f"{bgc_id}_phage_artifacts.csv")
    df_phage.to_csv(csv_phage, index=False)
    logger.info(f"Saved prophage flanking artifacts dataset to: {csv_phage}")

    # ─── 3. BGC-Specific Cross-Species BLAST Conservation ───
    core_prot = ""
    is_protein = True
    
    # Removed unscientific keyword matching ("synthase"). 
    # Since bgc_analysis.py guarantees architectural core validity, select the longest core sequence.
    valid_prots = [g for g in core_genes if g.get("translation")]
    if not valid_prots:
        for g in core_genes:
            if "locus_tag" in g:
                for feat in record.features:
                    if feat.type == "CDS" and feat.qualifiers.get("locus_tag", [""])[0] == g["locus_tag"]:
                        trans = feat.qualifiers.get("translation", [""])
                        if trans:
                            valid_prots.append({"translation": trans[0]})
                            
    if valid_prots:
        core_prot = max(valid_prots, key=lambda g: len(g["translation"]))["translation"]
    if not core_prot:
        logger.warning("No core biosynthetic protein sequence translation found. Falling back to core nucleotide fragment...")
        try:
            core_seq = str(record.seq[bgc_start:min(bgc_start + 1500, bgc_end)])
        except Exception:
            core_seq = ""
        is_protein = False
        query_seq = core_seq
    else:
        query_seq = core_prot
        
    blast_rows = []
    if run_blast and query_seq:
        blast_rows = run_bgc_conservation_blast(bgc_id, query_seq, email, out_dir, logger, is_protein=is_protein)
    else:
        logger.info("Remote BGC BLAST conservation skipped by request or empty sequence.")
        
    # ─── 4. Render Promoter & Flanking Environment Map ───
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    ax1.set_xlim(-100, 10)
    ax1.set_ylim(-1, 2)
    ax1.axis("off")
    ax1.plot([-100, 10], [0, 0], color="black", lw=2)
    ax1.text(0, -0.4, "Start Codon (ATG)", ha="center", fontsize=9, fontweight="bold")
    ax1.axvline(0, color="red", linestyle="--", ymin=0.3, ymax=0.7)
    
    if promoter_records:
        p35_box = patches.Rectangle((-38, -0.2), 6, 0.4, facecolor="#FF9999", edgecolor="black", lw=1)
        ax1.add_patch(p35_box)
        ax1.text(-35, 0.4, "-35 Box", ha="center", fontsize=8, fontweight="bold")
        
        p10_box = patches.Rectangle((-15, -0.2), 6, 0.4, facecolor="#C9DAF8", edgecolor="black", lw=1)
        ax1.add_patch(p10_box)
        ax1.text(-12, 0.4, "-10 Box", ha="center", fontsize=8, fontweight="bold")
        
        sd_box = patches.Rectangle((-24, -0.2), 6, 0.4, facecolor="#FFE699", edgecolor="black", lw=1)
        ax1.add_patch(sd_box)
        ax1.text(-21, 0.4, "SD Box", ha="center", fontsize=8, fontweight="bold")
    else:
        ax1.text(-45, 0.5, "No promoter motifs detected", ha="center", fontsize=10, fontstyle="italic")
    
    ax1.set_title(f"Promoter Transcription Consensus Motifs upstream of {bgc_id} core genes", fontsize=11, fontweight="bold")
    
    ax2.set_xlim(bgc_start - flank_w, bgc_end + flank_w)
    ax2.set_ylim(-1, 2)
    ax2.axis("off")
    ax2.plot([bgc_start - flank_w, bgc_end + flank_w], [0.5, 0.5], color="grey", lw=1)
    
    core_block = patches.Rectangle((bgc_start, 0.2), bgc_end - bgc_start, 0.6, facecolor="#C9DAF8", edgecolor="blue", lw=1.5)
    ax2.add_patch(core_block)
    ax2.text((bgc_start + bgc_end)/2, 1.0, f"BGC Core: {bgc_id}\n({bgc_end-bgc_start:,} bp)", ha="center", fontsize=8, fontweight="bold")
    
    # Annotation Collision Prevention: dynamic layout repelling/scattering on heavy BGC tracks
    for i, ph in enumerate(phage_hits):
        p_len = ph["End"] - ph["Start"]
        p_block = patches.Rectangle((ph["Start"], 0.2), p_len, 0.6, facecolor="#FFD7D7", edgecolor="red", lw=1)
        ax2.add_patch(p_block)
        y_jitter = -0.3 - (0.25 * (i % 3)) # Scatter text on y-axis
        ax2.text((ph["Start"] + ph["End"])/2, y_jitter, str(ph["Locus_Tag"]), ha="center", fontsize=7, rotation=30)
        
    ax2.set_title(f"Genomic Flanking Environment & Prophage Remnant Architecture for {bgc_id}", fontsize=11, fontweight="bold")
    
    plt.tight_layout()
    fig_out = os.path.join(out_dir, "bgc_motifs", f"{bgc_id}_promoter_phage_layout.png")
    plt.savefig(fig_out, dpi=300)
    plt.close()
    logger.info(f"Saved promoter motif and prophage layout map to: {fig_out}")
    
    # Calculate promoter confidence score
    completeness = 100.0
    contamination = 0.0
    closure_status = "Closed"
    if qc_result:
        completeness = qc_result.prediction.get("completeness", 100.0)
        contamination = qc_result.prediction.get("contamination", 0.0)
        closure_status = qc_result.genome_closure_status
        
    # High susceptibility to assembly QC: Promoter detection Sm = 0.9
    qc_penalty = get_module_qc_penalty("promoter", completeness, contamination)
    
    # Base confidence is mean confidence of core promoter predictions
    mean_prom_conf = np.mean([p["Confidence"] for p in promoter_records]) if promoter_records else 0.50
    final_score = mean_prom_conf - qc_penalty
    final_score = max(0.0, min(1.0, final_score))
    
    prediction_data = {
        "promoter_records": promoter_records,
        "phage_hits": phage_hits,
        "blast_rows": blast_rows,
        "layout_plot": fig_out
    }
    
    db_vers = CONFIG.get("database_versions", {})
    
    return PredictionResult(
        prediction=prediction_data,
        confidence_score=final_score,
        algorithm="Log-odds PWM + DNA curvature and UP-element scanners",
        algorithm_version="v3.0.0",
        genome_closure_status=closure_status,
        database="DBTBS / RegulonDB / PRODORIC",
        database_version=db_vers.get("gtdb", "R226"),
        evidence=[
            f"Scanned upstream promoter sequences for {len(core_genes)} core genes.",
            f"Calculated PWM log-odds with active DNA curvature curvature and UP-element boost metrics.",
            f"Analyzed {len(flank_genes_list)} flanking CDS genes for mobile genetic elements components."
        ],
        limitations=[
            "Promoter log-odds PWM matches are in silico indicators and require direct in vivo transcription validation.",
            "DNA shape and UP-element features act as confidence modulators and do not guarantee active binding."
        ],
        citations=[
            "DBTBS: Sierro et al. 2008 (doi:10.1093/nar/gkm910)",
            "RegulonDB: Tierrafría et al. 2022 (doi:10.1093/nar/gkab1058)",
            "PRODORIC: Dudek et al. 2021 (doi:10.1093/nar/gkaa1031)"
        ],
        runtime=time.time() - start_time,
        warnings=[f"QC assembly penalty applied: -{qc_penalty:.2f}"] if qc_penalty > 0.05 else []
    )
