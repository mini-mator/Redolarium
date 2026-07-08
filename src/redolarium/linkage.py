import os
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from collections import Counter
from Bio import SeqIO

# Standard host amino acid synthesis details (empirical ATP costs from Akashi & Gojobori 2002)
AA_METABOLIC_PATHWAY = {
    'A': {'name': 'Alanine', 'pathway': 'map00250: Alanine/glutamate metabolism', 'genes': 'alaS, alaT', 'atp_base': 11.7},
    'R': {'name': 'Arginine', 'pathway': 'map00220: Arginine biosynthesis', 'genes': 'argS, argA, argB, argC', 'atp_base': 27.3},
    'N': {'name': 'Asparagine', 'pathway': 'map00250: Alanine/glutamate metabolism', 'genes': 'asnS, ansA, asnO', 'atp_base': 14.7},
    'D': {'name': 'Aspartate', 'pathway': 'map00250: Alanine/glutamate metabolism', 'genes': 'aspS, aspB, pyc', 'atp_base': 12.7},
    'C': {'name': 'Cysteine', 'pathway': 'map00270: Cysteine/methionine metabolism', 'genes': 'cysS, cysK, cysE', 'atp_base': 24.7},
    'E': {'name': 'Glutamate', 'pathway': 'map00250: Alanine/glutamate metabolism', 'genes': 'gltX, gltA, gltB', 'atp_base': 15.3},
    'Q': {'name': 'Glutamine', 'pathway': 'map00250: Alanine/glutamate metabolism', 'genes': 'glnS, glnA', 'atp_base': 16.3},
    'G': {'name': 'Glycine', 'pathway': 'map00260: Glycine/serine/threonine metabolism', 'genes': 'glyS, glyA', 'atp_base': 11.7},
    'H': {'name': 'Histidine', 'pathway': 'map00340: Histidine metabolism', 'genes': 'hisS, hisA, hisB', 'atp_base': 38.3},
    'I': {'name': 'Isoleucine', 'pathway': 'map00290: Valine/leucine/isoleucine biosynthesis', 'genes': 'ileS, ilvB, ilvC', 'atp_base': 32.3},
    'L': {'name': 'Leucine', 'pathway': 'map00290: Valine/leucine/isoleucine biosynthesis', 'genes': 'leuS, leuA, leuB', 'atp_base': 27.3},
    'K': {'name': 'Lysine', 'pathway': 'map00300: Lysine biosynthesis', 'genes': 'lysS, lysA, dapA', 'atp_base': 30.3},
    'M': {'name': 'Methionine', 'pathway': 'map00270: Cysteine/methionine metabolism', 'genes': 'metS, metE, metC', 'atp_base': 34.3},
    'F': {'name': 'Phenylalanine', 'pathway': 'map00400: Phenylalanine/tyrosine/tryptophan biosynthesis', 'genes': 'pheS, pheA, aroG', 'atp_base': 52.0},
    'P': {'name': 'Proline', 'pathway': 'map00330: Arginine/proline metabolism', 'genes': 'proS, proA, proB', 'atp_base': 20.3},
    'S': {'name': 'Serine', 'pathway': 'map00260: Glycine/serine/threonine metabolism', 'genes': 'serS, serA, serB', 'atp_base': 11.7},
    'T': {'name': 'Threonine', 'pathway': 'map00260: Glycine/serine/threonine metabolism', 'genes': 'thrS, thrA, thrB', 'atp_base': 18.7},
    'W': {'name': 'Tryptophan', 'pathway': 'map00400: Phenylalanine/tyrosine/tryptophan biosynthesis', 'genes': 'trpS, trpA, trpB', 'atp_base': 74.3},
    'Y': {'name': 'Tyrosine', 'pathway': 'map00400: Phenylalanine/tyrosine/tryptophan biosynthesis', 'genes': 'tyrS, tyrA', 'atp_base': 50.0},
    'V': {'name': 'Valine', 'pathway': 'map00290: Valine/leucine/isoleucine biosynthesis', 'genes': 'valS, ilvB, ilvC', 'atp_base': 23.3}
}

def get_bioactive_peptide_product(desc):
    desc_l = desc.lower()
    if "surfactin" in desc_l:
        return "ELLVDLL", "Surfactin Heptapeptide Product"
    elif "fengycin" in desc_l:
        return "EYLDVTTEFV", "Fengycin Decapeptide Product"
    elif "plipastatin" in desc_l:
        return "EYLDVTTEFV", "Plipastatin Decapeptide Product"
    elif "iturin" in desc_l:
        return "NDNYEPV", "Iturin Heptapeptide Product"
    elif "nisin" in desc_l:
        return "ITSISLCTPGCKGALMGGCNMKTATCHCSIHVSK", "Nisin Lanthipeptide Product"
    elif "subtilin" in desc_l:
        return "WKSESLCTPGCKGALLGGCNMKTATCHCSIHVSK", "Subtilin Lanthipeptide Product"
    elif "subtilosin" in desc_l:
        return "NKGCATCSIGAACLVDGPIPDFIAGIMGGG", "Subtilosin Sactipeptide Product"
    elif "bacillibactin" in desc_l:
        return "DGT", "Bacillibactin Tripeptide Product"
    elif "bacilysin" in desc_l:
        return "AD", "Bacilysin Dipeptide Product"
    return None, None

def extract_bgc_product_peptide(region_genes, target_bgc=None):
    # 1. Scan region genes for structural translations (RiPPs)
    for g in region_genes:
        seq = g.get("Protein_Sequence", "")
        if seq and 20 < len(seq) < 150:
            # Found suitable precursor peptide translation
            return seq, "precursor_peptide"
            
    # 2. Try mapping compound name from target_bgc description/metadata (NRPS / others)
    if target_bgc:
        desc = target_bgc.get("BGC_Type", "")
        desc_str = str(desc)
        seq, label = get_bioactive_peptide_product(desc_str)
        if seq:
            return seq, "bioactive_peptide"
    
    # 3. If none found, concatenate core BGC enzyme translations as stoichiometric target sequence
    concat_seq = ""
    for g in region_genes:
        seq = g.get("Protein_Sequence", "")
        if seq:
            concat_seq += seq
    if concat_seq:
        return concat_seq, "biosynthetic_enzymes"
        
    # Return (None, None) instead of a fake sequence
    return None, None

def calculate_housekeeping_baseline(record, logger, out_dir, bgc_id):
    logger.info("Calculating thermodynamic baselines for core housekeeping proteins...")
    
    housekeeping_markers = ["recA", "gyrB", "rpoB", "gyrA", "dnaN", "pyrH", "dnaA", "pheS", "infB", "rpoA"]
    
    hk_rows = []
    for feat in record.features:
        if feat.type == "CDS":
            g_sym = feat.qualifiers.get("gene", [""])[0]
            prod = feat.qualifiers.get("product", [""])[0].lower()
            trans = feat.qualifiers.get("translation", [""])[0]
            ltag = feat.qualifiers.get("locus_tag", [""])[0]
            
            is_hk = any(m == g_sym for m in housekeeping_markers) or any(m in prod for m in ["ribosomal protein s", "ribosomal protein l", "rna polymerase subunit beta", "dna gyrase subunit"])
            
            if is_hk and trans:
                atp_cost = 0.0
                for aa in trans:
                    if aa in AA_METABOLIC_PATHWAY:
                        atp_cost += AA_METABOLIC_PATHWAY[aa]["atp_base"]
                    else:
                        atp_cost += 4.0
                processing_cost = len(trans) * 4.0
                total_cost = atp_cost + processing_cost
                
                hk_rows.append({
                    "Locus_Tag": ltag,
                    "Gene_Symbol": g_sym if g_sym else "NA",
                    "Product": feat.qualifiers.get("product", ["Predicted housekeeping protein"])[0],
                    "Length_aa": len(trans),
                    "Precursor_ATP_Cost": round(atp_cost, 2),
                    "Translation_ATP_Cost": round(processing_cost, 2),
                    "Total_ATP_Cost": round(total_cost, 2)
                })
                
    df_hk = pd.DataFrame(hk_rows)
    csv_hk = os.path.join(out_dir, "tabular_data", f"{bgc_id}_housekeeping_atp_stoichiometry.csv")
    df_hk.to_csv(csv_hk, index=False)
    logger.info(f"Saved housekeeping thermodynamic baseline stoichiometry to: {csv_hk}")
    return hk_rows

def run_linkage_pipeline(query_gb, region_genes, target_bgc, out_dir, logger):
    logger.info("Stage 9: Mapping Metabolic-BGC Linkage and Stoichiometry...")
    
    bgc_id = target_bgc["BGC_ID"]
    peptide_seq, target_type = extract_bgc_product_peptide(region_genes)
    
    if not peptide_seq:
        logger.warning(f"No suitable precursor peptide sequence or BGC enzymes found for linkage stoichiometry. Skipping linkage pipeline for BGC {bgc_id}.")
        return []
        
    # Calculate thermodynamic baseline for housekeeping genes
    record = SeqIO.read(query_gb, "genbank")
    calculate_housekeeping_baseline(record, logger, out_dir, bgc_id)
    
    # Calculate amino acid counts of the BGC peptide
    aa_counts = Counter(peptide_seq)
    logger.info(f"Prepropeptide parsed (length: {len(peptide_seq)} aa, type: {target_type}). Mapped AA demand: {dict(aa_counts)}")
    
    bgc_type = target_bgc.get("BGC_Type", "Other")
    is_nrps = "nrps" in bgc_type.lower()
    
    stoichiometry_rows = []
    total_precursor_atp = 0
    total_trans_atp = 0
    total_mod_atp = 0
    
    for aa, count in aa_counts.items():
        if aa in AA_METABOLIC_PATHWAY:
            info = AA_METABOLIC_PATHWAY[aa]
            precursor_atp = count * info["atp_base"]
            
            # Nonribosomal peptide synthesis (NRPS) bypasses ribosomal translation translation cost (4.0 ATP/aa elongation)
            translation_atp = 0.0 if is_nrps else count * 4.0
            
            # Dehydration post-translational modifications (consuming 2 ATP per Ser/Thr)
            mod_atp = 0
            if target_type == "precursor_peptide" and aa in ['S', 'T']:
                mod_atp = count * 2
                
            total_precursor_atp += precursor_atp
            total_trans_atp += translation_atp
            total_mod_atp += mod_atp
            
            stoichiometry_rows.append({
                "Amino_Acid": info["name"],
                "Code": aa,
                "Precursor_Count": count,
                "Host_Pathway": info["pathway"],
                "Candidate_Synthetases": info["genes"],
                "Precursor_ATP_Cost": precursor_atp,
                "Translation_ATP_Cost": translation_atp,
                "Modification_ATP_Cost": mod_atp,
                "Total_ATP_Cost": precursor_atp + translation_atp + mod_atp
            })
            
    df_st = pd.DataFrame(stoichiometry_rows)
    
    summary_row = {
        "Amino_Acid": "TOTAL STOICHIOMETRIC DEMAND",
        "Code": "-",
        "Precursor_Count": len(peptide_seq),
        "Host_Pathway": "All pathways active",
        "Candidate_Synthetases": "-",
        "Precursor_ATP_Cost": total_precursor_atp,
        "Translation_ATP_Cost": total_trans_atp,
        "Modification_ATP_Cost": total_mod_atp,
        "Total_ATP_Cost": total_precursor_atp + total_trans_atp + total_mod_atp
    }
    df_st = pd.concat([df_st, pd.DataFrame([summary_row])], ignore_index=True)
    
    csv_out = os.path.join(out_dir, "tabular_data", f"{bgc_id}_atp_stoichiometry.csv")
    df_st.to_csv(csv_out, index=False)
    logger.info(f"Saved ATP stoichiometry details to: {csv_out}")

    bgc_type = target_bgc.get("BGC_Type", "Other")
    if "Lantipeptide" in bgc_type or "RiPP" in bgc_type:
        synthesis_node = f"Ribosomal\nTranslation ({bgc_id})"
        mod_node = "Post-Translational\nModification (LanB/C)"
        export_node = "ABC Transporter\nExport (LanT)"
    elif "NRPS" in bgc_type:
        synthesis_node = f"Nonribosomal\nSynthesis ({bgc_id})"
        mod_node = "Backbone Assembly\n& Cyclization (TE)"
        export_node = "MFS/ABC Exporter\n(Efflux Pump)"
    elif "PKS" in bgc_type:
        synthesis_node = f"Polyketide\nSynthesis ({bgc_id})"
        mod_node = "Chain Assembly\n& Tailoring (TE)"
        export_node = "MFS/ABC Exporter\n(Efflux Pump)"
    else:
        synthesis_node = f"Core Scaffold\nSynthesis ({bgc_id})"
        mod_node = "Chemical Scaffold\nDecoration/Tailoring"
        export_node = "Active Transporter\nExport"

    # ==========================================
    # VISUAL: Systems Biology Pathway Network Map
    # ==========================================
    G = nx.DiGraph()
    
    pos = {
        "Glucose": (0, 4),
        "Glycolysis & TCA\n(Energy Pool)": (-2, 3),
        "ATP Generation\n(Substrate Phosphorylation)": (-3, 2),
        "Precursor\nCharged tRNAs": (3, 2),
        "tRNA Charging\nSynthetases": (0, 1.5),
        synthesis_node: (1.5, 0.5),
        mod_node: (-1.5, 0.5),
        export_node: (0, -0.5),
        "Extracellular Active\nBGC Product": (0, -1.5)
    }
    
    active_aas = [row for row in stoichiometry_rows if row.get("Precursor_Count", 0) > 0]
    num_aas = len(active_aas)
    for idx, aa_row in enumerate(active_aas):
        aa_name = aa_row["Amino_Acid"]
        x_pos = 1.0 + (idx / float(max(1, num_aas - 1))) * 2.0 if num_aas > 1 else 2.0
        y_pos = 3.2 + (0.3 if idx % 2 == 0 else -0.3)
        pos[aa_name] = (x_pos, y_pos)
        
    for node in pos:
        G.add_node(node)
        
    edges = [
        ("Glucose", "Glycolysis & TCA\n(Energy Pool)"),
        ("Glycolysis & TCA\n(Energy Pool)", "ATP Generation\n(Substrate Phosphorylation)"),
        ("tRNA Charging\nSynthetases", "Precursor\nCharged tRNAs"),
        ("ATP Generation\n(Substrate Phosphorylation)", "tRNA Charging\nSynthetases"),
        ("Precursor\nCharged tRNAs", synthesis_node),
        ("ATP Generation\n(Substrate Phosphorylation)", synthesis_node),
        (synthesis_node, mod_node),
        ("ATP Generation\n(Substrate Phosphorylation)", mod_node),
        (mod_node, export_node),
        ("ATP Generation\n(Substrate Phosphorylation)", export_node),
        (export_node, "Extracellular Active\nBGC Product")
    ]
    
    for aa_row in active_aas:
        aa_name = aa_row["Amino_Acid"]
        edges.append(("Glucose", aa_name))
        edges.append(("Glycolysis & TCA\n(Energy Pool)", aa_name))
        edges.append((aa_name, "Precursor\nCharged tRNAs"))
        
    G.add_edges_from(edges)
    
    # ==========================================
    # VISUAL: Figure 5 (Thermodynamic Precursor Burden Chart)
    # ==========================================
    try:
        import seaborn as sns
        plt.figure(figsize=(10, 6))
        sns.set_theme(style="whitegrid")
        
        # Filter out summary row and sort by abundance
        df_plot = df_st[df_st["Amino_Acid"] != "TOTAL STOICHIOMETRIC DEMAND"].copy()
        df_plot = df_plot.sort_values(by="Precursor_Count", ascending=False)
        
        ax1 = plt.gca()
        # Bars for count
        sns.barplot(x="Code", y="Precursor_Count", data=df_plot, ax=ax1, palette="cividis", hue="Code", legend=False)
        ax1.set_xlabel("Amino Acid Residue", fontsize=10, fontweight="bold")
        ax1.set_ylabel("Abundance Count in Peptide", fontsize=10, fontweight="bold", color="blue")
        ax1.tick_params(axis='y', labelcolor="blue")
        
        # Dual axis line for total cost
        ax2 = ax1.twinx()
        ax2.plot(df_plot["Code"], df_plot["Total_ATP_Cost"], color="red", marker="o", linewidth=2, label="Stoichiometric Cost (ATP)")
        ax2.set_ylabel("Total Translation Cost (ATP equivalents)", fontsize=10, fontweight="bold", color="red")
        ax2.tick_params(axis='y', labelcolor="red")
        
        plt.title(f"Thermodynamic Amino Acid Precursor Burden: {bgc_id} Peptide", fontsize=12, fontweight="bold", pad=15)
        plt.tight_layout()
        
        burden_fig_out = os.path.join(out_dir, "metabolic_pathways", f"{bgc_id}_precursor_burden.png")
        os.makedirs(os.path.join(out_dir, "metabolic_pathways"), exist_ok=True)
        from redolarium.utils import save_publication_plot
        save_publication_plot(burden_fig_out, dpi=300)
        plt.close()
        logger.info(f"Saved precursor burden dual-axis plot to: {burden_fig_out}")
    except Exception as ex:
        logger.warning(f"Failed to generate precursor burden chart: {ex}")

    plt.figure(figsize=(10, 8))
    nx.draw_networkx_nodes(G, pos, node_size=1800, node_color="#C9DAF8", edgecolors="black", linewidths=1)
    nx.draw_networkx_labels(G, pos, font_size=8, font_weight="bold", font_family="sans-serif")
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=1.5, edge_color="grey", arrowstyle="-|>", arrowsize=12)
    
    plt.title(f"Metabolic Linkage Map: Host pathways supplying precursors and energy to {bgc_id}", fontsize=11, fontweight="bold", pad=15)
    plt.axis("off")
    plt.tight_layout()
    
    fig_out = os.path.join(out_dir, "metabolic_pathways", f"{bgc_id}_metabolic_linkage.png")
    from redolarium.utils import save_publication_plot
    save_publication_plot(fig_out, dpi=300)
    plt.close()
    logger.info(f"Saved metabolic linkage map to: {fig_out}")
    
    return df_st.to_dict('records')

