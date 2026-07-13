import os
import time
import platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from Bio import SeqIO
from redolarium.utils import save_publication_plot
from redolarium.structures import PredictionResult
from redolarium.config import CONFIG
from redolarium.qc import get_module_qc_penalty

def run_phylogeny_pipeline(query_gb, ref_strains, identities, sim_matrix, out_dir, logger, bgc_blast_results=None, qc_result=None):
    logger.info("Stage 7: Conducting Molecular Clocking & Phylogenetics...")
    start_time = time.time()
    
    # Extract query organism
    query_record = max(list(SeqIO.parse(query_gb, "genbank")), key=lambda r: len(r.seq))
    query_org = query_record.annotations.get("organism", "Query Species")
    query_genus = query_org.split()[0] if len(query_org.split()) > 0 else "Bacterial"
    
    # ─── 1. BGC Presence/Absence Profile mapping ───
    present_accessions = set()
    if bgc_blast_results:
        for hit in bgc_blast_results:
            acc = hit.get("Accession", "")
            if acc:
                present_accessions.add(acc.split(".")[0].lower()) # Strip version numbers
                
    presence_profile = {}
    presence_profile[query_org] = "[+]" # Query self is always present
    for strain in ref_strains:
        name = strain.split(" (")[0]
        acc = strain.split(" (")[1].replace(")", "") if " (" in strain else strain
        acc_clean = acc.split(".")[0].lower()
        if acc_clean in present_accessions:
            presence_profile[name] = "[+]"
        else:
            presence_profile[name] = "[-]"
            
    # Save Presence/Absence Matrix
    pa_rows = []
    pa_rows.append({"Genome_Assembly": query_org, "Status": "Present", "Code": "[+]"})
    for strain in ref_strains:
        name = strain.split(" (")[0]
        code = presence_profile.get(name, "[-]")
        status = "Present" if code == "[+]" else "Absent"
        pa_rows.append({"Genome_Assembly": name, "Status": status, "Code": code})
        
    df_pa = pd.DataFrame(pa_rows)
    csv_pa = os.path.join(out_dir, "tabular_data", "bgc_presence_absence.csv")
    df_pa.to_csv(csv_pa, index=False)
    logger.info(f"Saved BGC presence/absence matrix to: {csv_pa}")

    # ─── 2. Calculate Divergence Speciation Times ───
    div_cfg = CONFIG.get("evolution", {}).get("divergence_model", {})
    lineage_rates = div_cfg.get("lineage_rate", {})
    
    rate_key = "default"
    query_org_lower = query_org.lower()
    if "symbiont" in query_org_lower or "buchnera" in query_org_lower or "wolbachia" in query_org_lower:
        rate_key = "endosymbiont"
    elif "intracellular" in query_org_lower or "chlamydia" in query_org_lower or "mycoplasma" in query_org_lower:
        rate_key = "intracellular"
    else:
        rate_key = "default"
        
    rate = lineage_rates.get(rate_key, 1.0)
    base_clock_multiplier = CONFIG.get("speciation_clock_multiplier", 15.0)
    clock_multiplier = base_clock_multiplier * rate
    
    # Peer-review defensible clock disclaimer
    clock_disclaimer = (
        "Divergence estimates are heuristic approximations and should not be interpreted "
        "as absolute evolutionary times because bacterial genomes frequently undergo "
        "horizontal gene transfer and lineage-specific rate variation."
    )
    
    divergence_events = []
    
    # Add root node (ancestral root based on lowest similarity reference)
    valid_idents = [idnt for idnt in identities if idnt > 0.0]
    min_ident = min(valid_idents) if valid_idents else 70.0
    root_time = round((100.0 - min_ident) * clock_multiplier, 2)
    
    divergence_events.append({
        "Event_ID": "EV_ROOT",
        "Node_Reference": f"Ancestral {query_genus} progenitor",
        "Estimated_Time_Mya (approx)": root_time,
        "Event_Type": "Ancestral Root",
        "Organism_Where": f"Ancestral {query_genus} niche",
        "Description": f"Estimated divergence root of ancestral {query_genus} lineages",
        "Calibration_Note": clock_disclaimer,
        "Clock_Disclaimer": clock_disclaimer
    })
    
    # Map top 5 reference strains to speciation nodes based on actual sequence identities
    added_nodes = 0
    for idx, (strain, ident) in enumerate(zip(ref_strains[:5], identities[:5])):
        if ident <= 0.0:
            continue
        name = strain.split(" (")[0]
        acc = strain.split(" (")[1].replace(")", "") if " (" in strain else strain
        div_time = round((100.0 - ident) * clock_multiplier, 2)
        
        divergence_events.append({
            "Event_ID": f"EV_{added_nodes+1:03d}",
            "Node_Reference": f"{name} split node",
            "Estimated_Time_Mya (approx)": div_time,
            "Event_Type": "Speciation",
            "Organism_Where": "Ecological diversification boundary",
            "Description": f"Divergence time of {query_org} progenitor from {name} ({acc}) clade",
            "Calibration_Note": clock_disclaimer,
            "Clock_Disclaimer": clock_disclaimer
        })
        added_nodes += 1
        
    divergence_events.append({
        "Event_ID": "EV_QUERY",
        "Node_Reference": f"{query_org} (Query)",
        "Estimated_Time_Mya (approx)": 0.0,
        "Event_Type": "Extant",
        "Organism_Where": "Query genome assembly",
        "Description": f"Extant isolate sequence of {query_org}",
        "Calibration_Note": "Present-day reference (0 Mya by definition)",
        "Clock_Disclaimer": "Present-day reference (0 Mya by definition)"
    })
    
    df_phy = pd.DataFrame(divergence_events)
    csv_out = os.path.join(out_dir, "tabular_data", "divergence_events.csv")
    df_phy.to_csv(csv_out, index=False)
    logger.info(f"Saved real divergence timeline dataset to: {csv_out}")

    # ─── 3. Visual: UPGMA Phylogenetic Placement Dendrogram ───
    tree_fig_out = os.path.join(out_dir, "phylogeny_trees", "phylogenetic_tree.png")
    try:
        dist_matrix = 100.0 - sim_matrix
        np.fill_diagonal(dist_matrix, 0.0)
        dist_matrix = (dist_matrix + dist_matrix.T) / 2.0
        condensed_dist = squareform(dist_matrix)
        linkage_matrix = sch.linkage(condensed_dist, method='average')
        
        tree_labels = []
        tree_labels.append(f"★ QUERY: {query_org[:22]} [+]")
        for strain in ref_strains:
            name = strain.split(" (")[0]
            status_code = presence_profile.get(name, "[-]")
            tree_labels.append(f"{name[:22]} {status_code}")
            
        plt.figure(figsize=(10, 8))
        sch.dendrogram(
            linkage_matrix,
            labels=tree_labels,
            orientation='left',
            leaf_font_size=8,
            color_threshold=0
        )
        plt.title("Phylogenetic Dendrogram (UPGMA based on Genomic AAI Similarity)\n[+] / [-] represents BGC presence/absence", fontsize=11, fontweight="bold", pad=15)
        plt.xlabel("Genomic Distance (100 - AAI %)", fontsize=10, fontweight="bold")
        plt.tight_layout()
        save_publication_plot(tree_fig_out, dpi=300)
        plt.close()
        logger.info(f"Saved UPGMA phylogenetic tree placement to: {tree_fig_out}")
    except Exception as ex:
        logger.warning(f"Failed to generate UPGMA phylogenetic dendrogram: {ex}")

    # ─── 4. Visual: Divergence Speciation Timeline Flowchart ───
    fig_out = os.path.join(out_dir, "phylogeny_trees", "divergence_timeline.png")
    try:
        G = nx.DiGraph()
        nodes_pos = {}
        
        n_root = f"Root\n({root_time} Mya)"
        nodes_pos[n_root] = (-root_time, 1.0)
        G.add_node(n_root)
        
        prev_node = n_root
        for idx, event in enumerate(divergence_events[1:-1]):
            div_time = event["Estimated_Time_Mya (approx)"]
            label = event["Node_Reference"].split(" split")[0]
            n_name = f"{label[:15]}\n({div_time} Mya)"
            nodes_pos[n_name] = (-div_time, 1.0 - (idx * 0.3))
            G.add_node(n_name)
            G.add_edge(prev_node, n_name)
            prev_node = n_name
            
        n_query = f"{query_org[:15]}\n(0 Mya)"
        nodes_pos[n_query] = (0.0, 0.5)
        G.add_node(n_query)
        G.add_edge(prev_node, n_query)
        
        plt.figure(figsize=(10, 6))
        nx.draw_networkx_nodes(G, nodes_pos, node_size=1500, node_color="#FFE699", edgecolors="black", linewidths=1)
        nx.draw_networkx_labels(G, nodes_pos, font_size=8, font_weight="bold")
        nx.draw_networkx_edges(G, nodes_pos, width=1.5, edge_color="black", arrowstyle="->", arrowsize=12)
        
        plt.title(f"Molecular Clocking & Speciation Timeline: {query_org}", fontsize=11, fontweight="bold", pad=15)
        plt.xlabel("Evolutionary Divergence Time (Millions of Years Ago - Mya)", fontsize=10, fontweight="bold")
        plt.axvline(0, color="grey", linestyle=":")
        plt.xlim(-root_time - 15, 10)
        plt.gca().get_yaxis().set_visible(False)
        plt.tight_layout()
        save_publication_plot(fig_out, dpi=300)
        plt.close()
        logger.info(f"Saved speciation timeline flowchart to: {fig_out}")
    except Exception as ex:
        logger.warning(f"Failed to generate Speciation Timeline Flowchart: {ex}")
        
    # Calculate phylogeny module confidence score
    completeness = 100.0
    contamination = 0.0
    closure_status = "Closed"
    if qc_result:
        completeness = qc_result.prediction.get("completeness", 100.0)
        contamination = qc_result.prediction.get("contamination", 0.0)
        closure_status = qc_result.genome_closure_status
        
    # Penalty calculation
    qc_penalty = get_module_qc_penalty("hgt", completeness, contamination)
    phy_score = 0.95 - qc_penalty
    phy_score = max(0.0, min(1.0, phy_score))
    
    prediction_data = {
        "divergence_events": divergence_events,
        "tree_image": tree_fig_out,
        "timeline_image": fig_out
    }
    
    db_vers = CONFIG.get("database_versions", {})
    
    return PredictionResult(
        prediction=prediction_data,
        confidence_score=phy_score,
        algorithm="Evolutionary Divergence Heuristic",
        algorithm_version="v3.0.0",
        genome_closure_status=closure_status,
        database="NCBI Taxonomy / TreeOfLife",
        database_version=db_vers.get("gtdb", "R226"),
        evidence=[
            f"Lineage divergence multiplier adjusted for rate factor: {rate:.2f}x",
            f"Calculated root node divergence time: {root_time} Mya",
            f"Speciation timeline built with {len(divergence_events)} nodes"
        ],
        limitations=[
            "Divergence estimates are heuristic approximations and should not be interpreted as absolute evolutionary times.",
            "Bacterial rate variation and recombination events violate standard molecular clock assumptions."
        ],
        citations=[
            "Ochman et al. 2000 (doi:10.1038/35012500)",
            "Ho & Duchene 2014 (doi:10.1093/molbev/msu221)"
        ],
        runtime=time.time() - start_time,
        metadata={"parameters": {"lineage_rate": rate, "clock_multiplier": clock_multiplier}}
    )
