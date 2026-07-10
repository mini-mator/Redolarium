# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Redolarium Contributors, IIT Kanpur
# See LICENSE and THIRD_PARTY_LICENSES.md for full licence details.

import os
import re
import csv
import numpy as np
import pandas as pd
from Bio import SeqIO
from redolarium.config import CONFIG

# ==========================================
# CUE MAPPING: Internal Dictionary
# ==========================================
# Robust, curated static rule-based internal dictionary targeting specific Sigma factors
# and primary transcription factors, avoiding external live DB API calls.

CUE_MAPPING = {
    # Sigma Factors
    "sigma-b": "General Stress / Environmental Stress",
    "sigma-h": "Heat Shock / Oxidative Stress",
    "sigma-e": "Extracytoplasmic / Envelope Stress",
    "sigma-w": "Alkaline Shock / Toxin Resistance",
    "sigma-m": "Cell Wall Stress",
    "sigma-f": "Flagellar Synthesis / Chemotaxis",
    "sigma-n": "Nitrogen Limitation (sigma-54)",
    "sigma-s": "Stationary Phase / Starvation (RpoS)",
    "sigma-70": "Constitutive (Housekeeping / Exponential growth)",
    "sigma-a": "Constitutive (Housekeeping / Exponential growth)",
    
    # Transcription Factors
    "fur": "Iron Limitation",
    "oxyr": "Oxidative Stress (Peroxides)",
    "soxr": "Oxidative Stress (Superoxide)",
    "lexa": "DNA Damage (SOS Response)",
    "fnr": "Anaerobic / Hypoxia",
    "arca": "Anaerobic / Respiratory Shift",
    "phob": "Phosphate Limitation",
    "cpxr": "Envelope Stress / Protein Misfolding",
    "ompr": "Osmotic Stress",
    "zur": "Zinc Limitation"
}

def extract_promoter_cues(promoter_features):
    """
    Parses a list of identified promoter/regulatory features and maps them to environmental cues.
    Returns a list of mapped cues with evidence.
    """
    mapped_cues = []
    
    for feature in promoter_features:
        name = feature.get("Name", "").lower()
        desc = feature.get("Description", "").lower()
        
        # Search dictionary for hits
        for key, cue in CUE_MAPPING.items():
            if key in name or key in desc:
                mapped_cues.append({
                    "Regulatory_Element": key.upper(),
                    "Feature_Name": feature.get("Name", "Unknown"),
                    "Environmental_Cue": cue,
                    "Location": feature.get("Location", "Unknown"),
                    "Score": feature.get("Score", 0.0)
                })
                
    return mapped_cues

def scan_upstream_promoters(bgc_start, bgc_end, strand, sequence, window=500):
    """
    Scans the immediate upstream region of the entire BGC operon cluster.
    """
    seq_len = len(sequence)
    promoters = []
    
    # Define upstream region based on operon strand
    if strand == "+":
        up_start = max(0, bgc_start - window)
        up_end = bgc_start
        up_seq = sequence[up_start:up_end]
    else:
        up_start = bgc_end
        up_end = min(seq_len, bgc_end + window)
        from Bio.Seq import Seq
        up_seq = str(Seq(sequence[up_start:up_end]).reverse_complement())
        
    # Basic motif scanning (heuristic)
    # -10 and -35 boxes for generic sigma-70 / sigma-A
    # This is a simplified regex matching for standard motifs.
    # Reference: Hawley & McClure 1983 (sigma70)
    
    motifs = {
        "sigma-70": [r"TTGACA.{15,19}TATAAT", "Constitutive"],
        "sigma-b": [r"GTTT.{13,17}GGGTAT", "Stress Response"],
        "sigma-n": [r"TGGCAC.{4,6}TTGCW", "Nitrogen Limitation"]
    }
    
    for sigma, (pattern, default_cue) in motifs.items():
        for match in re.finditer(pattern, up_seq, re.IGNORECASE):
            loc_offset = up_start + match.start() if strand == "+" else up_end - match.end()
            promoters.append({
                "Name": sigma,
                "Description": f"Putative {sigma} promoter motif",
                "Location": loc_offset,
                "Score": 0.85
            })
            
    return promoters

def run_promoter_prediction(query_gb, target_bgc, out_dir, logger):
    """
    Genome-Aware Promoter Identification and Environmental Cue Mapping.
    """
    logger.info("Stage 10: Genome-Aware Promoter Identification and Environmental Cue Mapping...")
    
    bgc_id = target_bgc["BGC_ID"]
    bgc_start = target_bgc["Core_Start"]
    bgc_end = target_bgc["Core_End"]
    strand = "+" # Assume forward for the entire cluster or take from major operon
    
    try:
        record = SeqIO.read(query_gb, "genbank")
    except Exception:
        try:
            record = SeqIO.read(query_gb, "fasta")
        except Exception as e:
            logger.warning(f"Failed to read genome for promoter prediction: {e}")
            return []
            
    sequence = str(record.seq)
    
    # 1. Scan immediate upstream region
    promoters = scan_upstream_promoters(bgc_start, bgc_end, strand, sequence)
    
    if not promoters:
        logger.info(f"No strong promoter motifs detected upstream of BGC {bgc_id}.")
        
    # 2. Map to Environmental Cues
    cues = extract_promoter_cues(promoters)
    
    if cues:
        os.makedirs(os.path.join(out_dir, "regulatory_networks"), exist_ok=True)
        csv_out = os.path.join(out_dir, "regulatory_networks", f"{bgc_id}_environmental_cues.csv")
        df_cues = pd.DataFrame(cues)
        df_cues.to_csv(csv_out, index=False)
        logger.info(f"Saved Environmental Cue mappings to: {csv_out}")
    
    return cues
