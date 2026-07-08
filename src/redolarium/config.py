# Central configuration and threshold constants for Redolarium Pipeline
CONFIG = {
    # Keywords for secondary metabolite BGC mining
    "bgc_keywords": [
        "lantibiotic", "lanthipeptide", "ripp", "mersacidin",
        "sublancin", "subtilin", "ericin", "paenibacillin",
        "lanthionine", "sactipeptide", "thiopeptide",
        "nisin", "epidermin", "polyketide", "nonribosomal",
        "nrps", "pks", "bacillaene", "fengycin", "plipastatin",
        "surfactin", "iturin", "mycosubtilin", "siderophore"
    ],

    # Flanking window sizes for BGC definition
    "bgc_flank_window": 10000,          # bp flanking coordinates around BGC core
    "promoter_search_upstream": 300,    # bp upstream of ATG for promoter motifs

    # Phage and mobile genetic elements keywords
    "phage_keywords": [
        "phage", "prophage", "integrase", "transposase", "recombinase",
        "capsid", "tail fiber", "terminase", "lysin", "holin",
        "excisionase", "att site", "insertion sequence", "is element",
        "mobile element", "conjugal transfer", "relaxase", "mobilization",
        "transposon", "tnp", "tns"
    ],

    # Bacterial regulatory consensus binding motifs with per-sigma spacer ranges
    # References: Hawley & McClure 1983 (σ70); Helmann 1995 (σA/σW/σB Bacillus);
    #             Boylan et al. 1992 (σB); Bose et al. 2008 (σ54/σN);
    #             Potvin et al. 2008 (σ54 Pseudomonas); Huang et al. 2017 (σS).
    # spacer_min/spacer_max: distance between end of –35 box and start of –10 box
    "sigma_motifs": {
        "sigmaA (housekeeping)": {
            "minus35": r"TTGAC[AT]",
            "minus10": r"TAT[AT]AT|TACTAT",
            "spacer_min": 16,
            "spacer_max": 18
        },
        "sigmaW (envelope stress)": {
            "minus35": r"TGGAA[ACT]",
            "minus10": r"CGT[AT]T",
            "spacer_min": 17,
            "spacer_max": 19
        },
        "sigmaB (general stress)": {
            "minus35": r"GTTTTAA",
            "minus10": r"GGGTAT",
            "spacer_min": 13,
            "spacer_max": 15
        },
        "sigmaS (stationary phase / RpoS)": {
            "minus35": r"TTGAC[AT]",
            "minus10": r"CTATAC[T]",
            "spacer_min": 16,
            "spacer_max": 18
        },
        "sigma54 (nitrogen/stress / RpoN)": {
            # sigma54 uses -12/-24 architecture, NOT -10/-35.
            # minus35 field here represents the -24 GG element;
            # minus10 field represents the -12 GC element.
            # Spacer is measured between the two conserved elements.
            "minus35": r"TGGCAC",
            "minus10": r"TTGC[AT]",
            "spacer_min": 10,
            "spacer_max": 12
        }
    },

    # Maximum inter-gene gap (bp) for grouping core BGC hits into a single cluster
    # Based on antiSMASH default: 10 kb between core biosynthetic genes
    # Reference: Blin et al. 2023 (antiSMASH 7.0, doi:10.1093/nar/gkad344)
    "bgc_clustering_gap_bp": 10000,

    # HGT sliding window thresholds (Langille & Bhatt 2006, Teeling et al. 2004)
    "hgt": {
        "gc_window": 5000,
        "gc_step": 500,
        "gc_zscore_threshold": 2.0,
        "tnf_window": 5000,
        "tnf_step": 2500
    },
    "speciation_clock_multiplier": 15.0,

    # Reference database citations with DOIs (complying with formatting requirements)
    "citations": [
        {"module": "Pan-Genomics", "tool": "NCBI Entrez Search", "ref": "Coordinators, N. R. (2018). Database resources of the National Center for Biotechnology Information.", "doi": "10.1093/nar/gky1111"},
        {"module": "Pan-Genomics", "tool": "Biopython Align", "ref": "Cock, P. J. et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics.", "doi": "10.1093/bioinformatics/btp163"},
        {"module": "BGC Annotation", "tool": "antiSMASH 7.0", "ref": "Blin, K. et al. (2023). antiSMASH 7.0: new and improved predictions for BGC detection.", "doi": "10.1093/nar/gkad344"},
        {"module": "HGT GC-bias", "tool": "GC-Profiler", "ref": "Langille, M. G. & Bhatt, I. (2006). Sliding window GC-profiling for horizontal transfers.", "doi": "10.1186/1471-2164-7-142"},
        {"module": "HGT TNF-bias", "tool": "Genomic Signature", "ref": "Teeling, H. et al. (2004). Application of tetranucleotide frequencies for fragment clustering.", "doi": "10.1111/j.1462-2920.2004.00638.x"},
        {"module": "HGT Scores", "tool": "HGT-Composite", "ref": "Ravenhall, M. et al. (2015). Inferring horizontal gene transfer using genomic signatures.", "doi": "10.1371/journal.pcbi.1004523"},
        {"module": "Phylogeny", "tool": "Bio.Phylo", "ref": "Talevich, N. et al. (2012). Bio.Phylo: phylogeny framework in Biopython.", "doi": "10.1186/1471-2105-13-209"},
        {"module": "Molecular Docking", "tool": "RCSB PDB API", "ref": "Berman, H. M. et al. (2000). The Protein Data Bank.", "doi": "10.1093/nar/28.1.235"},
        {"module": "Molecular Docking", "tool": "ChEMBL API", "ref": "Mendez, D. et al. (2019). ChEMBL: towards direct deposition of bioactivity data.", "doi": "10.1093/nar/gky1075"},
        {"module": "Metabolic Linkage", "tool": "KEGG REST API", "ref": "Kanehisa, M. et al. (2023). KEGG database for taxonomy-based pathways.", "doi": "10.1093/nar/gkac963"},
        {"module": "CAZyme Screening", "tool": "dbCAN3", "ref": "Huang, L. et al. (2023). dbCAN3: automated carbohydrate-active enzyme annotation.", "doi": "10.1093/nar/gkad343"},
        {"module": "AMR CARD", "tool": "CARD-RGI", "ref": "Alcock, B. P. et al. (2023). CARD 2023: Comprehensive Antibiotic Resistance Database.", "doi": "10.1093/nar/gkac920"},
        {"module": "Virulence VFDB", "tool": "VFDB", "ref": "Chen, L. et al. (2016). VFDB 2016: hierarchical dataset for virulence factors.", "doi": "10.1093/nar/gkv1239"}
    ],
    
    # Metabolism & KEGG REST configuration
    "metabolism": {
        "kegg_cache_dir": "resources/kegg_cache",
        "kegg_cache_expiry_days": 7,
        "kegg_timeout": 30,
        "kegg_proxy": None,
        "mode": "fast"  # Options: "fast" (KEGG REST reconstruction), "deep" (gapseq/MetaPathPredict)
    },

    # Screening & Biosafety configurations
    "screening": {
        "weights": {
            "identity": 0.3,
            "coverage": 0.3,
            "bitscore": 0.4
        },
        "thresholds": {
            "min_identity": 40.0,    # 40%
            "max_evalue": 1e-5
        },
        "hmm_db_path": "resources/essential_bgc.hmm",
        "manifest_path": "resources/manifest.json",
        "risk_matrix_path": "config/risk_matrix.json",
        "qs_signatures_path": "config/qs_signatures.json",
        "vf_signatures_path": "config/vf_signatures.json"
    }
}


# openpyxl stylesheet color palette
EXCEL_PALETTE = {
    "header_fill": "1F3864",     # Navy Blue
    "alt_row_fill": "E8F0FE",    # Light Blue-Grey
    "positive_fill": "E2F0CB",   # Present/Green
    "negative_fill": "FFE0E0",   # Absent/Red
    "warning_fill": "FF9999",    # Outlier/Flagged
    "critical_fill": "FFD7D7",   # Phage/MGE (light pink)
    "highlight_fill": "FFE699",  # Summary composite row (yellow)
    "bgc_core_fill": "C9DAF8"    # Core BGC genes (soft blue)
}

# Excel tab color associations
TAB_COLORS = {
    "Summary": "1F3864",
    "BGC_Gene_Architecture": "2E75B6",
    "Promoter_Analysis": "70AD47",
    "Phage_Artifacts": "FF0000",
    "HGT_Evidence": "ED7D31",
    "GC_Profile": "A9D18E",
    "Biosynthetic_Pathway": "4472C4",
    "BLAST_Homologs": "7030A0",
    "Contig_BLAST": "8FCE00",
    "Flanking_Context": "92D050",
    "Methods_References": "595959",
    "Metabolic_Integration": "008080"
}

# Dynamically override CONFIG with user values from config/config.yaml if it exists
import os
import yaml
config_yaml_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "config", "config.yaml"))
if os.path.exists(config_yaml_path):
    try:
        with open(config_yaml_path, "r", encoding="utf-8") as f:
            user_config = yaml.safe_load(f)
            if user_config and isinstance(user_config, dict):
                for k, v in user_config.items():
                    if isinstance(v, dict) and k in CONFIG and isinstance(CONFIG[k], dict):
                        CONFIG[k].update(v)
                    else:
                        CONFIG[k] = v
    except Exception:
        pass

