# Central configuration and threshold constants for Redolarium Pipeline
CONFIG = {
    # Flanking window sizes for BGC definition
    # Reference: Blin et al. 2023 (antiSMASH 7.0, doi:10.1093/nar/gkad344); Cimermancic et al. 2014 (doi:10.1016/j.cell.2014.06.034)
    "bgc_flank_window": 10000,          # bp flanking coordinates around BGC core
    "promoter_search_upstream": 600,    # bp upstream of ATG for promoter motifs

    # Bacterial regulatory consensus binding motifs structured by taxonomic phylum
    # References: Hawley & McClure 1983 (E. coli σ70); Helmann 1995 (Bacillus σA/σW/σB);
    #             Potvin et al. 2008 (Pseudomonas σ54); Bourn & Babb 1995 (Streptomyces HrdB);
    #             Bayley et al. 2000 (Bacteroides consensus).
    # spacer_min/spacer_max: distance between end of –35 box and start of –10 box
    "sigma_motifs": {
        "Firmicutes": {
            "sigmaA (housekeeping)": {
                "instances_35": ["TTGACA", "TTGATA", "TTGAAA", "TTGACT", "TTGACC"],
                "instances_10": ["TATAAT", "TATACT", "TATATT", "TACTAT"],
                "spacer_min": 16,
                "spacer_max": 18,
                "induction": "Constitutive / Growth-phase dependent"
            },
            "sigmaW (envelope stress)": {
                "instances_35": ["TGGAAA", "TGGAAC", "TGGAAT", "CGGAAA"],
                "instances_10": ["CGTAAT", "CGTTAT", "CGTCAT"],
                "spacer_min": 17,
                "spacer_max": 19,
                "induction": "Alkaline / Envelope Stress"
            },
            "sigmaB (general stress)": {
                "instances_35": ["GTTTTAA", "GTTTTAA", "GTTTTAT"],
                "instances_10": ["GGGTAT", "GGGTAT"],
                "spacer_min": 13,
                "spacer_max": 15,
                "induction": "General Environmental Stress"
            },
            "ComA/AgrA Box (QS)": {
                "instances_35": ["TTGCCT", "TTGCAT", "TTGCAA"],
                "instances_10": ["TATAAT", "TACAAT"],
                "spacer_min": 15,
                "spacer_max": 18,
                "induction": "Peptide-mediated Quorum Sensing"
            }
        },
        "Pseudomonadota": {
            "sigma70 (housekeeping)": {
                "instances_35": ["TTGACA", "TTGATA", "TTGAAA", "TTGACT", "TTGACC"],
                "instances_10": ["TATAAT", "TATACT", "TATATT", "TACTAT"],
                "spacer_min": 16,
                "spacer_max": 18,
                "induction": "Constitutive / Growth-phase dependent"
            },
            "sigmaS (stationary phase)": {
                "instances_35": ["TTGACA", "TTGATA", "TTGACT"],
                "instances_10": ["CTATACT", "CTATACA"],
                "spacer_min": 16,
                "spacer_max": 18,
                "induction": "Stationary Phase / General Stress"
            },
            "sigma54 (nitrogen stress)": {
                "instances_35": ["TGGCAC", "TGGCTC"],
                "instances_10": ["TTGCA", "TTGCT"],
                "spacer_min": 10,
                "spacer_max": 12,
                "induction": "Nitrogen / Carbon Starvation"
            },
            "LuxR/LasR Box (QS)": {
                "instances_35": ["ACCTGC", "ACCTGTA", "ACCTGAT"],
                "instances_10": ["GCAGGT", "GCAGGC"],
                "spacer_min": 4,
                "spacer_max": 8,
                "induction": "AHL-mediated Quorum Sensing"
            }
        },
        "Actinomycetota": {
            "HrdB (housekeeping)": {
                "instances_35": ["TCGACC", "TTGACA", "TCGACT"],
                "instances_10": ["TAGAAT", "TACAAT", "TATAAT", "GAGAAT"],
                "spacer_min": 16,
                "spacer_max": 19,
                "induction": "Constitutive / Growth-phase dependent"
            },
            "SARP/ArpA Box (QS)": {
                "instances_35": ["TCGAGA", "TCGAGC"],
                "instances_10": ["CTCAGA", "CTCAGC"],
                "spacer_min": 10,
                "spacer_max": 22,
                "induction": "Gamma-butyrolactone Signaling / Secondary Metabolism"
            }
        },
        "Bacteroidota": {
            "Bf-consensus": {
                "instances_35": ["TTG", "TTG", "TAG"],
                "instances_10": ["TAAAAT", "TATAAT", "TACAAT"],
                "spacer_min": 17,
                "spacer_max": 20,
                "induction": "Constitutive / Growth-phase dependent"
            }
        }
    },

    # Reference: Blin et al. 2023 (antiSMASH 7.0, doi:10.1093/nar/gkad344)
    "bgc_clustering_gap_bp": 10000,

    # Horizontal Gene Transfer (HGT) evolutionary divergence thresholds
    # Reference: Lawrence & Ochman 1997 (doi:10.1007/s002399900038) - GC-content amelioration and anomaly detection
    # Z-score >= 2.0 is the standard 95% statistical confidence interval for compositional anomaly
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
        {"module": "HGT TNF-bias", "tool": "Genomic Signature", "ref": "Teeling, H. et al. (2004). Application of tetranucleotide frequencies for fragment clustering.", "doi": "10.1093/bioinformatics/bth054"},
        {"module": "HGT Scores", "tool": "HGT-Composite", "ref": "Ravenhall, M. et al. (2015). Inferring horizontal gene transfer using genomic signatures.", "doi": "10.1371/journal.pcbi.1004167"},
        {"module": "Phylogeny", "tool": "Bio.Phylo", "ref": "Talevich, N. et al. (2012). Bio.Phylo: phylogeny framework in Biopython.", "doi": "10.1186/1471-2105-13-209"},
        {"module": "Molecular Docking", "tool": "RCSB PDB API", "ref": "Berman, H. M. et al. (2000). The Protein Data Bank.", "doi": "10.1093/nar/28.1.235"},
        {"module": "Molecular Docking", "tool": "ChEMBL API", "ref": "Mendez, D. et al. (2019). ChEMBL: towards direct deposition of bioactivity data.", "doi": "10.1093/nar/gky1075"},
        {"module": "Metabolic Linkage", "tool": "KEGG REST API", "ref": "Kanehisa, M. et al. (2023). KEGG database for taxonomy-based pathways.", "doi": "10.1093/nar/gkac963"},
        {"module": "CAZyme Screening", "tool": "dbCAN3", "ref": "Huang, L. et al. (2023). dbCAN3: automated carbohydrate-active enzyme annotation.", "doi": "10.1093/nar/gkad328"},
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
        # Weights empirically optimized via MIBiG grid search (see benchmarks/weight_optimization.py)
        "weights": {
            "identity": 0.3,
            "coverage": 0.3,
            "bitscore": 0.4
        },
        "thresholds": {
            # Reference: Rost 1999 (doi:10.1093/protein/12.2.85) for homology detection zone
            "min_identity": 40.0,
            # Reference: Pearson 2013 (doi:10.1002/0471250953.bi0301s42); Eddy 2011 (doi:10.1371/journal.pcbi.1002195)
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

