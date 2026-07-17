import streamlit as st
import sys
import os

_src_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)
    
from front_end.utils.state_manager import init_session_state

st.set_page_config(
    page_title="Redolarium | Settings",
    page_icon="",
    layout="wide"
)

init_session_state()

# Initialize default advanced settings
if "adv_bgc_gap" not in st.session_state:
    st.session_state.adv_bgc_gap = 10000
if "adv_bgc_flank" not in st.session_state:
    st.session_state.adv_bgc_flank = 10000
if "adv_hgt_zscore" not in st.session_state:
    st.session_state.adv_hgt_zscore = 2.0
if "adv_hgt_tnf_window" not in st.session_state:
    st.session_state.adv_hgt_tnf_window = 5000
if "adv_promoter_window" not in st.session_state:
    st.session_state.adv_promoter_window = 600
if "adv_homology_identity" not in st.session_state:
    st.session_state.adv_homology_identity = 40.0
if "adv_homology_evalue" not in st.session_state:
    st.session_state.adv_homology_evalue = 1e-5
if "adv_clock_multiplier" not in st.session_state:
    st.session_state.adv_clock_multiplier = 15.0
if "adv_metabolic_mode" not in st.session_state:
    st.session_state.adv_metabolic_mode = "fast"
if "adv_cargo_identity" not in st.session_state:
    st.session_state.adv_cargo_identity = 0.3
if "adv_cargo_coverage" not in st.session_state:
    st.session_state.adv_cargo_coverage = 0.3
if "adv_cargo_bitscore" not in st.session_state:
    st.session_state.adv_cargo_bitscore = 0.4

st.markdown("""
<style>
    .reportview-container {
        background: #f8f9fa;
    }
    h1, h2, h3 {
        color: #0b3d91;
        font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
    }
    .stButton>button {
        background-color: #0b3d91;
        color: white;
        border-radius: 4px;
        border: none;
        padding: 10px 24px;
        font-weight: bold;
    }
    .stButton>button:hover {
        background-color: #082d6b;
        color: white;
    }
    .tooltip-text {
        font-size: 0.9em;
        color: #666;
        font-style: italic;
    }
</style>
""", unsafe_allow_html=True)

st.title("Analysis Settings")
st.markdown("Select the downstream analysis modules you wish to execute.")

# 1. Complete Genome Annotation
st.session_state.analysis_genome_annotation = st.checkbox(
    "Complete Genome Annotation", 
    value=st.session_state.analysis_genome_annotation,
    help="Runs structural and functional annotation across the entire genome."
)

# 2. Query vs Reference Comparison
st.session_state.analysis_query_vs_ref = st.checkbox(
    "Query vs Reference Comparison", 
    value=st.session_state.analysis_query_vs_ref,
    help="Triggers remote BLASTp homolog searches and comparative genomics against reference sequences."
)

# 3. CAZymes Analysis
st.session_state.analysis_cazymes = st.checkbox(
    "CAZymes Analysis", 
    value=st.session_state.analysis_cazymes,
    help="Screens for Carbohydrate-Active enZymes (CAZymes) in the annotated genome."
)

# 4. BGC Analysis
st.session_state.analysis_bgc = st.checkbox(
    "BGC (Biosynthetic Gene Cluster) Analysis", 
    value=st.session_state.analysis_bgc,
    help="Predicts secondary metabolite regions using antiSMASH integration and structural HMM profiling."
)

if st.session_state.analysis_bgc:
    st.info("BGC Execution Mode: Analyze all predicted BGCs automatically (Required for Cloud/Batch execution)")
    st.session_state.bgc_option = "Option A: Analyze all predicted BGCs automatically"

st.markdown("---")
st.markdown("### Advanced Sensitivity Parameters (Proceed with Caution)")
st.markdown("Redolarium enforces peer-reviewed biological thresholds. Modifying these values changes the sensitivity of the pipeline. **Citations for default parameters are provided below.**")

with st.expander("Configure Advanced Settings"):
    st.markdown("#### BGC & Genomic Island Delineation")
    st.session_state.adv_bgc_gap = st.number_input(
        "BGC Gene Clustering Gap (bp)", 
        value=st.session_state.adv_bgc_gap,
        help="Default based on antiSMASH 7.0 (Blin et al. 2023, doi:10.1093/nar/gkad344)"
    )
    st.session_state.adv_bgc_flank = st.number_input(
        "BGC Flanking Window Size (bp)", 
        value=st.session_state.adv_bgc_flank,
        help="Boundary definitions for genomic islands (Cimermancic et al. 2014, doi:10.1016/j.cell.2014.06.034)"
    )
    
    st.markdown("#### Horizontal Gene Transfer & Promoters")
    st.session_state.adv_hgt_zscore = st.number_input(
        "HGT GC Z-Score Cutoff", 
        value=st.session_state.adv_hgt_zscore,
        step=0.1,
        help="Default derived from horizontal transfer mapping (Langille & Bhatt 2006, doi:10.1186/1471-2164-7-142)"
    )
    st.session_state.adv_hgt_tnf_window = st.number_input(
        "HGT TNF Sliding Window (bp)", 
        value=st.session_state.adv_hgt_tnf_window,
        help="Tetranucleotide frequency resolution limit (Teeling et al. 2004, doi:10.1093/bioinformatics/bth054)"
    )
    st.session_state.adv_promoter_window = st.number_input(
        "Promoter Search Window Upstream (bp)", 
        value=st.session_state.adv_promoter_window,
        help="Standard upstream regulatory window for bacterial promoters."
    )
    
    st.markdown("#### Phylogeny & Homology")
    st.session_state.adv_homology_identity = st.number_input(
        "Homology Min Identity (%)", 
        value=st.session_state.adv_homology_identity,
        step=1.0,
        help="Homology detection zone (Rost 1999, doi:10.1093/protein/12.2.85)"
    )
    st.session_state.adv_homology_evalue = st.number_input(
        "Homology Max E-value (Format: 1e-5)", 
        value=float(st.session_state.adv_homology_evalue),
        format="%.1e",
        help="Sequence similarity statistical significance (Pearson 2013, doi:10.1002/0471250953.bi0301s42)"
    )
    st.session_state.adv_clock_multiplier = st.number_input(
        "Speciation Clock Multiplier (Mya / 1% div)", 
        value=st.session_state.adv_clock_multiplier,
        step=1.0,
        help="Bacterial molecular clock calibration (Ochman et al. 2000, doi:10.1016/S0092-8674(00)80405-8)"
    )
    
    st.markdown("#### Metabolism & Cargo")
    st.session_state.adv_metabolic_mode = st.selectbox(
        "Metabolic Pipeline Mode", 
        ["fast", "deep"],
        index=0 if st.session_state.adv_metabolic_mode == "fast" else 1,
        help="KEGG metabolic reference network (Kanehisa et al. 2023, doi:10.1093/nar/gkac963)"
    )
    st.markdown("Cargo Scoring Weights (Identity / Coverage / Bitscore) - Empirically optimized via MIBiG database.")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.session_state.adv_cargo_identity = st.number_input("Identity Wt", value=st.session_state.adv_cargo_identity, step=0.1)
    with col2:
        st.session_state.adv_cargo_coverage = st.number_input("Coverage Wt", value=st.session_state.adv_cargo_coverage, step=0.1)
    with col3:
        st.session_state.adv_cargo_bitscore = st.number_input("Bitscore Wt", value=st.session_state.adv_cargo_bitscore, step=0.1)

st.markdown("---")
if st.button("Proceed to Terminal & Execution"):
    if not (st.session_state.analysis_genome_annotation or 
            st.session_state.analysis_query_vs_ref or 
            st.session_state.analysis_cazymes or 
            st.session_state.analysis_bgc):
        st.error("Validation Error: Please select at least one analysis module to proceed. Please change it.")
        st.stop()
    else:
        st.session_state.process_started = False
        st.session_state.cloud_triggered = False
        if 'subprocess_log_list' in st.session_state:
            del st.session_state.subprocess_log_list
        if 'subprocess_status_dict' in st.session_state:
            del st.session_state.subprocess_status_dict
        if 'subprocess_obj' in st.session_state:
            del st.session_state.subprocess_obj
        st.session_state.logs = []
        
        st.switch_page("pages/2_Terminal.py")
