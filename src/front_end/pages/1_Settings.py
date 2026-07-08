import streamlit as st
import sys
import os

_src_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)
    
from front_end.utils.state_manager import init_session_state

st.set_page_config(
    page_title="Redolarium | Settings",
    page_icon="⚙️",
    layout="wide"
)

init_session_state()

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
    st.markdown("### BGC Execution Mode")
    st.markdown("<p class='tooltip-text'>Note: If you select 'Pause after prediction', you must remain on the Terminal page to select specific BGCs once the initial detection phase completes.</p>", unsafe_allow_html=True)
    
    bgc_opt = st.radio(
        "Select BGC run strategy:",
        (
            "Option A: Analyze all predicted BGCs automatically",
            "Option B: Pause after prediction to select specific BGCs for downstream analysis"
        ),
        index=0 if "Option A" in st.session_state.bgc_option else 1
    )
    st.session_state.bgc_option = bgc_opt

st.markdown("---")
if st.button("Proceed to Terminal & Execution"):
    st.success("Settings saved! Please select '2_Terminal' from the sidebar to start the pipeline.")
