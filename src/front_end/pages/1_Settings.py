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
if st.button("Proceed to Terminal & Execution"):
    if not (st.session_state.analysis_genome_annotation or 
            st.session_state.analysis_query_vs_ref or 
            st.session_state.analysis_cazymes or 
            st.session_state.analysis_bgc):
        st.error("Validation Error: Please select at least one analysis module to proceed. Please change it.")
        st.stop()
    else:
        st.switch_page("pages/2_Terminal.py")
