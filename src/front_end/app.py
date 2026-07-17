import streamlit as st
import os
import sys

# Ensure project root is in sys path to allow redolarium imports
_src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)

from front_end.utils.state_manager import init_session_state
from front_end.utils.validators import validate_ncbi_accession, parse_comma_separated_accessions

st.set_page_config(
    page_title="Redolarium | Input Configuration",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

init_session_state()

# Professional Academic Aesthetic CSS (NCBI / Elsevier inspired)
st.markdown("""
<style>
    .reportview-container {
        background: #ffffff;
    }
    h1, h2, h3 {
        color: #222222;
        font-family: 'Arial', sans-serif;
        border-bottom: 1px solid #eaeaea;
        padding-bottom: 5px;
    }
    h1 { color: #0b3d91; }
    .stButton>button {
        background-color: #0b3d91;
        color: white;
        border-radius: 3px;
        border: 1px solid #082d6b;
        padding: 6px 16px;
        font-size: 14px;
        font-weight: 600;
    }
    .stButton>button:hover {
        background-color: #082d6b;
        color: white;
        border-color: #061e47;
    }
    .disclaimer {
        font-size: 0.85em;
        color: #444;
        background-color: #f8f9fa;
        border: 1px solid #dee2e6;
        border-left: 4px solid #0b3d91;
        padding: 12px;
        margin-bottom: 20px;
    }
</style>
""", unsafe_allow_html=True)

st.title("Redolarium Pipeline Setup")
st.markdown("""
Welcome to **Redolarium v1.0.0**, a universal genomic, metabolic, and BGC analysis pipeline.
Please configure your analysis parameters below.
""")

st.markdown("""
<div class='disclaimer'>
    <strong>Data Retention Policy:</strong> All results and submitted sequences are kept for 7 days and then automatically deleted from our cloud servers.
    <br>
    <a href="https://github.com/mini-mator/Redolarium/blob/main/LICENSE" target="_blank" style="color: #0b3d91; text-decoration: none;">View Software License</a>
</div>
""", unsafe_allow_html=True)

with st.sidebar:
    st.markdown("### Developer Options")
    st.session_state.local_mode = st.checkbox(
        "Run Locally (Bypass Cloud)", 
        value=st.session_state.local_mode,
        help="If checked, the pipeline will execute natively on your machine instead of using GitHub Actions."
    )
    st.markdown("---")

# 1. Query Sequence
st.header("1. Query Sequence")
query_type = st.radio("Select Query Input Method:", ("Input Sequence File", "NCBI Accession ID"), index=0 if st.session_state.query_type == "Input Sequence File" else 1)
st.session_state.query_type = query_type

if query_type == "Input Sequence File":
    uploaded_query = st.file_uploader("Upload Query GenBank File (.gb)", type=["gb", "gbk"])
    if uploaded_query:
        st.session_state.query_file_content = uploaded_query.getvalue()
        st.session_state.query_file_name = uploaded_query.name
elif query_type == "NCBI Accession ID":
    query_acc = st.text_input("Enter NCBI Accession ID (e.g., CP031675.1)", value=st.session_state.query_accession)
    if query_acc:
        if validate_ncbi_accession(query_acc):
            st.success("Valid NCBI Accession Format.")
            st.session_state.query_accession = query_acc
        else:
            st.error("Invalid NCBI Accession format. Please check and try again.")
            st.session_state.query_accession = ""

# 2. Reference Sequence
st.header("2. Reference Sequence")
ref_type = st.radio("Select Reference Input Method:", ("Local Reference Sequence", "NCBI Accession ID", "BLAST"), index=["Local Reference Sequence", "NCBI Accession ID", "BLAST"].index(st.session_state.ref_type))
st.session_state.ref_type = ref_type

if ref_type == "Local Reference Sequence":
    uploaded_refs = st.file_uploader("Upload Reference GenBank File(s) (.gb)", type=["gb", "gbk"], accept_multiple_files=True)
    if uploaded_refs:
        st.session_state.ref_file_contents = [(f.name, f.getvalue()) for f in uploaded_refs]
elif ref_type == "NCBI Accession ID":
    ref_accs = st.text_input("Enter multiple comma-separated Accession IDs", value=st.session_state.ref_accession)
    if ref_accs:
        valid, invalid = parse_comma_separated_accessions(ref_accs)
        if invalid:
            st.warning(f"The following IDs appear invalid: {', '.join(invalid)}")
        if valid:
            st.success(f"Valid IDs ready for use: {', '.join(valid)}")
            st.session_state.ref_accession = ",".join(valid)
elif ref_type == "BLAST":
    st.info("BLAST mode activated. No explicit reference file required; the pipeline will perform remote homolog searches.")

# 3. System Parameters
st.header("3. System Parameters")
col1, col2 = st.columns(2)
with col1:
    email = st.text_input("Email ID (Required for NCBI queries & Results Delivery)", value=st.session_state.email_id)
    if email:
        st.session_state.email_id = email
with col2:
    cores = st.number_input("Number of CPU Cores to utilize", min_value=1, max_value=64, value=st.session_state.num_cores, step=1)
    st.session_state.num_cores = cores

st.markdown("---")
if st.button("Proceed to Analysis Settings"):
    errors = []
    if not st.session_state.email_id:
        errors.append("An Email ID is required to proceed.")
    if query_type == "Input Sequence File" and not st.session_state.query_file_content:
        errors.append("Please upload a Query Sequence file.")
    if query_type == "NCBI Accession ID" and not st.session_state.query_accession:
        errors.append("Please enter a valid Query NCBI Accession ID.")
    if ref_type == "Local Reference Sequence" and not st.session_state.ref_file_contents:
        errors.append("Please upload at least one Reference Sequence file.")
    if ref_type == "NCBI Accession ID" and not st.session_state.ref_accession:
        errors.append("Please enter valid Reference NCBI Accession IDs.")
    
    if errors:
        for err in errors:
            st.error(f"Validation Error: {err} Please change it.")
        st.stop()
    else:
        st.switch_page("pages/1_Settings.py")
