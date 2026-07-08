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

# NCBI Academic Aesthetic CSS
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
    .disclaimer {
        font-size: 0.85em;
        color: #555;
        border-left: 3px solid #0b3d91;
        padding-left: 10px;
        margin-bottom: 20px;
    }
</style>
""", unsafe_allow_html=True)

st.title("Redolarium Pipeline Setup")
st.markdown("""
Welcome to **Redolarium v3.0.0**, a universal genomic, metabolic, and BGC analysis pipeline.
Please configure your analysis parameters below.
""")

st.markdown("""
<div class='disclaimer'>
    <strong>Data Retention Policy:</strong> All results and submitted sequences are strictly kept for 7 days and then automatically and permanently deleted from our servers.
    <br>
    <a href="https://github.com/Redolarium/Redolarium/blob/main/LICENSE" target="_blank">View Software License</a>
</div>
""", unsafe_allow_html=True)

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
    if not st.session_state.email_id:
        st.error("An Email ID is required to proceed.")
    elif query_type == "Input Sequence File" and not st.session_state.query_file_content:
        st.error("Please upload a Query Sequence file.")
    elif query_type == "NCBI Accession ID" and not st.session_state.query_accession:
        st.error("Please enter a valid Query NCBI Accession ID.")
    elif ref_type == "Local Reference Sequence" and not st.session_state.ref_file_contents:
        st.error("Please upload at least one Reference Sequence file.")
    elif ref_type == "NCBI Accession ID" and not st.session_state.ref_accession:
        st.error("Please enter valid Reference NCBI Accession IDs.")
    else:
        st.success("Configuration saved! Please select '1_Settings' from the sidebar to continue.")
