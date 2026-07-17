import streamlit as st

def init_session_state():
    """Initializes the session state variables required across the multipage app."""
    
    # Page 1: Input Configuration State
    if 'query_type' not in st.session_state:
        st.session_state.query_type = "Input Sequence File"
    if 'query_accession' not in st.session_state:
        st.session_state.query_accession = ""
    if 'query_file_content' not in st.session_state:
        st.session_state.query_file_content = None
    if 'query_file_name' not in st.session_state:
        st.session_state.query_file_name = ""

    if 'ref_type' not in st.session_state:
        st.session_state.ref_type = "Local Reference Sequence"
    if 'ref_accession' not in st.session_state:
        st.session_state.ref_accession = ""
    if 'ref_file_contents' not in st.session_state:
        st.session_state.ref_file_contents = [] # List of tuples (filename, content)
        
    if 'email_id' not in st.session_state:
        st.session_state.email_id = ""
    if 'num_cores' not in st.session_state:
        st.session_state.num_cores = 4
    if 'local_mode' not in st.session_state:
        st.session_state.local_mode = False

    # Page 2: Analysis Settings State
    if 'analysis_genome_annotation' not in st.session_state:
        st.session_state.analysis_genome_annotation = True
    if 'analysis_query_vs_ref' not in st.session_state:
        st.session_state.analysis_query_vs_ref = False
    if 'analysis_cazymes' not in st.session_state:
        st.session_state.analysis_cazymes = True
    if 'analysis_bgc' not in st.session_state:
        st.session_state.analysis_bgc = True
    if 'bgc_option' not in st.session_state:
        st.session_state.bgc_option = "Option A: Analyze all predicted BGCs automatically"

    # Page 3: Terminal / Subprocess State
    if 'process_started' not in st.session_state:
        st.session_state.process_started = False
    if 'logs' not in st.session_state:
        st.session_state.logs = []
    if 'needs_input' not in st.session_state:
        st.session_state.needs_input = False
    if 'process_finished' not in st.session_state:
        st.session_state.process_finished = False
    if 'tmp_out_dir' not in st.session_state:
        st.session_state.tmp_out_dir = ""
    if 'bgc_choices_presented' not in st.session_state:
        st.session_state.bgc_choices_presented = []

    # Page 4: Results State
    if 'zip_path' not in st.session_state:
        st.session_state.zip_path = ""
