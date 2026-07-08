import streamlit as st
import sys
import os
import time
import tempfile

_src_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)

from front_end.utils.state_manager import init_session_state
from front_end.utils.execution import build_command_args, start_subprocess
from redolarium.main import download_genbank_from_ncbi

st.set_page_config(
    page_title="Redolarium | Terminal",
    page_icon="🖥️",
    layout="wide"
)

init_session_state()

st.markdown("""
<style>
    .terminal-box {
        background-color: #1e1e1e;
        color: #00ff00;
        font-family: 'Courier New', Courier, monospace;
        padding: 15px;
        border-radius: 5px;
        height: 500px;
        overflow-y: scroll;
        white-space: pre-wrap;
    }
    h1, h2, h3 {
        color: #0b3d91;
    }
</style>
""", unsafe_allow_html=True)

st.title("Execution Terminal")

if not st.session_state.email_id:
    st.error("Missing required inputs. Please go back to the Input Configuration page.")
    st.stop()

# 1. Setup Temporary Directory and Input Files
if not st.session_state.tmp_out_dir:
    st.session_state.tmp_out_dir = tempfile.mkdtemp(prefix="redo_")
    
if not st.session_state.process_started:
    st.session_state.process_started = True
    
    with st.spinner("Preparing input files..."):
        # Handle Query
        if st.session_state.query_type == "Input Sequence File":
            query_path = os.path.join(st.session_state.tmp_out_dir, st.session_state.query_file_name)
            with open(query_path, "wb") as f:
                f.write(st.session_state.query_file_content)
            st.session_state.query_file_path = query_path
        else:
            # NCBI Accession ID
            acc = st.session_state.query_accession
            st.info(f"Downloading Query Accession {acc} from NCBI...")
            downloaded_path = download_genbank_from_ncbi(acc, st.session_state.email_id, st.session_state.tmp_out_dir)
            if not downloaded_path:
                st.error("Failed to download Query from NCBI. Please check your connection or Accession ID.")
                st.stop()
            st.session_state.query_file_path = downloaded_path

        # Handle Reference
        if st.session_state.ref_type == "Local Reference Sequence":
            ref_paths = []
            for name, content in st.session_state.ref_file_contents:
                rp = os.path.join(st.session_state.tmp_out_dir, name)
                with open(rp, "wb") as f:
                    f.write(content)
                ref_paths.append(rp)
            st.session_state.ref_file_path = ",".join(ref_paths)
        elif st.session_state.ref_type == "NCBI Accession ID":
            accs = [p.strip() for p in st.session_state.ref_accession.split(",")]
            ref_paths = []
            for acc in accs:
                st.info(f"Downloading Reference Accession {acc} from NCBI...")
                rp = download_genbank_from_ncbi(acc, st.session_state.email_id, st.session_state.tmp_out_dir)
                if rp:
                    ref_paths.append(rp)
            if not ref_paths:
                st.error("Failed to download References from NCBI.")
                st.stop()
            st.session_state.ref_file_path = ",".join(ref_paths)
        else:
            st.session_state.ref_file_path = "" # BLAST mode

    # Build command and start subprocess
    cmd = build_command_args(st.session_state)
    st.session_state.logs.append(f"$ {' '.join(cmd)}\n")
    
    process, log_list, status_dict = start_subprocess(cmd, st.session_state.tmp_out_dir)
    
    # Store references in session state so we can access them across reruns
    st.session_state.subprocess_obj = process
    st.session_state.subprocess_log_list = log_list
    st.session_state.subprocess_status_dict = status_dict

# Ensure we have the references
if 'subprocess_log_list' in st.session_state:
    log_list = st.session_state.subprocess_log_list
    status_dict = st.session_state.subprocess_status_dict
    process = st.session_state.subprocess_obj
else:
    st.error("Subprocess state lost. Please restart.")
    st.stop()

# 2. Render Terminal output
terminal_placeholder = st.empty()

# Helper function to render logs
def render_terminal():
    display_text = "".join(log_list[-100:]) # Show last 100 lines to prevent lag
    terminal_placeholder.markdown(f"<div class='terminal-box'>{display_text}</div>", unsafe_allow_html=True)

# 3. Handle Interactive BGC Input Phase
if status_dict.get('needs_input', False):
    render_terminal()
    st.warning("Pipeline Paused: Please select BGC(s) for downstream analysis.")
    
    # Simple form to take input
    with st.form("bgc_input_form"):
        user_input = st.text_input("Enter BGC number(s) (e.g., '1', '1, 2', 'A' for all, 'Q' to quit):")
        submitted = st.form_submit_button("Submit")
        if submitted:
            if user_input.strip():
                # Write to stdin
                try:
                    process.stdin.write(f"{user_input.strip()}\n")
                    process.stdin.flush()
                    # Reset flag
                    status_dict['needs_input'] = False
                    log_list.append(f"> {user_input.strip()}\n")
                    st.rerun()
                except Exception as e:
                    st.error(f"Error writing to process: {e}")
            else:
                st.error("Input cannot be empty.")
else:
    # 4. Auto-refresh Loop for live logs
    if not status_dict.get('process_finished', False):
        render_terminal()
        time.sleep(1)
        st.rerun()
    else:
        # Process Finished
        render_terminal()
        if status_dict.get('returncode') == 0:
            st.success("Pipeline Execution Completed Successfully!")
            if st.button("Proceed to Results & Download"):
                st.switch_page("pages/3_Results.py")
        else:
            st.error(f"Pipeline Terminated with errors (Exit Code: {status_dict.get('returncode')})")
            if st.button("Proceed to Results (Partial)"):
                st.switch_page("pages/3_Results.py")
