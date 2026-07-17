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
    page_title="Redolarium | Execution",
    page_icon="",
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

st.title("Pipeline Execution")

if not st.session_state.email_id:
    st.error("Missing required inputs. Please go back to the Input Configuration page.")
    st.stop()

# 1. Setup Temporary Directory and Input Files
if not st.session_state.tmp_out_dir:
    st.session_state.tmp_out_dir = tempfile.mkdtemp(prefix="redo_")
    
if not st.session_state.process_started:
    st.session_state.process_started = True
    st.session_state.job_id = str(os.urandom(8).hex())
    
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

    if st.session_state.local_mode:
        st.info("Running in Local Developer Mode. Check terminal logs below.")
        cmd = build_command_args(st.session_state)
        st.session_state.logs.append(f"$ {' '.join(cmd)}\n")
        
        process, log_list, status_dict = start_subprocess(cmd, st.session_state.tmp_out_dir)
        st.session_state.subprocess_obj = process
        st.session_state.subprocess_log_list = log_list
        st.session_state.subprocess_status_dict = status_dict
    else:
        st.info("Initiating Cloud Compute via GitHub Actions...")
        from front_end.utils.execution import trigger_github_action
        success, msg = trigger_github_action(st.session_state)
        if not success:
            st.error(f"Failed to trigger Cloud Compute: {msg}")
            st.stop()
        st.session_state.cloud_triggered = True

# --- RENDER LOGIC ---
if st.session_state.local_mode:
    # 2. Render Terminal output (Local Mode)
    if 'subprocess_log_list' not in st.session_state:
        st.error("Subprocess state lost. Please restart.")
        st.stop()
        
    log_list = st.session_state.subprocess_log_list
    status_dict = st.session_state.subprocess_status_dict
    
    terminal_placeholder = st.empty()
    display_text = "".join(log_list[-100:])
    terminal_placeholder.markdown(f"<div class='terminal-box'>{display_text}</div>", unsafe_allow_html=True)
    
    if not status_dict.get('process_finished', False):
        time.sleep(1)
        st.rerun()
    else:
        if status_dict.get('returncode') == 0:
            st.success("Local Pipeline Execution Completed Successfully!")
            if st.button("Proceed to Results Dashboard"):
                st.switch_page("pages/3_Results.py")
        else:
            st.error(f"Pipeline Terminated with errors (Exit Code: {status_dict.get('returncode')})")
            if st.button("Proceed to Results (Partial)"):
                st.switch_page("pages/3_Results.py")
else:
    # 3. Cloud Mode Polling
    st.markdown("### Cloud Execution in Progress")
    st.markdown("Your job has been dispatched to GitHub Actions. This page will automatically poll for results.")
    st.markdown("**Please keep this window open.** Depending on your genome size, antiSMASH can take anywhere from 10 minutes to 3 hours.")
    
    status_placeholder = st.empty()
    steps_placeholder = st.empty()
    cloud_terminal_placeholder = st.empty()
    
    # Check for completion via the Artifact API
    import requests
    import toml
    
    token = ""
    secrets_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".streamlit", "secrets.toml")
    if os.path.exists(secrets_path):
        with open(secrets_path, "r") as f:
            token = toml.load(f).get("GITHUB_TOKEN", "")
            
    job_id = st.session_state.job_id
    artifact_ready = False
    
    if token:
        headers = {"Accept": "application/vnd.github.v3+json", "Authorization": f"Bearer {token}"}
        
        try:
            # Poll for the Artifact
            artifacts_req = requests.get(f"https://api.github.com/repos/mini-mator/Redolarium/actions/artifacts?name=redolarium_results_{job_id}", headers=headers)
            if artifacts_req.status_code == 200:
                artifacts_data = artifacts_req.json()
                if artifacts_data.get("total_count", 0) > 0:
                    artifact = artifacts_data["artifacts"][0]
                    if not artifact.get("expired"):
                        artifact_ready = True
                        result_url = artifact["archive_download_url"]
                        status_placeholder.success(f"Job Finished! Downloading Artifact...")
                        
                        # Automatically download and extract
                        st.info("Extracting results package and cleaning up...")
                        from front_end.utils.execution import download_github_artifact, cleanup_github_jobs_branch
                        
                        if download_github_artifact(result_url, token, st.session_state.tmp_out_dir):
                            cleanup_github_jobs_branch(job_id, token)
                            st.success("Ready!")
                            if st.button("Proceed to Results Dashboard"):
                                st.switch_page("pages/3_Results.py")
                        else:
                            st.error("Failed to download or extract results package.")
        except Exception as e:
            st.error(f"Error querying artifacts: {e}")
            
        if not artifact_ready:
            # If not finished, query GitHub API to get real-time step progress!
            status_placeholder.info(f"Job {job_id} is running... waiting for completion marker.")
            
            if token:
                try:
                    # Get latest run
                    runs_req = requests.get("https://api.github.com/repos/mini-mator/Redolarium/actions/runs?per_page=1", headers=headers)
                    if runs_req.status_code == 200:
                        runs_data = runs_req.json()
                        if runs_data.get("workflow_runs"):
                            latest_run_id = runs_data["workflow_runs"][0]["id"]
                            
                            # Get jobs for this run
                            jobs_req = requests.get(f"https://api.github.com/repos/mini-mator/Redolarium/actions/runs/{latest_run_id}/jobs", headers=headers)
                            if jobs_req.status_code == 200:
                                jobs_data = jobs_req.json()
                                if jobs_data.get("jobs"):
                                    gh_job_id = jobs_data["jobs"][0]["id"]
                                    steps = jobs_data["jobs"][0].get("steps", [])
                                    
                                    # Render steps
                                    steps_md = "#### Pipeline Stages:\n"
                                    for step in steps:
                                        name = step.get("name", "Unknown Step")
                                        # Skip GitHub internal setup steps to keep UI clean
                                        if name in ["Set up job", "Complete job", "Post Checkout Repository", "Post Set up Python 3.10"]:
                                            continue
                                            
                                        status = step.get("status")
                                        conclusion = step.get("conclusion")
                                        
                                        if status == "completed":
                                            if conclusion == "success":
                                                steps_md += f"- [x] {name}\n"
                                            else:
                                                steps_md += f"- [ ] **FAILED**: {name}\n"
                                        elif status == "in_progress":
                                            steps_md += f"- [ ] **(In Progress...)** {name}\n"
                                        else:
                                            steps_md += f"- [ ] {name}\n"
                                    
                                    steps_placeholder.markdown(steps_md)
                                    
                                    # Fetch real-time logs
                                    logs_req = requests.get(f"https://api.github.com/repos/mini-mator/Redolarium/actions/jobs/{gh_job_id}/logs", headers=headers, allow_redirects=True)
                                    if logs_req.status_code == 200:
                                        lines = logs_req.text.split('\n')
                                        # Filter out gh action timestamps to make it look like a normal terminal
                                        clean_lines = []
                                        for line in lines[-200:]:
                                            import re
                                            # Remove GitHub timestamp "2026-07-17T09:05:44.3313829Z "
                                            clean_line = re.sub(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d+Z\s", "", line)
                                            clean_lines.append(clean_line)
                                            
                                        display_text = "\n".join(clean_lines)
                                        cloud_terminal_placeholder.markdown(f"<div class='terminal-box'>{display_text}</div>", unsafe_allow_html=True)
                except Exception as e:
                    pass # Silently fallback if API polling fails
            
            time.sleep(10)
            st.rerun()
