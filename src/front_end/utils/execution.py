import subprocess
import threading
import os
import sys
import time
import queue

def build_command_args(state):
    """
    Translates the Streamlit session state into CLI arguments for main.py.
    """
    # Locate main.py relative to execution.py (src/front_end/utils/execution.py -> src/redolarium/main.py)
    frontend_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    src_dir = os.path.dirname(frontend_dir)
    main_script_path = os.path.join(src_dir, "redolarium", "main.py")
    
    cmd = [sys.executable, main_script_path]
    
    cmd.extend(["--query", state.query_file_path])
        
    if state.ref_type != "BLAST":
        cmd.extend(["--ref", state.ref_file_path])
        
    cmd.extend(["--out", state.tmp_out_dir])
    
    if state.email_id:
        cmd.extend(["--email", state.email_id])
        
    cmd.extend(["--cores", str(state.num_cores)])
    
    if state.analysis_query_vs_ref:
        cmd.append("--run-blast")
        
    if state.analysis_bgc:
        if "Analyze all" in state.bgc_option:
            cmd.extend(["--target-bgc", "all"])
            
    return cmd

def log_reader_thread(process, log_list, status_dict, tmp_out_dir):
    """
    Background thread to read stdout/stderr from the process and populate the log_list.
    Implements fault tolerance by avoiding blocking indefinitely if the process crashes silently.
    """
    try:
        # iter() with readline reads until EOF. By using stderr=STDOUT, errors are captured here too.
        for line in iter(process.stdout.readline, ''):
            if not line:
                break
                
            status_dict['last_log_time'] = time.time()
                
            # Sanitize path to prevent leaking system paths
            if tmp_out_dir and tmp_out_dir in line:
                line = line.replace(tmp_out_dir, "[SECURE_OUTPUT_DIR]")
                
            log_list.append(line)
            
            # Detect Deterministic Interactive BGC Prompt
            if "[INPUT_REQUIRED]" in line or "Select BGC number(s) to analyze" in line:
                status_dict['needs_input'] = True
                
        # Wait for the process to actually terminate to get the return code safely
        try:
            process.wait(timeout=10) # Give it 10 seconds to wrap up after stdout closes
        except subprocess.TimeoutExpired:
            log_list.append("\n[WARNING] Process stdout closed but process did not exit. Forcing termination.\n")
            process.kill()
            process.wait()
            
        status_dict['process_finished'] = True
        status_dict['returncode'] = process.returncode
        
        if process.returncode != 0:
            log_list.append(f"\n[CRITICAL ERROR] Subprocess exited with code {process.returncode}\n")
            
    except Exception as e:
        log_list.append(f"\n[SRE THREAD ERROR] Log reader thread crashed: {e}\n")
        status_dict['process_finished'] = True
        status_dict['returncode'] = -1

def kill_process(process, log_list, status_dict):
    """
    Safely kills the subprocess if it hangs or the user forces a stop.
    """
    try:
        process.kill()
        log_list.append("\n[SRE SYSTEM] Process forcibly terminated by user/timeout.\n")
        status_dict['process_finished'] = True
        status_dict['returncode'] = -9
    except Exception as e:
        log_list.append(f"\n[SRE SYSTEM] Failed to kill process: {e}\n")

def start_subprocess(cmd, tmp_out_dir):
    """
    Starts the subprocess and the reader thread with strict error boundaries.
    """
    log_list = []
    status_dict = {
        'needs_input': False, 
        'process_finished': False, 
        'returncode': None,
        'last_log_time': time.time(),
        'error_start': None
    }
    
    try:
        # Use unbuffered output and merge stderr into stdout for unified parsing
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stdin=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        
        thread = threading.Thread(
            target=log_reader_thread, 
            args=(process, log_list, status_dict, tmp_out_dir),
            daemon=True
        )
        thread.start()
        
        return process, log_list, status_dict
        
    except Exception as e:
        # Hard failure to spawn the process (e.g. invalid executable, permissions)
        log_list.append(f"[SRE CRITICAL FAILURE] Failed to launch backend pipeline subprocess: {e}\n")
        status_dict['process_finished'] = True
        status_dict['returncode'] = -1
        return None, log_list, status_dict

def upload_file_externally(file_path):
    """Uploads a file to transfer.sh to bypass GitHub branch protection timeouts."""
    import requests
    import os
    import gzip
    
    filename = os.path.basename(file_path) + ".gz"
    try:
        with open(file_path, 'rb') as f:
            content = f.read()
        compressed = gzip.compress(content)
        
        r = requests.put(f"https://transfer.sh/{filename}", data=compressed)
        if r.status_code == 200:
            return r.text.strip()
        else:
            print(f"External upload failed: {r.status_code} {r.text}")
    except Exception as e:
        print(f"upload exception: {e}")
    return None

def trigger_github_action(state):
    """Triggers the remote GitHub Actions workflow via repository_dispatch."""
    import requests
    import streamlit as st
    import os
    
    if getattr(state, 'query_type', "") != "NCBI Accession ID" and (not hasattr(state, 'query_file_path') or not os.path.exists(state.query_file_path)):
        return False, "Query sequence file is missing."
        
    token = None
    try:
        token = st.secrets.get("GITHUB_TOKEN", None)
    except FileNotFoundError:
        return False, "secrets.toml not found. Please create .streamlit/secrets.toml with your GITHUB_TOKEN."
        
    if not token:
        return False, "GITHUB_TOKEN is missing in .streamlit/secrets.toml"
        
    query_url = ""
    query_acc = ""
    if state.get('query_type', "") == "NCBI Accession ID":
        query_acc = state.get('query_accession', "")
    else:
        query_url = upload_file_externally(state.get('query_file_path', ""))
        if not query_url:
            return False, "Failed to upload query sequence externally."
            
    ref_acc = ""
    ref_url = ""
    if state.get('ref_type', "") == "NCBI Accession ID":
        ref_acc = state.get('ref_accession', "")
    elif state.get('ref_type', "") == "Local Reference Sequence":
        import zipfile
        import os
        zip_path = os.path.join(state.get('tmp_out_dir', ""), "references.zip")
        with zipfile.ZipFile(zip_path, 'w') as z:
            for rp in state.get('ref_file_path', "").split(","):
                if rp and os.path.exists(rp):
                    z.write(rp, os.path.basename(rp))
        ref_url = upload_file_externally(zip_path)
        if not ref_url:
            return False, "Failed to upload reference sequences externally."
        
    repo = "mini-mator/Redolarium"
    url = f"https://api.github.com/repos/{repo}/dispatches"
    
    headers = {
        "Accept": "application/vnd.github+json",
        "Authorization": f"Bearer {token}",
        "X-GitHub-Api-Version": "2022-11-28"
    }
    
    opts = []
    if state.get('adv_bgc_gap'): opts.extend(["--adv-bgc-gap", str(state.get('adv_bgc_gap'))])
    if state.get('adv_bgc_flank'): opts.extend(["--adv-bgc-flank", str(state.get('adv_bgc_flank'))])
    if state.get('adv_hgt_zscore'): opts.extend(["--adv-hgt-zscore", str(state.get('adv_hgt_zscore'))])
    if state.get('adv_hgt_tnf_window'): opts.extend(["--adv-hgt-tnf-window", str(state.get('adv_hgt_tnf_window'))])
    if state.get('adv_promoter_window'): opts.extend(["--adv-promoter-window", str(state.get('adv_promoter_window'))])
    if state.get('adv_homology_identity'): opts.extend(["--adv-homology-identity", str(state.get('adv_homology_identity'))])
    if state.get('adv_homology_evalue'): opts.extend(["--adv-homology-evalue", str(state.get('adv_homology_evalue'))])
    if state.get('adv_clock_multiplier'): opts.extend(["--adv-clock-multiplier", str(state.get('adv_clock_multiplier'))])
    if state.get('adv_metabolic_mode'): opts.extend(["--adv-metabolic-mode", str(state.get('adv_metabolic_mode'))])
    if state.get('adv_cargo_identity'): opts.extend(["--adv-cargo-identity", str(state.get('adv_cargo_identity'))])
    if state.get('adv_cargo_coverage'): opts.extend(["--adv-cargo-coverage", str(state.get('adv_cargo_coverage'))])
    if state.get('adv_cargo_bitscore'): opts.extend(["--adv-cargo-bitscore", str(state.get('adv_cargo_bitscore'))])
    
    payload = {
        "event_type": "run_redolarium",
        "client_payload": {
            "job_id": state.get('job_id', ""),
            "query_url": query_url,
            "query_acc": query_acc,
            "ref_acc": ref_acc,
            "ref_url": ref_url,
            "email": state.get('email_id', ""),
            "target_bgc": "all" if state.get('analysis_bgc', False) else "none",
            "run_blast": "true" if state.get('analysis_query_vs_ref', False) else "false",
            "opts": " ".join(opts)
        }
    }
    
    try:
        r = requests.post(url, json=payload, headers=headers)
        if r.status_code == 204:
            return True, "Success"
        else:
            return False, f"GitHub API Error: {r.status_code} - {r.text}"
    except Exception as e:
        return False, str(e)

def download_github_artifact(artifact_url, token, tmp_out_dir):
    """Downloads and extracts an artifact zip from GitHub Actions."""
    import requests
    import zipfile
    import io
    
    headers = {
        "Authorization": f"Bearer {token}",
        "Accept": "application/vnd.github+json"
    }
    
    try:
        r = requests.get(artifact_url, headers=headers, stream=True)
        if r.status_code == 200:
            z = zipfile.ZipFile(io.BytesIO(r.content))
            z.extractall(tmp_out_dir)
            return True
        else:
            print(f"Failed to download artifact: {r.status_code} {r.text}")
    except Exception as e:
        print(f"Error extracting artifact: {e}")
    return False

def cleanup_github_jobs_branch(job_id, token):
    """(Deprecated) Cleanup function. Left for compatibility."""
    return True