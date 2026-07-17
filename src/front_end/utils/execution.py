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

def upload_to_github_jobs_branch(file_path, job_id, token):
    """Uploads a file directly to the GitHub jobs branch."""
    import requests
    import base64
    import os
    
    try:
        with open(file_path, 'rb') as f:
            content = f.read()
        
        b64_content = base64.b64encode(content).decode('utf-8')
        filename = f"job_{job_id}_query.gbk"
        url = f"https://api.github.com/repos/mini-mator/Redolarium/contents/{filename}"
        
        headers = {
            "Authorization": f"Bearer {token}",
            "Accept": "application/vnd.github.v3+json"
        }
        
        payload = {
            "message": f"Upload query for job {job_id}",
            "content": b64_content,
            "branch": "jobs"
        }
        
        r = requests.put(url, json=payload, headers=headers)
        if r.status_code in [200, 201]:
            data = r.json()
            return data["content"]["download_url"]
        else:
            err_msg = f"GitHub upload failed: {r.status_code} {r.text}"
            print(err_msg)
            with open("upload_error.txt", "w") as f_err:
                f_err.write(err_msg)
    except Exception as e:
        err_msg = f"upload exception: {e}"
        print(err_msg)
        with open("upload_error.txt", "w") as f_err:
            f_err.write(err_msg)
    return None

def trigger_github_action(state):
    """Triggers the remote GitHub Actions workflow via repository_dispatch."""
    import requests
    import streamlit as st
    
    if not hasattr(state, 'query_file_path') or not os.path.exists(state.query_file_path):
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
    if getattr(state, 'query_type', "") == "NCBI Accession ID":
        query_acc = state.query_accession
    else:
        query_url = upload_to_github_jobs_branch(state.query_file_path, state.job_id, token)
        if not query_url:
            return False, "Failed to upload query sequence to GitHub"
            
    ref_acc = ""
    if getattr(state, 'ref_type', "") == "NCBI Accession ID":
        ref_acc = state.ref_accession
        
    repo = "mini-mator/Redolarium"
    url = f"https://api.github.com/repos/{repo}/dispatches"
    
    headers = {
        "Accept": "application/vnd.github+json",
        "Authorization": f"Bearer {token}",
        "X-GitHub-Api-Version": "2022-11-28"
    }
    
    payload = {
        "event_type": "run_redolarium",
        "client_payload": {
            "job_id": state.job_id,
            "query_url": query_url,
            "query_acc": query_acc,
            "ref_acc": ref_acc,
            "email": state.email_id,
            "target_bgc": "all" if state.analysis_bgc else "none",
            "run_blast": "true" if state.analysis_query_vs_ref else "false",
            "adv_bgc_gap": getattr(state, 'adv_bgc_gap', 10000),
            "adv_bgc_flank": getattr(state, 'adv_bgc_flank', 10000),
            "adv_hgt_zscore": getattr(state, 'adv_hgt_zscore', 2.0),
            "adv_hgt_tnf_window": getattr(state, 'adv_hgt_tnf_window', 5000),
            "adv_promoter_window": getattr(state, 'adv_promoter_window', 600),
            "adv_homology_identity": getattr(state, 'adv_homology_identity', 40.0),
            "adv_homology_evalue": getattr(state, 'adv_homology_evalue', 1e-5),
            "adv_clock_multiplier": getattr(state, 'adv_clock_multiplier', 15.0),
            "adv_metabolic_mode": getattr(state, 'adv_metabolic_mode', 'fast'),
            "adv_cargo_identity": getattr(state, 'adv_cargo_identity', 0.3),
            "adv_cargo_coverage": getattr(state, 'adv_cargo_coverage', 0.3),
            "adv_cargo_bitscore": getattr(state, 'adv_cargo_bitscore', 0.4)
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
    """Deletes the temporary query file from the jobs branch to keep the repo clean."""
    import requests
    
    filename = f"job_{job_id}_query.gbk"
    url = f"https://api.github.com/repos/mini-mator/Redolarium/contents/{filename}?ref=jobs"
    
    headers = {
        "Authorization": f"Bearer {token}",
        "Accept": "application/vnd.github.v3+json"
    }
    
    # 1. Get the SHA of the file (required for deletion)
    r = requests.get(url, headers=headers)
    if r.status_code == 200:
        sha = r.json()["sha"]
        
        # 2. Delete the file
        del_payload = {
            "message": f"Cleanup query for job {job_id}",
            "sha": sha,
            "branch": "jobs"
        }
        del_r = requests.delete(f"https://api.github.com/repos/mini-mator/Redolarium/contents/{filename}", json=del_payload, headers=headers)
        if del_r.status_code == 200:
            return True
    return False