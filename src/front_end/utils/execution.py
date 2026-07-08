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
