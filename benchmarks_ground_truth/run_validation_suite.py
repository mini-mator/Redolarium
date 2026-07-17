#!/usr/bin/env python3
import os
import sys
import glob
import time
import shutil
import hashlib
import pandas as pd
from rich.live import Live
from rich.table import Table
from rich.progress import Progress
import concurrent.futures
import subprocess
from multiprocessing import Manager
import json

BASE_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(BASE_DIR, "data")
MANIFEST = os.path.join(DATA_DIR, "ground_truth_manifest.csv")
LEDGER = os.path.join(DATA_DIR, "accession_ledger.csv")
RESULTS_DIR = os.path.join(DATA_DIR, "results")
HASH_FILE = os.path.join(DATA_DIR, "code_hash.txt")
SRC_DIR = os.path.join(os.path.dirname(BASE_DIR), "src", "redolarium")

def get_code_hash():
    """Hash all .py files in src/redolarium to detect code changes."""
    hasher = hashlib.md5()
    py_files = sorted(glob.glob(os.path.join(SRC_DIR, "**/*.py"), recursive=True))
    for f in py_files:
        if "network_retry.py" in f: continue
        with open(f, "rb") as file:
            hasher.update(file.read())
    return hasher.hexdigest()

def verify_code_state():
    """Wipe results if the source code was modified between runs."""
    current_hash = get_code_hash()
    if os.path.exists(HASH_FILE):
        with open(HASH_FILE, "r") as f:
            old_hash = f.read().strip()
        if old_hash != current_hash:
            if os.path.exists(RESULTS_DIR):
                if os.name == "nt":
                    os.system(f'rmdir /s /q "{RESULTS_DIR}"')
                else:
                    os.system(f'rm -rf "{RESULTS_DIR}"')
            os.makedirs(RESULTS_DIR, exist_ok=True)
            with open(HASH_FILE, "w") as f:
                f.write(current_hash)
            return True
    else:
        os.makedirs(RESULTS_DIR, exist_ok=True)
        with open(HASH_FILE, "w") as f:
            f.write(current_hash)
    return False

def run_genome(worker_id, row, shared_status):
    acc = row["Genome_Accession"]
    category = row.get("Category", "Unknown")
    gbk_path = os.path.join(DATA_DIR, "genomes", f"{acc}.gbk")
    out_dir = os.path.join(RESULTS_DIR, category, acc)
    
    if not os.path.exists(gbk_path):
        shared_status[worker_id] = f"[red]Missing file for {acc}[/red]"
        return False

    final_report = os.path.join(out_dir, "final_report.csv")
    if os.path.exists(final_report):
        shared_status[worker_id] = f"[green]Finished {acc} (Cached)[/green]"
        return True

    shared_status[worker_id] = f"[cyan]Starting {acc}...[/cyan]"
    
    # Run the Redolarium pipeline via subprocess
    cmd = [
        sys.executable, "-u", os.path.join(SRC_DIR, "main.py"),
        "--query", gbk_path,
        "--target-bgc", "all",
        "--out", out_dir,
        "--no-download",
        "--email", "benchmark_redolarium@iitk.ac.in"
    ]
    
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    
    full_output = []
    for line in process.stdout:
        line = line.strip()
        if line:
            full_output.append(line)
        if "Stage" in line:
            shared_status[worker_id] = f"[yellow]Running {acc}: {line}[/yellow]"
            
    process.wait()
    
    # Save the raw execution log for FAIR reproducibility
    log_file = os.path.join(out_dir, "pipeline_execution.log")
    os.makedirs(out_dir, exist_ok=True)
    with open(log_file, "w") as f:
        f.write(f"COMMAND: {' '.join(cmd)}\n\n")
        f.write("\n".join(full_output))
    
    if process.returncode == 0:
        shared_status[worker_id] = f"[green]Completed {acc}[/green]"
        return True
    else:
        shared_status[worker_id] = f"[red]Failed {acc}[/red]"
        error_log = os.path.join(DATA_DIR, "validation_failures.log")
        with open(error_log, "a") as f:
            f.write(f"\n{'='*50}\n")
            f.write(f"FAILURE: {acc} (Worker {worker_id})\n")
            f.write(f"{'='*50}\n")
            f.write("\n".join(full_output[-50:])) # Log the last 50 lines
            f.write(f"\n{'='*50}\n")
        return False

def generate_table(shared_status, overall_progress_str):
    table = Table(show_header=True, header_style="bold magenta", expand=True)
    table.add_column("Redolarium Ground Truth Validation Dashboard", justify="left")
    
    table.add_row(f"[bold blue]{overall_progress_str}[/bold blue]")
    table.add_row("")
    
    for i in range(5):
        status = shared_status.get(i, "Idle")
        table.add_row(f"Instance {i+1}: {status}")
        
    return table

def main():
    if not os.path.exists(LEDGER):
        print("Accession ledger not found. Please run harvest_golden_genomes.py first.")
        sys.exit(1)
        
    df = pd.read_csv(LEDGER)
    restarted = verify_code_state()
    
    manager = Manager()
    shared_status = manager.dict()
    for i in range(5):
        shared_status[i] = "Initializing..."
        
    with Live(generate_table(shared_status, "Starting..."), refresh_per_second=4) as live:
        if restarted:
            live.console.print("[bold red]Code changes detected! Wiped old validation cache and restarting all genomes.[/bold red]")
            time.sleep(2)
            
        with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:
            futures = {}
            for i, (_, row) in enumerate(df.iterrows()):
                worker_id = i % 5
                futures[executor.submit(run_genome, worker_id, row, shared_status)] = row["Genome_Accession"]
                
            not_done = set(futures.keys())
            completed = 0
            total = len(df)
            start_time = time.time()
            
            while not_done:
                done, not_done = concurrent.futures.wait(not_done, timeout=0.25)
                
                for future in done:
                    completed += 1
                
                elapsed = time.time() - start_time
                avg_time = elapsed / completed if completed > 0 else 0
                remaining = avg_time * (total - completed)
                
                prog_str = f"Progress: [{completed}/{total}] | Elapsed: {elapsed/60.0:.1f}m | ETA: {remaining/60.0:.1f}m"
                live.update(generate_table(shared_status, prog_str))
                
    print("\nPhase 2 Complete. Pipeline execution finished across golden genomes.")
    print("Ready for generate_truth_analytics.py to compile the Confusion Matrix.")

if __name__ == "__main__":
    main()
