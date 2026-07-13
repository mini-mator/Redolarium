import os
import sys
import glob
import time
import shutil
import hashlib
import pandas as pd
from rich.live import Live
from rich.table import Table
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, TimeRemainingColumn
import concurrent.futures
import subprocess
from multiprocessing import Manager

BENCHMARK_DIR = "F:/IITk/Redolarium/redolarium_github/benchmarks/data"
MANIFEST = os.path.join(BENCHMARK_DIR, "benchmark_manifest.csv")
RESULTS_DIR = os.path.join(BENCHMARK_DIR, "results")
HASH_FILE = os.path.join(BENCHMARK_DIR, "code_hash.txt")
SRC_DIR = "F:/IITk/Redolarium/redolarium_github/src/redolarium"

def get_code_hash():
    """Hash all .py files in src/redolarium to detect code changes."""
    hasher = hashlib.md5()
    py_files = sorted(glob.glob(os.path.join(SRC_DIR, "**/*.py"), recursive=True))
    for f in py_files:
        if "network_retry.py" in f: continue # ignore retry wrapper itself
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
            # Code changed! Restart everything!
            if os.path.exists(RESULTS_DIR):
                shutil.rmtree(RESULTS_DIR)
            os.makedirs(RESULTS_DIR, exist_ok=True)
            with open(HASH_FILE, "w") as f:
                f.write(current_hash)
            return True # Restarted
    else:
        os.makedirs(RESULTS_DIR, exist_ok=True)
        with open(HASH_FILE, "w") as f:
            f.write(current_hash)
    return False

def run_genome(worker_id, row, shared_status):
    acc = row["Accession"]
    gbk_path = row["File_Path"]
    out_dir = os.path.join(RESULTS_DIR, acc)
    
    # Pausable check
    final_report = os.path.join(out_dir, "final_report.csv")
    if os.path.exists(final_report):
        shared_status[worker_id] = "[green]Finished {acc} (Cached)[/green]".format(acc=acc)
        return True

    shared_status[worker_id] = "[cyan]Starting {acc}...[/cyan]".format(acc=acc)
    
    # Run the Redolarium pipeline via subprocess
    cmd = [
        sys.executable, "-u", "F:/IITk/Redolarium/redolarium_github/src/redolarium/main.py",
        "--query", gbk_path,
        "--target-bgc", "all",   # Bypasses the interactive console prompt in Stage 4
        "--out", out_dir,
        "--no-download",
        "--email", "benchmark_redolarium@iitk.ac.in"
    ]
    
    # We pipe stdout to intercept the stages to update the dashboard dynamically
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    
    full_output = []
    for line in process.stdout:
        line = line.strip()
        if line:
            full_output.append(line)
        if "Stage" in line:
            shared_status[worker_id] = "[yellow]Running {acc}: {line}[/yellow]".format(acc=acc, line=line)
            
    process.wait()
    
    if process.returncode == 0:
        shared_status[worker_id] = "[green]Completed {acc}[/green]".format(acc=acc)
        return True
    else:
        shared_status[worker_id] = "[red]Failed {acc}[/red]".format(acc=acc)
        # Append the full traceback/output to a master error log
        error_log = os.path.join(BENCHMARK_DIR, "benchmark_failures.log")
        with open(error_log, "a") as f:
            f.write("\n{0}\n".format('='*50))
            f.write("FAILURE: {acc} (Worker {worker_id})\n".format(acc=acc, worker_id=worker_id))
            f.write("{0}\n".format('='*50))
            f.write("\n".join(full_output))
            f.write("\n{0}\n".format('='*50))
        return False

def generate_table(shared_status, overall_progress_str):
    table = Table(show_header=True, header_style="bold magenta", expand=True)
    table.add_column("100-Genome Stress Test Dashboard", justify="left")
    
    table.add_row("[bold blue]{0}[/bold blue]".format(overall_progress_str))
    table.add_row("")
    
    for i in range(5):
        status = shared_status.get(i, "Idle")
        table.add_row("Instance {0}: {1}".format(i+1, status))
        
    return table

def main():
    if not os.path.exists(MANIFEST):
        print("Manifest not found. Please run scrape_benchmark_genomes.py first.")
        sys.exit(1)
        
    df = pd.read_csv(MANIFEST)
    
    # Cap at 100 genomes for the stress test
    if len(df) > 100:
        df = df.head(100)
    
    restarted = verify_code_state()
    
    manager = Manager()
    shared_status = manager.dict()
    for i in range(5):
        shared_status[i] = "Initializing..."
        
    with Live(generate_table(shared_status, "Starting..."), refresh_per_second=4) as live:
        if restarted:
            live.console.print("[bold red]Code changes detected! Wiped old cache and restarting all genomes.[/bold red]")
            time.sleep(2)
            
        with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:
            futures = {}
            for i, (_, row) in enumerate(df.iterrows()):
                worker_id = i % 5
                futures[executor.submit(run_genome, worker_id, row, shared_status)] = row["Accession"]
                
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
                
                prog_str = "Progress: [{completed}/{total}] | Elapsed: {elapsed:.1f}m | ETA: {remaining:.1f}m".format(
                    completed=completed, total=total, elapsed=elapsed/60.0, remaining=remaining/60.0)
                
                live.update(generate_table(shared_status, prog_str))

if __name__ == "__main__":
    main()
