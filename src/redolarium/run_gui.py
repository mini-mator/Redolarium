import os
import sys
import json
import time
import queue
import logging
import threading
import traceback
import webbrowser
from http.server import HTTPServer, SimpleHTTPRequestHandler
from Bio import SeqIO

# Ensure parent directory is in sys.path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from redolarium.utils import setup_logging, compile_bgc_excel_report, write_word_report
from redolarium.annotation import run_annotation_pipeline
from redolarium.metabolism import run_metabolism_pipeline
from redolarium.screening import run_screening_pipeline
from redolarium.bgc_analysis import detect_all_bgcs, evaluate_targeted_bgc
from redolarium.target_bgc_analysis import run_target_bgc_analysis
from redolarium.evolution import run_evolution_pipeline
from redolarium.phylogeny import run_phylogeny_pipeline
from redolarium.docking import run_docking_pipeline
from redolarium.linkage import run_linkage_pipeline

# Global App State
STATE = {
    "status": "idle",       # idle, running, completed, failed
    "stage": "Not Started", # Annotation, Metabolism, etc.
    "percentage": 0,        # 0 to 100
    "error": None,
    "elapsed": 0,
    "estimated_remaining": 0
}

# Logs state
LOGS = []
LOGS_EVENT = threading.Event()
PRESETS_FILE = os.path.join(parent_dir, "gui_presets.json")

class WebLogHandler(logging.Handler):
    def emit(self, record):
        msg = self.format(record)
        LOGS.append(msg)
        LOGS_EVENT.set()

# Add a lock for thread-safe state modification
STATE_LOCK = threading.Lock()

def update_state(status=None, stage=None, percentage=None, error=None, elapsed=None, remaining=None):
    with STATE_LOCK:
        if status is not None: STATE["status"] = status
        if stage is not None: STATE["stage"] = stage
        if percentage is not None: STATE["percentage"] = percentage
        if error is not None: STATE["error"] = error
        if elapsed is not None: STATE["elapsed"] = elapsed
        if remaining is not None: STATE["estimated_remaining"] = remaining

def load_presets():
    if os.path.exists(PRESETS_FILE):
        try:
            with open(PRESETS_FILE, "r") as f:
                return json.load(f)
        except Exception:
            pass
    
    # Default presets
    return {
        "Default Quick Run": {
            "query": "CP031675.1.gb",
            "ref": "AL009126.3.gb",
            "out": "results",
            "email": "researcher@example.com",
            "cores": 4,
            "dock_runs": 1,
            "run_blast": False,
            "target_bgc": "BGC_01"
        },
        "Streptomyces BGC Focus": {
            "query": "AL009126.3.gb",
            "ref": "",
            "out": "results_streptomyces",
            "email": "researcher@example.com",
            "cores": 4,
            "dock_runs": 3,
            "run_blast": True,
            "target_bgc": "BGC_01"
        }
    }

def save_preset(name, config):
    presets = load_presets()
    presets[name] = config
    try:
        with open(PRESETS_FILE, "w") as f:
            json.dump(presets, f, indent=4)
        return True
    except Exception:
        return False

def run_pipeline_task(config):
    query = config.get("query")
    ref = config.get("ref") or None
    out_dir = config.get("out", "results")
    email = config.get("email", "researcher@example.com")
    cores = int(config.get("cores", 4))
    dock_runs = int(config.get("dock_runs", 1))
    run_blast = bool(config.get("run_blast", False))
    target_bgc_id = config.get("target_bgc", "BGC_01")
    
    start_time = time.time()
    
    try:
        # Create output directories
        subdirs = [
            "docking_images", "phylogeny_trees", "tabular_data", 
            "metabolic_pathways", "comparative_genomics", "hgt_evolution", 
            "screening_cazymes", "bgc_motifs"
        ]
        for sub in subdirs:
            os.makedirs(os.path.join(out_dir, sub), exist_ok=True)
            
        logger = logging.getLogger("redolarium")
        logger.setLevel(logging.INFO)
        
        # Attach our custom web logger handler if not already attached
        handler_exists = False
        for h in logger.handlers:
            if isinstance(h, WebLogHandler):
                handler_exists = True
                break
        
        if not handler_exists:
            wh = WebLogHandler()
            wh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
            logger.addHandler(wh)
            
        logger.info("Initializing Redolarium Web UI Pipeline run...")
        logger.info(f"Query file: {query}")
        logger.info(f"Reference file: {ref}")
        
        # 1. Annotation & Comparative Genomics
        update_state(status="running", stage="Annotation & Comparative Genomics", percentage=10, remaining=180)
        logger.info("Stage 1/7: Running Annotation & Comparative Genomics...")
        orth_map, ref_strains, identities, sim_matrix = run_annotation_pipeline(
            query, ref, cores, email, out_dir, logger
        )
        
        # 2. Metabolic Reconstruction
        update_state(stage="Metabolic Pathway Reconstruction", percentage=25, remaining=140)
        logger.info("Stage 2/7: Running Metabolic Pathway Reconstruction...")
        metab_data = run_metabolism_pipeline(query, orth_map, out_dir, logger)
        
        # 3. Screening & Cazyme Mapping
        update_state(stage="CAZyme & Cargo Screening", percentage=38, remaining=110)
        logger.info("Stage 3/7: Running Screening pipelines...")
        screened_cargo = run_screening_pipeline(query, orth_map, out_dir, logger)
        
        # 4. BGC Delineation
        update_state(stage="BGC Detection", percentage=50, remaining=85)
        logger.info("Stage 4/7: Detecting Biosynthetic Gene Clusters...")
        bgc_regions = detect_all_bgcs(query, logger)
        
        if not bgc_regions:
            logger.error("No Biosynthetic Gene Clusters detected in genome!")
            update_state(status="failed", stage="Error", percentage=100, error="No BGCs found in input genome.", remaining=0)
            return
            
        # Select target BGC
        selected_bgc = None
        for bgc in bgc_regions:
            if bgc["BGC_ID"].lower() == target_bgc_id.lower():
                selected_bgc = bgc
                break
                
        if not selected_bgc:
            selected_bgc = bgc_regions[0]
            logger.warning(f"Target BGC {target_bgc_id} not found. Defaulting to first detected: {selected_bgc['BGC_ID']}.")
            
        bgc_id = selected_bgc["BGC_ID"]
        logger.info(f"Targeting cluster for analysis: {bgc_id}...")
        
        # 5a. Delineate core & flanking layout
        update_state(stage="Target BGC Delineation", percentage=60, remaining=65)
        region_genes, flanking_genes, _, _ = evaluate_targeted_bgc(
            query, ref, orth_map, selected_bgc, out_dir, logger
        )
        
        # 5b. Promoter & prophage analysis
        update_state(stage="Promoter & Phage Hits Analysis", percentage=70, remaining=50)
        promoter_records, phage_hits, blast_rows = run_target_bgc_analysis(
            query, selected_bgc, run_blast, email, out_dir, logger
        )
        
        # 5c. Evolutionary HGT signatures
        update_state(stage="Evolutionary HGT Tracing", percentage=78, remaining=35)
        hgt_results = run_evolution_pipeline(query, selected_bgc, out_dir, logger)
        
        # 5d. Phylogenetics
        update_state(stage="Marker-Gene Phylogenetics", percentage=85, remaining=20)
        run_phylogeny_pipeline(query, ref_strains, identities, sim_matrix, out_dir, logger, bgc_blast_results=blast_rows)
        
        # 5e. Molecular Docking
        update_state(stage="Molecular Docking Simulation", percentage=90, remaining=10)
        # Perform actual run_docking_pipeline
        run_docking_pipeline(query, selected_bgc, region_genes, out_dir, logger)
        
        # 5f. Precursor ATP stoichiometry
        update_state(stage="Stoichiometry & Cost Linkage Mapping", percentage=95, remaining=5)
        stoichiometry_data = run_linkage_pipeline(query, region_genes, selected_bgc, out_dir, logger)
        
        # 6. Report generation
        update_state(stage="Report Compilation", percentage=98, remaining=2)
        xls_filename = f"{bgc_id}_BGC_Analysis_Metabolism_Integrated.xlsx"
        xls_path = os.path.join(out_dir, xls_filename)
        query_record = max(list(SeqIO.parse(query, "genbank")), key=lambda r: len(r.seq))
        query_org = query_record.annotations.get("organism", "Query Isolate")
        
        logger.info("Compiling integrated Excel workbook...")
        compile_bgc_excel_report(
            selected_bgc, orth_map, ref_strains, metab_data,
            region_genes, flanking_genes, promoter_records, phage_hits,
            hgt_results, stoichiometry_data, xls_path, query_org, logger,
            bgc_blast_results=blast_rows
        )
        
        doc_filename = f"{bgc_id}_BGC_Report.docx"
        logger.info("Drafting Word report document...")
        sections = [
            ("1. Executive Summary", [
                f"This document presents the complete genomic, metabolic, and secondary metabolism Biosynthetic "
                f"Gene Cluster (BGC) analysis reconstructed for cluster {bgc_id} in {query_org}."
            ]),
            ("2. Biosynthetic Gene Cluster Core Architecture", [
                f"The core region of {bgc_id} contains {len(region_genes)} core genes. Ortholog mapping indicates key biosynthetic roles."
            ])
        ]
        write_word_report(doc_filename, bgc_id, query_org, sections, out_dir)
        
        logger.info("Redolarium web pipeline execution completed successfully!")
        elapsed = round(time.time() - start_time, 1)
        update_state(status="completed", stage="Completed", percentage=100, elapsed=elapsed, remaining=0)
        
    except Exception as e:
        err_msg = f"Pipeline Exception: {e}\n{traceback.format_exc()}"
        logging.getLogger("redolarium").error(err_msg)
        update_state(status="failed", stage="Error", percentage=100, error=str(e), remaining=0)


class RedolariumHandler(SimpleHTTPRequestHandler):
    def translate_path(self, path):
        # Override to serve from static/ folder
        if path.startswith("/static/") or path == "/":
            relative_path = path.lstrip("/")
            if not relative_path:
                relative_path = "static/index.html"
            return os.path.join(parent_dir, relative_path)
        return super().translate_path(path)

    def do_GET(self):
        # Custom API Endpoints
        if self.path == "/" or self.path.startswith("/static/"):
            return super().do_GET()
            
        elif self.path == "/api/status":
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.end_headers()
            self.wfile.write(json.dumps(STATE).encode("utf-8"))
            
        elif self.path == "/api/presets":
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.end_headers()
            self.wfile.write(json.dumps(load_presets()).encode("utf-8"))
            
        elif self.path == "/api/sysinfo":
            # Detect CPU cores and list GenBank/Fasta files in workspace root
            import multiprocessing
            try:
                cores = multiprocessing.cpu_count()
            except Exception:
                cores = 4
                
            files = []
            for f in os.listdir(parent_dir):
                if f.endswith((".gb", ".gbk", ".fasta", ".fa", ".fna")):
                    files.append(f)
                    
            sysinfo = {
                "cpu_cores": cores,
                "detected_files": files
            }
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.end_headers()
            self.wfile.write(json.dumps(sysinfo).encode("utf-8"))
            
        elif self.path == "/api/logs":
            # Server-Sent Events Endpoint
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.send_header("Cache-Control", "no-cache")
            self.send_header("Connection", "keep-alive")
            self.end_headers()
            
            # Send current logs
            last_index = 0
            for log in LOGS:
                self.wfile.write(f"data: {json.dumps(log)}\n\n".encode("utf-8"))
                last_index += 1
            self.wfile.flush()
            
            # Stream incoming logs
            while STATE["status"] == "running" or LOGS_EVENT.wait(timeout=2.0):
                if STATE["status"] != "running" and last_index >= len(LOGS):
                    break
                LOGS_EVENT.clear()
                while last_index < len(LOGS):
                    log = LOGS[last_index]
                    try:
                        self.wfile.write(f"data: {json.dumps(log)}\n\n".encode("utf-8"))
                        self.wfile.flush()
                    except (ConnectionError, BrokenPipeError):
                        return # Client disconnected
                    last_index += 1
                    
        elif self.path.startswith("/api/results"):
            # Explore output results folder dynamically
            query_params = self.path.split("?")
            out_dir = "results"
            if len(query_params) > 1:
                parts = query_params[1].split("=")
                if len(parts) > 1 and parts[0] == "out":
                    out_dir = parts[1]
            
            res_dir = os.path.join(parent_dir, out_dir)
            if not os.path.exists(res_dir):
                self.send_response(404)
                self.end_headers()
                self.wfile.write(b"Results folder not found")
                return
                
            files_dict = {}
            for root, dirs, files in os.walk(res_dir):
                rel_path = os.path.relpath(root, res_dir)
                if rel_path == ".":
                    rel_path = "root"
                files_dict[rel_path] = [f for f in files if not f.startswith(".")]
                
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.end_headers()
            self.wfile.write(json.dumps({"out_dir": out_dir, "structure": files_dict}).encode("utf-8"))
            
        elif self.path.startswith("/api/files/"):
            # Download/serve result file
            # Format: /api/files/<out_dir>/<rel_path>
            file_path = self.path[len("/api/files/"):]
            # basic protection against path traversal
            abs_path = os.path.abspath(os.path.join(parent_dir, file_path))
            if not abs_path.startswith(os.path.abspath(parent_dir)):
                self.send_response(403)
                self.end_headers()
                return
                
            if not os.path.exists(abs_path) or os.path.isdir(abs_path):
                self.send_response(404)
                self.end_headers()
                return
                
            # Detect content type
            content_type = "application/octet-stream"
            if abs_path.endswith(".png"): content_type = "image/png"
            elif abs_path.endswith(".svg"): content_type = "image/svg+xml"
            elif abs_path.endswith(".jpg") or abs_path.endswith(".jpeg"): content_type = "image/jpeg"
            elif abs_path.endswith(".xlsx"): content_type = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            elif abs_path.endswith(".docx"): content_type = "application/vnd.openxmlformats-officedocument.wordprocessingml.document"
            elif abs_path.endswith(".txt") or abs_path.endswith(".json") or abs_path.endswith(".log"): content_type = "text/plain"
            
            self.send_response(200)
            self.send_header("Content-Type", content_type)
            # Add attachment header for excel/word files
            if abs_path.endswith((".xlsx", ".docx")):
                self.send_header("Content-Disposition", f'attachment; filename="{os.path.basename(abs_path)}"')
            self.end_headers()
            
            with open(abs_path, "rb") as f:
                self.wfile.write(f.read())
        else:
            self.send_response(404)
            self.end_headers()

    def do_POST(self):
        content_length = int(self.headers['Content-Length'])
        post_data = self.rfile.read(content_length)
        
        if self.path == "/api/run":
            if STATE["status"] == "running":
                self.send_response(400)
                self.end_headers()
                self.wfile.write(b"Pipeline already running")
                return
                
            config = json.loads(post_data.decode("utf-8"))
            # Clear logs
            LOGS.clear()
            LOGS_EVENT.set()
            
            # Start pipeline thread
            t = threading.Thread(target=run_pipeline_task, args=(config,))
            t.daemon = True
            t.start()
            
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.end_headers()
            self.wfile.write(json.dumps({"message": "Started"}).encode("utf-8"))
            
        elif self.path == "/api/presets":
            data = json.loads(post_data.decode("utf-8"))
            name = data.get("name")
            config = data.get("config")
            if name and config:
                success = save_preset(name, config)
                self.send_response(200 if success else 500)
                self.send_header("Content-Type", "application/json")
                self.end_headers()
                self.wfile.write(json.dumps({"success": success}).encode("utf-8"))
            else:
                self.send_response(400)
                self.end_headers()
        else:
            self.send_response(404)
            self.end_headers()

def main():
    port = 8000
    server_address = ('127.0.0.1', port)
    httpd = HTTPServer(server_address, RedolariumHandler)
    print(f"Redolarium Web GUI server started at http://localhost:{port}")
    
    # Auto-open browser in a thread
    def open_browser():
        time.sleep(1.5)
        webbrowser.open(f"http://localhost:{port}/")
        
    threading.Thread(target=open_browser, daemon=True).start()
    
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down server gracefully...")
    except Exception as e:
        print(f"\nFATAL: Unexpected GUI server crash: {e}")
        import traceback
        traceback.print_exc()
    finally:
        httpd.server_close()

if __name__ == "__main__":
    main()
