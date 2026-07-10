import os
import sys
import time
import queue
import threading
import logging
import streamlit as st
from Bio import SeqIO
from redolarium.utils import setup_logging, compile_bgc_excel_report, write_word_report
from redolarium.annotation import run_annotation_pipeline
from redolarium.metabolism import run_metabolism_pipeline
from redolarium.screening import run_screening_pipeline
from redolarium.bgc_analysis import detect_all_bgcs, evaluate_targeted_bgc
from redolarium.target_bgc_analysis import run_target_bgc_analysis
from redolarium.evolution import run_evolution_pipeline
from redolarium.phylogeny import run_phylogeny_pipeline
from redolarium.linkage import run_linkage_pipeline

class QueueLoggingHandler(logging.Handler):
    """Thread-safe logging handler that writes log messages to a queue."""
    def __init__(self, log_queue):
        super().__init__()
        self.log_queue = log_queue

    def emit(self, record):
        msg = self.format(record)
        self.log_queue.put(msg)

def run_pipeline_thread(query, ref, out_dir, email, cores, run_blast, target_bgc_id, log_queue, status_event):
    """Worker function that runs the complete pipeline steps in a background thread."""
    try:
        # Set up folders
        subdirs = [
            "phylogeny_trees", "tabular_data", 
            "metabolic_pathways", "comparative_genomics", "hgt_evolution", 
            "screening_cazymes", "bgc_motifs"
        ]
        for sub in subdirs:
            os.makedirs(os.path.join(out_dir, sub), exist_ok=True)
            
        logger = logging.getLogger("redolarium")
        logger.setLevel(logging.INFO)
        
        # Attach Queue logger handler
        qh = QueueLoggingHandler(log_queue)
        qh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(qh)
        
        logger.info("Starting background pipeline thread...")
        
        # 1. Annotation & Comparative Genomics
        orth_map, ref_strains, identities, sim_matrix = run_annotation_pipeline(
            query, ref, cores, email, out_dir, logger
        )
        
        # 2. Metabolic Reconstruction
        metab_data = run_metabolism_pipeline(query, orth_map, out_dir, logger)
        
        # 3. Screening & Profit Mapping
        screened_cargo = run_screening_pipeline(query, orth_map, out_dir, logger)
        
        # 4. BGC Delineation
        bgc_regions = detect_all_bgcs(query, logger)
        
        if not bgc_regions:
            logger.error("No Biosynthetic Gene Clusters detected in genome!")
            status_event["error"] = "No BGCs found."
            status_event["complete"] = True
            return
            
        # Select target BGC
        selected_bgc = None
        for bgc in bgc_regions:
            if bgc["BGC_ID"].lower() == target_bgc_id.lower():
                selected_bgc = bgc
                break
                
        if not selected_bgc:
            selected_bgc = bgc_regions[0]
            logger.warning(f"Target BGC {target_bgc_id} not found. Defaulting to first BGC: {selected_bgc['BGC_ID']}.")
            
        bgc_id = selected_bgc["BGC_ID"]
        logger.info(f"Targeting cluster for analysis: {bgc_id}...")
        
        # 5a. Delineate core & flanking layout
        region_genes, flanking_genes, _, _ = evaluate_targeted_bgc(
            query, ref, orth_map, selected_bgc, out_dir, logger
        )
        
        # 5b. Promoter & prophage analysis
        result_prom = run_target_bgc_analysis(
            query, selected_bgc, run_blast, email, out_dir, logger
        )
        promoter_records = result_prom.prediction.get("promoter_records", [])
        phage_hits = result_prom.prediction.get("phage_hits", [])
        blast_rows = result_prom.prediction.get("blast_rows", [])
        
        # 5c. Evolutionary HGT signatures
        hgt_results = run_evolution_pipeline(query, selected_bgc, out_dir, logger)
        
        # 5d. Phylogenetics
        run_phylogeny_pipeline(query, ref_strains, identities, sim_matrix, out_dir, logger, bgc_blast_results=blast_rows)
        
        # 5e. Precursor ATP stoichiometry
        stoichiometry_data = run_linkage_pipeline(query, region_genes, selected_bgc, out_dir, logger)
        
        # 6. Report generation
        xls_filename = f"{bgc_id}_BGC_Analysis_Metabolism_Integrated.xlsx"
        xls_path = os.path.join(out_dir, xls_filename)
        query_record = SeqIO.read(query, "genbank")
        query_org = query_record.annotations.get("organism", "Query Isolate")
        
        logger.info("Compiling Excel workbook...")
        compile_bgc_excel_report(
            selected_bgc, orth_map, ref_strains, metab_data,
            region_genes, flanking_genes, promoter_records, phage_hits,
            hgt_results, stoichiometry_data, xls_path, query_org, logger,
            bgc_blast_results=blast_rows
        )
        
        # Word report
        doc_filename = f"{bgc_id}_BGC_Report.docx"
        logger.info("Drafting Word report...")
        sections = [
            ("1. Executive Summary", [
                f"This document presents the complete genomic, metabolic, and secondary metabolism Biosynthetic "
                f"Gene Cluster (BGC) analysis reconstructed for cluster {bgc_id} in {query_org}."
            ]),
            ("2. Biosynthetic Gene Cluster Core Architecture", [
                f"The core region of {bgc_id} contains {len(region_genes)} core genes. Ortholog mapping indicates key biosynthetic roles."
            ])
        ]
        if not region_genes and not flanking_genes:
            logger.warning(f"Aborting Word report generation for {bgc_id} due to missing structural cluster inputs. Falling back to Excel workbook only.")
        else:
            write_word_report(doc_filename, bgc_id, query_org, sections, out_dir)
            logger.info("Word report drafted successfully.")
        
        logger.info("Pipeline thread completed successfully.")
        status_event["complete"] = True
        
    except Exception as e:
        import traceback
        err_msg = f"Thread exception: {e}\n{traceback.format_exc()}"
        logging.getLogger("redolarium").error(err_msg)
        status_event["error"] = str(e)
        status_event["complete"] = True

def main():
    st.set_page_config(page_title="Redolarium Dashboard", layout="wide")
    st.title("🧬 Redolarium: Genomic and Secondary Metabolism Delineation")
    
    st.sidebar.header("Execution Parameters")
    query_file = st.sidebar.text_input("Query GenBank Path", "sequence_ksm_I.fasta")
    ref_file = st.sidebar.text_input("Reference GenBank Path (Optional)", "")
    out_dir = st.sidebar.text_input("Output Directory", "results_gui")
    email = st.sidebar.text_input("Entrez Email Address", "researcher@example.com")
    
    cores = st.sidebar.slider("CPU Cores for Orthology Alignment", 1, 8, 4)
    run_blast = st.sidebar.checkbox("Run Remote BGC BLAST Conservation", value=False)
    target_bgc = st.sidebar.text_input("Target BGC ID (e.g. BGC_01)", "BGC_01")
    
    # Initialize session status
    if "log_queue" not in st.session_state:
        st.session_state.log_queue = queue.Queue()
    if "status_event" not in st.session_state:
        st.session_state.status_event = {"complete": False, "error": None}
    if "log_text" not in st.session_state:
        st.session_state.log_text = ""
    if "running" not in st.session_state:
        st.session_state.running = False
        
    run_btn = st.sidebar.button("Launch Analysis Pipeline")
    
    if run_btn:
        st.session_state.log_queue = queue.Queue()
        st.session_state.status_event = {"complete": False, "error": None}
        st.session_state.log_text = "Initializing pipeline thread...\n"
        st.session_state.running = True
        
        # Start background thread
        t = threading.Thread(
            target=run_pipeline_thread,
            args=(
                query_file, None if not ref_file else ref_file, out_dir, email,
                cores, run_blast, target_bgc, st.session_state.log_queue, st.session_state.status_event
            )
        )
        t.start()
        st.info("Pipeline started running in background thread. Monitoring logs...")
        
        # Rerun to update the state immediately
        if hasattr(st, "rerun"):
            st.rerun()
        else:
            st.experimental_rerun()
        
    # Logging window & monitoring
    log_area = st.empty()
    
    if st.session_state.running or st.session_state.log_text:
        # Drain queue
        new_logs = False
        while not st.session_state.log_queue.empty():
            try:
                msg = st.session_state.log_queue.get_nowait()
                st.session_state.log_text += msg + "\n"
                new_logs = True
            except queue.Empty:
                break
                
        if st.session_state.log_text:
            log_area.code(st.session_state.log_text)
            
        if st.session_state.status_event["complete"]:
            st.session_state.running = False
            if st.session_state.status_event["error"]:
                st.error(f"Pipeline crashed: {st.session_state.status_event['error']}")
            else:
                st.success("Redolarium analysis run finished successfully! Output saved in: " + out_dir)
        elif st.session_state.running:
            time.sleep(0.5)
            if hasattr(st, "rerun"):
                st.rerun()
            else:
                st.experimental_rerun()

if __name__ == "__main__":
    main()
