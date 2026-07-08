# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Redolarium Contributors, IIT Kanpur
# See LICENSE and THIRD_PARTY_LICENSES.md for full licence details.
# type: ignore
import os
import sys

# --- Path bootstrap: ensure the project root is on sys.path when run directly ---
_project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)
# ---------------------------------------------------------------------------------
import argparse
import json
import logging
from Bio import SeqIO
from redolarium.utils import setup_logging, write_word_report, compile_bgc_excel_report
from redolarium.qc import run_qc_pipeline
from redolarium.annotation import run_annotation_pipeline
from redolarium.metabolism import run_metabolism_pipeline
from redolarium.bgc_analysis import detect_all_bgcs, evaluate_targeted_bgc, run_bgc_pipeline
from redolarium.evolution import run_evolution_pipeline
from redolarium.phylogeny import run_phylogeny_pipeline
from redolarium.docking import run_docking_pipeline
from redolarium.linkage import run_linkage_pipeline
from redolarium.screening import run_screening_pipeline
from redolarium.target_bgc_analysis import run_target_bgc_analysis

def download_genbank_from_ncbi(accession_id, email, dest_dir="."):
    print(f"\nConnecting to NCBI Entrez to download Accession '{accession_id}'...")
    import urllib.request
    import urllib.parse
    
    encoded_id = urllib.parse.quote(accession_id)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={encoded_id}&rettype=gbwithparts&retmode=text&email={email}"
    
    dest_path = os.path.join(dest_dir, f"{accession_id}.gb")
    
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as response:
            content = response.read().decode('utf-8')
            
        if not content.strip() or "Error" in content or "Nothing to return" in content:
            print(f"Error: NCBI did not return data for accession '{accession_id}'.")
            return None
            
        if not content.startswith("LOCUS"):
            snippet = content[:200].replace("\n", " ")
            print(f"Error: Invalid file format returned from NCBI: {snippet}")
            return None
            
        with open(dest_path, "w", encoding="utf-8") as f:
            f.write(content)
            
        print(f"Success! Downloaded and saved GenBank file as '{dest_path}'")
        return dest_path
        
    except Exception as e:
        print(f"Network error downloading from NCBI: {e}")
        return None

def get_genbank_file_interactive(prompt_label, current_email, required=True):
    while True:
        print(f"\n--- Specify {prompt_label} GenBank File ---")
        print("[1] Use local file path")
        print("[2] Download from NCBI by Accession ID")
        if not required:
            print("[3] Skip (no reference)")
        print("[Q] Exit Setup / Quit")
        
        choice = input("Enter choice: ").strip().lower()
        if choice == 'q':
            print("Exiting setup wizard.")
            sys.exit(0)
        elif choice == '1':
            while True:
                path = input(f"Enter path to local {prompt_label} GenBank file(s) or directory (comma-separated if multiple, or 'q' to quit): ").strip().strip("\"'")
                if not path:
                    break
                if path.lower() == 'q':
                    print("Exiting setup wizard.")
                    sys.exit(0)
                if os.path.exists(path):
                    return path, current_email
                else:
                    print(f"Error: File '{path}' not found. Please try again or press Enter to go back.")
        elif choice == '2':
            while True:
                acc = input(f"Enter NCBI Accession ID (e.g., CP031675.1) (or 'q' to quit): ").strip()
                if not acc:
                    break
                if acc.lower() == 'q':
                    print("Exiting setup wizard.")
                    sys.exit(0)
                
                if not current_email or current_email == "researcher@example.com":
                    email_input = input("Enter email for remote NCBI queries (required by NCBI): ").strip()
                    if email_input:
                        current_email = email_input
                
                downloaded_path = download_genbank_from_ncbi(acc, current_email)
                if downloaded_path:
                    return downloaded_path, current_email
                else:
                    print("Error: NCBI download failed. Please double check the accession or try again.")
        elif not required and choice == '3':
            return None, current_email
        else:
            print("Invalid choice. Please choose one of the options listed.")

def generate_snakefile(out_dir, config_dict):
    """
    Generates a standard Snakemake Snakefile representing the execution DAG checkpoints.
    Allows swapping reference databases seamlessly without modifying pipeline code.
    """
    snakefile_path = os.path.join(out_dir, "Snakefile")
    content = f"""# Redolarium Snakemake Execution DAG
# Generated dynamically for reproducibility and audit trails

config = {json.dumps(config_dict, indent=4)}

rule all:
    input:
        os.path.join(config["out_dir"], "tabular_data", "evidence_graph.json")

# Rule 0: Assembly Quality Control
rule qc:
    input:
        query = config["query"]
    output:
        qc_json = os.path.join(config["out_dir"], "tabular_data", "qc_metrics.json")
    run:
        from redolarium.qc import run_qc_pipeline
        from redolarium.utils import setup_logging
        logger = setup_logging(config["out_dir"])
        qc_res = run_qc_pipeline(input.query, config["out_dir"], logger)
        import pickle
        with open(output.qc_json, "wb") as f:
            pickle.dump(qc_res, f)
"""
    with open(snakefile_path, "w", encoding="utf-8") as f:
        f.write(content)
    return snakefile_path

def main():
    parser = argparse.ArgumentParser(
        description="Redolarium v3.0.0 — Universal Genomic, Metabolic & BGC Analysis Pipeline\n"
                    "MIT Licence | Third-party tools retain their own licences (see THIRD_PARTY_LICENSES.md)\n"
                    "KEGG API usage is restricted to academic/non-commercial use (see KEGG ToS).\n"
                    "antiSMASH (AGPL-3.0) is an optional external dependency invoked via Docker.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--query", default=None, help="Path to query GenBank file")
    parser.add_argument("--ref", default=None, help="Path to reference GenBank file, comma-separated list of files, or directory (optional)")
    parser.add_argument("--out", default="results", help="Directory path to save results")
    # NCBI ToS require a valid email address for Entrez/BLAST requests.
    # A placeholder is intentionally not set as default to force users to provide one.
    parser.add_argument("--email", default=None, help="Your email address for NCBI Entrez/BLAST queries (required by NCBI ToS)")
    parser.add_argument("--cores", type=int, default=4, help="Number of CPU cores for parallel execution")
    parser.add_argument("--dock-runs", type=int, default=1, help="Number of molecular docking iterations")
    parser.add_argument("--run-blast", action="store_true", help="Trigger remote BLASTp homolog searches")
    parser.add_argument("--run-docking", action="store_true", help="Trigger molecular docking analysis (Speculative, Experimental & Under Development)")
    parser.add_argument("--target-bgc", default=None, help="Force analysis of a specific BGC ID (e.g. BGC_01)")
    parser.add_argument("--no-download", action="store_true", help="Disable automatic download of reference genomes on genus mismatch")
    parser.add_argument("--force-kegg-refresh", action="store_true", help="Force refresh of cached KEGG pathway, EC, and KO databases")
    parser.add_argument("--use-snakemake", action="store_true", help="Trigger execution using Snakemake API")
    
    args = parser.parse_args()

    # --- NCBI ToS compliance: validate email ---
    _PLACEHOLDER_EMAIL = "researcher@example.com"
    if args.email is None:
        args.email = _PLACEHOLDER_EMAIL
    if args.email == _PLACEHOLDER_EMAIL:
        print(
            "\n[WARNING] No --email supplied. Using placeholder 'researcher@example.com'.\n"
            "  NCBI's Entrez/BLAST Terms of Service require a valid contact email.\n"
            "  Remote NCBI queries (Entrez fetch, BLASTp) may be throttled or blocked.\n"
            "  Please re-run with: --email your@institution.edu\n"
        )

    # If query is not specified, run the interactive setup wizard
    if not args.query:
        print("\n=== Redolarium: Interactive Setup Wizard ===")
        print("Please provide the required parameters below. (Type 'Q' to exit at any time)")
        
        args.query, args.email = get_genbank_file_interactive("Query", args.email, required=True)
        args.ref, args.email = get_genbank_file_interactive("Reference", args.email, required=False)
            
        out_val = input("Enter directory path to save results [default: results, or 'q' to quit]: ").strip()
        if out_val.lower() == 'q':
            print("Exiting setup wizard.")
            sys.exit(0)
        if out_val:
            args.out = out_val.strip("\"'")
            
        if not args.email or args.email == "researcher@example.com":
            email_val = input("Enter email for remote NCBI queries [default: researcher@example.com, or 'q' to quit]: ").strip()
            if email_val.lower() == 'q':
                print("Exiting setup wizard.")
                sys.exit(0)
            if email_val:
                args.email = email_val
            
        cores_val = input("Enter number of CPU cores for parallel execution [default: 4, or 'q' to quit]: ").strip()
        if cores_val.lower() == 'q':
            print("Exiting setup wizard.")
            sys.exit(0)
        if cores_val:
            try:
                args.cores = int(cores_val)
            except ValueError:
                print("Using default (4) due to invalid core count.")
                
        dock_choice = input("Trigger molecular docking analysis? (y/n) [default: n, or 'q' to quit]: ").strip().lower()
        if dock_choice == 'q':
            print("Exiting setup wizard.")
            sys.exit(0)
        args.run_docking = (dock_choice == 'y')
        
        if args.run_docking:
            dock_val = input("Enter number of molecular docking iterations [default: 1, or 'q' to quit]: ").strip()
            if dock_val.lower() == 'q':
                print("Exiting setup wizard.")
                sys.exit(0)
            if dock_val:
                try:
                    args.dock_runs = int(dock_val)
                except ValueError:
                    print("Using default (1) due to invalid iteration count.")
                
        blast_val = input("Trigger remote BLASTp homolog searches? (y/n) [default: n, or 'q' to quit]: ").strip().lower()
        if blast_val == 'q':
            print("Exiting setup wizard.")
            sys.exit(0)
        args.run_blast = (blast_val == 'y')
            
        print("\n=== Interactive Setup Completed ===\n")
    
    out_dir = args.out
    subdirs = [
        "docking_images", "phylogeny_trees", "tabular_data", 
        "metabolic_pathways", "comparative_genomics", "hgt_evolution", 
        "screening_cazymes", "bgc_motifs"
    ]
    for sub in subdirs:
        os.makedirs(os.path.join(out_dir, sub), exist_ok=True)
        
    logger = setup_logging(out_dir)
    logger.info("Initializing Redolarium analytical run...")
    logger.info(f"Query file: {args.query}")
    logger.info(f"Reference file: {args.ref}")
    logger.info(f"Output directory structure initialized under: {out_dir}")
    
    # Generate Snakefile for reproducibility
    config_dict = {
        "query": args.query,
        "ref": args.ref,
        "out_dir": out_dir,
        "email": args.email,
        "cores": args.cores,
        "no_download": args.no_download
    }
    generate_snakefile(out_dir, config_dict)
    logger.info(f"Dynamically generated Snakefile for workflow reproduction at: {os.path.join(out_dir, 'Snakefile')}")
    
    # Snakemake execution option
    if args.use_snakemake:
        try:
            import snakemake
            logger.info("Triggering workflow execution DAG using Snakemake's Python API...")
            success = snakemake.snakemake(
                os.path.join(out_dir, "Snakefile"),
                config=config_dict
            )
            if success:
                logger.info("[PASS] Snakemake DAG execution finished successfully.")
            else:
                logger.warning("Snakemake execution returned non-success status. Falling back to sequential Python execution.")
        except ImportError:
            logger.warning("[FALLBACK ACTIVE] Snakemake package not found in current environment. Defaulting to Python sequential execution DAG.")

    # ─── Run Stage 0: Assembly Quality Control ───
    qc_result = None
    try:
        qc_result = run_qc_pipeline(args.query, out_dir, logger)
    except Exception as e:
        logger.error(f"[ERROR] Stage 0 (QC) failed: {e}")
    
    # ─── Run Stage 1: Annotation & Comparative Genomics ───
    result_anno = None
    orth_map = None
    ref_strains = []
    identities = {}
    sim_matrix = None
    try:
        result_anno = run_annotation_pipeline(
            args.query, args.ref, args.cores, args.email, out_dir, logger, no_download=args.no_download, qc_result=qc_result
        )
        orth_map = result_anno.prediction.get("ortholog_mapping")
        ref_strains = result_anno.prediction.get("ref_strains", [])
        identities = result_anno.prediction.get("identities", {})
        sim_matrix = result_anno.prediction.get("sim_matrix")
    except Exception as e:
        logger.error(f"[ERROR] Stage 1 (Annotation) failed: {e}")
    
    # ─── Run Stage 2: Metabolic Reconstruction ───
    metab_data = None
    if orth_map is not None:
        try:
            metab_data = run_metabolism_pipeline(args.query, orth_map, out_dir, logger, force_refresh=args.force_kegg_refresh)
        except Exception as e:
            logger.error(f"[ERROR] Stage 2 (Metabolism) failed: {e}")
    else:
        logger.warning("[WARNING] Skipping Stage 2 (Metabolism) due to missing orth_map from Stage 1.")
    
    # ─── Run Stage 3: Application Screening & Phenotypic Profit Catalog ───
    result_screen = None
    screened_cargo = None
    if orth_map is not None:
        try:
            result_screen = run_screening_pipeline(args.query, orth_map, out_dir, logger, qc_result=qc_result)
            screened_cargo = result_screen.prediction
        except Exception as e:
            logger.error(f"[ERROR] Stage 3 (Screening) failed: {e}")
    else:
        logger.warning("[WARNING] Skipping Stage 3 (Screening) due to missing orth_map from Stage 1.")
    
    # ─── Run Stage 4: Biosynthetic Gene Cluster Delineation ───
    result_bgc = None
    bgc_regions = []
    if orth_map is not None:
        try:
            result_bgc = run_bgc_pipeline(args.query, args.ref, orth_map, out_dir, logger, qc_result=qc_result)
            bgc_regions = result_bgc.prediction.get("bgc_list", [])
            
            if not bgc_regions:
                logger.error("[ERROR] No Biosynthetic Gene Clusters detected in query genome!")
            else:
                logger.info(f"Detected {len(bgc_regions)} candidate Biosynthetic Gene Clusters:")
                for bgc in bgc_regions:
                    logger.info(f" - {bgc['BGC_ID']}: {bgc['BGC_Type']} at {bgc['Start_Coord']}-{bgc['End_Coord']} ({bgc['Size_bp']:,} bp, {len(bgc['Hits'])} core genes)")
        except Exception as e:
            logger.error(f"[ERROR] Stage 4 (BGC Delineation) failed: {e}")
    else:
        logger.warning("[WARNING] Skipping Stage 4 (BGC Delineation) due to missing orth_map from Stage 1.")
        
    if not bgc_regions:
        logger.info("Pipeline complete (no BGCs to process downstream). Exiting.")
        sys.exit(0)
        
    # ─── Run Stage 5: BGC Select Loop & BGC-Specific Analyses ───
    while True:
        bgcs_to_process = []
        if args.target_bgc:
            if args.target_bgc.lower() in ("all", "a"):
                bgcs_to_process = bgc_regions
            elif "," in args.target_bgc:
                parts = [p.strip().lower() for p in args.target_bgc.split(",") if p.strip()]
                for part in parts:
                    for bgc in bgc_regions:
                        if bgc["BGC_ID"].lower() == part:
                            bgcs_to_process.append(bgc)
                            break
            else:
                for bgc in bgc_regions:
                    if bgc["BGC_ID"].lower() == args.target_bgc.lower():
                        bgcs_to_process.append(bgc)
                        break
            if not bgcs_to_process:
                logger.error(f"Target BGC '{args.target_bgc}' not found in detected BGCs!")
                sys.exit(1)
        else:
            print("\n=== Detected Biosynthetic Gene Clusters (BGCs) ===")
            for i, bgc in enumerate(bgc_regions):
                print(f"[{i+1}] {bgc['BGC_ID']}: {bgc['BGC_Type']} (Coordinates: {bgc['Start_Coord']}-{bgc['End_Coord']})")
            print("[A] Run All BGCs")
            print("[Q] Quit Pipeline")
            sys.stdout.flush()

            # WSL/docker subprocess calls can corrupt the inherited Windows console stdin handle.
            # Re-open a direct connection to the console input device before reading input.
            try:
                if os.name == "nt":
                    _con = open("CONIN$", "r")
                    choice = _con.readline().strip()
                    _con.close()
                else:
                    choice = input("\n[INPUT_REQUIRED] Select BGC number(s) to analyze (e.g. '1', '1, 2, 3', 'A' for all, or Q to quit): ").strip()
            except OSError:
                choice = input("\n[INPUT_REQUIRED] Select BGC number(s) to analyze (e.g. '1', '1, 2, 3', 'A' for all, or Q to quit): ").strip()
            if choice.lower() == 'q':
                logger.info("User requested exit. Terminating.")
                sys.exit(0)
            
            if choice.lower() in ('a', 'all'):
                bgcs_to_process = bgc_regions
            elif ',' in choice:
                parts = [p.strip() for p in choice.split(',') if p.strip()]
                for part in parts:
                    try:
                        idx = int(part) - 1
                        if 0 <= idx < len(bgc_regions):
                            bgcs_to_process.append(bgc_regions[idx])
                        else:
                            print(f"Warning: Selection {part} is out of range. Skipping.")
                    except ValueError:
                        print(f"Warning: Invalid selection '{part}'. Skipping.")
            else:
                try:
                    idx = int(choice) - 1
                    if 0 <= idx < len(bgc_regions):
                        bgcs_to_process.append(bgc_regions[idx])
                    else:
                        print(f"Invalid selection. Please choose between 1 and {len(bgc_regions)}.")
                        continue
                except ValueError:
                    print("Invalid input. Please enter a number, comma-separated numbers, 'A', or 'Q'.")
                    continue
            
            if not bgcs_to_process:
                print("No valid BGCs selected. Please try again.")
                continue
        
        for selected_bgc in bgcs_to_process:
            bgc_id = selected_bgc["BGC_ID"]
            logger.info(f"Starting downstream analysis for target: {bgc_id} ({selected_bgc['BGC_Type']})")
            
            # 5a. Delineate core and flanking genes layout
            region_genes, flanking_genes = [], []
            try:
                region_genes, flanking_genes, _, _ = evaluate_targeted_bgc(
                    args.query, args.ref, orth_map, selected_bgc, out_dir, logger
                )
            except Exception as e:
                logger.error(f"[ERROR] Stage 5a (BGC Evaluation) failed for {bgc_id}: {e}")
            
            # 5b. Run target BGC promoter, prophage, and conservation BLAST
            result_prom = None
            promoter_records = []
            phage_hits = []
            blast_rows = []
            try:
                result_prom = run_target_bgc_analysis(
                    args.query, selected_bgc, args.run_blast, args.email, out_dir, logger, qc_result=qc_result
                )
                promoter_records = result_prom.prediction.get("promoter_records", [])
                phage_hits = result_prom.prediction.get("phage_hits", [])
                blast_rows = result_prom.prediction.get("blast_rows", [])
            except Exception as e:
                logger.error(f"[ERROR] Stage 5b (Promoter/BLAST) failed for {bgc_id}: {e}")
            
            # 5c. Run Evolutionary HGT signatures (GC/TNF deviations)
            result_evol = None
            hgt_results = {}
            try:
                result_evol = run_evolution_pipeline(args.query, selected_bgc, out_dir, logger, ortholog_mapping=orth_map, qc_result=qc_result)
                hgt_results = result_evol.prediction if result_evol else {}
            except Exception as e:
                logger.error(f"[ERROR] Stage 5c (HGT signatures) failed for {bgc_id}: {e}")
            
            # 5d. Molecular Clocking & Phylogenetics of top 25 reference strains
            result_phy = None
            try:
                result_phy = run_phylogeny_pipeline(args.query, ref_strains, identities, sim_matrix, out_dir, logger, bgc_blast_results=blast_rows, qc_result=qc_result)
            except Exception as e:
                logger.error(f"[ERROR] Stage 5d (Phylogeny) failed for {bgc_id}: {e}")
            
            # 5e. Database-driven molecular docking coordinate download, ChEMBL Ki lookup & PyMOL subprocess (Speculative & Optional)
            result_dock = None
            from redolarium.structures import PredictionResult
            
            if args.run_docking:
                try:
                    result_dock = run_docking_pipeline(args.query, selected_bgc, region_genes, out_dir, logger, qc_result=qc_result)
                except Exception as e:
                    logger.error(f"[ERROR] Stage 5e (Molecular Docking) failed for {bgc_id}: {e}")
            else:
                logger.info("Molecular docking skipped. (Note: Molecular docking is speculative, experimental, and under active development.)")
            
            if not result_dock:
                result_dock = PredictionResult(
                    prediction=None,
                    confidence_score=0.0,
                    algorithm="Smina / Vinardo (Skipped/Failed)",
                    algorithm_version="N/A",
                    evidence=["Docking skipped or failed."],
                    limitations=["Molecular docking was bypassed or crashed."],
                    runtime=0.0
                )
            
            # 5f. Precursor stoichiometric ATP synthesis cost calculation
            stoichiometry_data = {}
            try:
                stoichiometry_data = run_linkage_pipeline(args.query, region_genes, selected_bgc, out_dir, logger)
            except Exception as e:
                logger.error(f"[ERROR] Stage 5f (Stoichiometry) failed for {bgc_id}: {e}")
            
            # ─── Synthesize Evidence Graph DAG ───
            from redolarium.utils import build_and_draw_evidence_graph
            try:
                build_and_draw_evidence_graph(
                    qc_result, result_anno, result_bgc, result_evol,
                    result_prom, result_dock, result_phy, out_dir, logger
                )
            except Exception as e:
                logger.error(f"[ERROR] Evidence Graph synthesis failed for {bgc_id}: {e}")
            
            # ─── Compile 12-sheet Excel Workbook ───
            xls_filename = f"{bgc_id}_BGC_Analysis_Metabolism_Integrated.xlsx"
            xls_path = os.path.join(out_dir, xls_filename)
            query_record = SeqIO.read(args.query, "genbank")
            query_org = query_record.annotations.get("organism", "Query Isolate")
            
            logger.info(f"Compiling 12-sheet Excel workbook for {bgc_id}...")
            compile_bgc_excel_report(
                selected_bgc, orth_map, ref_strains, metab_data,
                region_genes, flanking_genes, promoter_records, phage_hits,
                hgt_results, stoichiometry_data, xls_path, query_org, logger,
                bgc_blast_results=blast_rows,
                qc_result=qc_result,
                result_anno=result_anno,
                result_bgc=result_bgc,
                result_evol=result_evol,
                result_prom=result_prom,
                result_dock=result_dock,
                result_phy=result_phy
            )
            
            # ─── Compile Word Report ───
            doc_filename = f"{bgc_id}_BGC_Report.docx"
            logger.info(f"Drafting comprehensive Word report for {bgc_id}...")
            
            # Concat uncertainties (limitations) across all nodes safely
            all_unc = []
            for res in [qc_result, result_anno, result_bgc, result_evol, result_prom, result_dock, result_phy]:
                if res and getattr(res, 'limitations', None):
                    all_unc.extend(res.limitations)
            unc_paras = [f"- {u}" for u in all_unc]
            
            # Check for low-confidence modules requiring manual review
            flagged_modules = []
            for name, res in [("Assembly Quality Control", qc_result), ("BGC Delineation", result_bgc), ("Molecular Docking Validation", result_dock)]:
                if res and getattr(res, "manual_intervention_required", False):
                    flagged_modules.append(name)
            
            sections = [
                ("1. Executive Summary", [
                    f"This document presents the complete genomic, metabolic, and secondary metabolism Biosynthetic "
                    f"Gene Cluster (BGC) analysis reconstructed for cluster {bgc_id} in {query_org}. The automated "
                    f"workflow evaluates sequence similarity, stoichiometric precursors, molecular docking, and "
                    f"horizontal gene transfer signatures to evaluate biological potential and evolutionary context."
                ]),
                ("2. Biosynthetic Gene Cluster Core Architecture", [
                    f"The core region of {bgc_id} spans {selected_bgc['Core_Start']:,} to {selected_bgc['Core_End']:,} bp "
                    f"and contains {len(region_genes)} core genes. Ortholog mapping indicates key biosynthetic roles "
                    f"and sequence alignments are detailed in sheet 'BGC_Gene_Architecture' of the Excel workbook."
                ]),
                ("3. Promoter & Transcription Regulation", [
                    f"Upstream regions of core genes were scanned for consensus Sigma boxes and Shine-Dalgarno motifs. "
                    f"The promoter motif alignments suggest stress-associated transcriptional coordination. The detailed "
                    f"promoter motif sequences and positions are saved in results under the 'bgc_motifs' directory."
                ]),
                ("4. Flanking Prophage Remnants & MGE Screening", [
                    (
                        f"No mobile genetic elements were detected in the 30 kb flanking coordinates of {bgc_id}."
                        if len(phage_hits) == 0 else
                        f"Scanning of the 30 kb flanking coordinates revealed {len(phage_hits)} mobile genetic element(s) "
                        f"(integrase, transposase, or prophage-related genes). "
                        f"This is consistent with, but does not conclusively demonstrate, horizontal gene transfer. "
                        f"Experimental validation (e.g., comparative genomics, transcriptional analysis) is required "
                        f"to confirm the acquisition mechanism of this BGC."
                    ),
                    f"The flanking genetic environment is mapped in the Excel sheet 'Phage_Artifacts'."
                ]),
                ("5. Horizontal Gene Transfer (HGT) Signatures", [
                    f"GC content deviations and tetranucleotide frequency (TNF) distances trace the horizontal transfer history of the BGC. "
                    f"Sliding-window genomic anomalies are plotted in results/hgt_evolution/{bgc_id}_hgt_signatures.png.",
                    f"Detailed window Z-scores and composite HGT scores are saved in the Excel sheet 'GC_Profile'."
                ]),
                ("6. Molecular Clocking & Divergence Timeline", [
                    f"Phylogenetic placement and divergence times are modeled relative to baseline reference nodes. "
                    f"The divergence timeline is plotted in results/phylogeny_trees/divergence_timeline.png.",
                    f"IMPORTANT: Divergence time estimates use the Ochman et al. 2000 bacterial molecular clock "
                    f"(1% divergence per 15 Mya). This clock has ~50% uncertainty without fossil calibration "
                    f"and was calibrated for Enterobacteria. Treat as approximate only.",
                    f"Event markers are saved in results/tabular_data/divergence_events.csv."
                ]),
                ("7. Structural Molecular Docking Validation", [
                    f"Molecular docking of structural peptides against target receptor molecules was conducted. A PyMOL rendering PML script "
                    f"and 2D contact distance map are saved under results/docking_images/.",
                    f"Binding energies and contact residue lists are exported to the Excel sheet 'BLAST_Homologs'."
                ]),
                ("8. Precursor Stoichiometry and Metabolic Linkage", [
                    f"Ribosomal translation ATP demand and tRNA charging stoichiometry were calculated. The metabolic precursor flow "
                    f"and energy pool linkage network are plotted in results/metabolic_pathways/{bgc_id}_metabolic_linkage.png.",
                    f"Stoichiometric details and ATP costs are saved in the Excel sheet 'Metabolic_Integration'."
                ]),
                ("9. BGC Cross-Species Conservation", [
                    f"Cross-species conservation analysis of the core biosynthetic sequence of {bgc_id} was performed using remote NCBI BLAST. "
                    f"The alignment results, including E-values and sequence identities across homologous bacterial genomes, are exported to "
                    f"results/tabular_data/{bgc_id}_conservation_blast.csv."
                ]),
                ("10. Systemic Uncertainties & Limitations", [
                    "The following limitations and analytical caveats were identified across the execution checkpoints. "
                    "These uncertainties should be reported honestly in peer review:"
                ] + unc_paras)
            ]
            
            if flagged_modules:
                sections.insert(0, ("[FLAGGED FOR REVIEW] Expert Action Required", [
                    f"WARNING: The pipeline has detected low confidence scores (below 0.3) in the following key modules: {', '.join(flagged_modules)}.",
                    "An expert must manually inspect results/tabular_data/evidence_graph.json to audit and verify these findings before publication."
                ]))
                
            write_word_report(doc_filename, bgc_id, query_org, sections, out_dir)
            logger.info(f"Word report compiled successfully: {os.path.join(out_dir, doc_filename)}")
            
        if args.target_bgc:
            logger.info("Target BGC finished. Exiting pipeline run.")
            break
            
        again = input("\nDo you want to analyze another BGC or selection? (y/n): ").strip().lower()
        if again != 'y':
            logger.info("Pipeline complete. Exiting.")
            break

if __name__ == "__main__":
    main()
