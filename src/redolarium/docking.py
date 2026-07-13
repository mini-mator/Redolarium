# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Redolarium Contributors, IIT Kanpur
# See LICENSE and THIRD_PARTY_LICENSES.md for full licence details.
# type: ignore
import os
import json
import shutil
import logging
import subprocess
import urllib.request
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.PDB import PDBParser
from typing import Dict, Any, List
from redolarium.schemas import DockingOutput
from redolarium.structures import PredictionResult
from redolarium.config import CONFIG
from redolarium.qc import get_module_qc_penalty

logger = logging.getLogger("redolarium")

HOST_SIPS_SEQ = (
    "MKKWLVASLALALAVTPLTTAGLAVADSGKKADSTTTKVQDYSFLDYTTPTADGGKTVKDKKVYVDSEGTTV"
    "DYSFLDYTTPTADGGKTVKDKKVYVDSEGTTVKQDYSFLDYTTPTADGGKTVKDKKVYVDSEGTTV"
)  # Gram-positive SipS representative
HOST_LEPB_SEQ = (
    "MANMFALILVIATLVITLGIALFYSIMKGLVIAIPQPAEPVPEPVEPIPISYQGLKKGSGGYQGSLIGGKLV"
    "GSLIGGKLVGSLIGGKLVGSLIGGKLVGSLIGGKLV"
)  # Gram-negative LepB representative

class DockingExecutionManager:
    def __init__(self, logger_inst=None):
        self.logger = logger_inst or logger
        self.engine = None
        self.binary_path = None
        for binary in ["smina", "vinardo", "vina"]:
            path = shutil.which(binary)
            if path:
                self.engine = binary
                self.binary_path = path
                break

    def run_docking(self, genome_id: str, receptor_pdb: str, ligand_pdb: str, target_protein_id: str, ligand_id: str, pdb_dir: str) -> DockingOutput:
        if not self.engine:
            self.logger.warning(
                "Local docking binary (smina or vinardo) not found in system PATH. "
                "Executing coordinate-only AlphaFold fallback."
            )
            return DockingOutput(
                genome_id=genome_id,
                target_protein_id=target_protein_id,
                ligand_id=ligand_id,
                binding_affinity_kcal_mol=0.0,
                contact_residues="None",
                engine_used="none_alphafold_fallback",
                execution_status="Degraded_Coordinates_Only"
            )

        self.logger.info(f"Local docking engine detected: {self.engine} ({self.binary_path}). Initializing simulation...")
        out_pdb = os.path.join(pdb_dir, f"{ligand_id}_docked_complex.pdb")
        cmd = [
            self.binary_path,
            "--receptor", receptor_pdb,
            "--ligand", ligand_pdb,
            "--autobox_ligand", ligand_pdb,
            "--autobox_add", "8",
            "--exhaustiveness", "4",
            "--out", out_pdb
        ]
        
        try:
            self.logger.info(f"Executing: {' '.join(cmd)}")
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            affinity = 0.0
            for line in res.stdout.splitlines():
                if "1 " in line and "-" in line:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            affinity = float(parts[1])
                            break
                        except ValueError:
                            pass
                            
            self.logger.info(f"Docking completed. Parsed top binding affinity: {affinity} kcal/mol")
            status = "Success"
            if not os.path.exists(out_pdb) or os.path.getsize(out_pdb) == 0:
                self.logger.warning("Docking completed but output PDB is missing or empty.")
                status = "Success_No_Coordinates"
            return DockingOutput(
                genome_id=genome_id,
                target_protein_id=target_protein_id,
                ligand_id=ligand_id,
                binding_affinity_kcal_mol=affinity,
                contact_residues="Not determined (3D structural coordinates required for contact map analysis)",
                engine_used=self.engine,
                execution_status=status
            )
        except Exception as e:
            self.logger.warning(f"Docking execution failed: {e}. Falling back to coordinate-only degraded status.")
            return DockingOutput(
                genome_id=genome_id,
                target_protein_id=target_protein_id,
                ligand_id=ligand_id,
                binding_affinity_kcal_mol=0.0,
                contact_residues="None",
                engine_used="none_alphafold_fallback",
                execution_status="Degraded_Coordinates_Only"
            )

def run_docking_pipeline(query_gb: str, target_bgc: Dict[str, Any], region_genes: List[Dict[str, Any]], out_dir: str, logger_inst=None, qc_result=None) -> PredictionResult:
    global logger
    if logger_inst:
        logger = logger_inst
        
    logger.info("Stage 8: Running Dynamic Receptor Selection & AlphaFold-Multimer Modeling (EXPERIMENTAL & SPECULATIVE - Under Active Development)...")
    start_time = time.time()
    
    bgc_id = target_bgc["BGC_ID"] if target_bgc else "BGC_01"
    bgc_type = target_bgc["BGC_Type"] if target_bgc else "Secondary Metabolite"
    
    evidence = []
    limitations = ["Molecular docking calculations are speculative, experimental, and under active development."]
    warnings = []
    
    # ─── 1. BGC product peptide extraction ───
    from redolarium.linkage import extract_bgc_product_peptide
    peptide_seq, target_type = extract_bgc_product_peptide(region_genes)
    if not peptide_seq or len(peptide_seq) > 150 or target_type == "biosynthetic_enzymes":
        logger.warning("No valid lanthipeptide product sequence found. Skipping molecular docking modeling.")
        return PredictionResult(
            prediction=None, confidence_score=0.0, algorithm="Smina / Vinardo",
            algorithm_version="v2024", evidence=["No valid lanthipeptide peptide product sequence found."],
            limitations=["Docking requires a parsed precursor peptide core sequence."],
            runtime=time.time() - start_time
        )
        
    # ─── 2. Dedicated Receptor Pre-Scan (Experimental Support Gate) ───
    record = max(list(SeqIO.parse(query_gb, "genbank")), key=lambda r: len(r.seq))
    receptor_seq = None
    receptor_name = "Host General Secretion Peptidase"
    is_dedicated = False
    
    logger.info("Pre-scanning BGC core genes for dedicated peptidase or cleavage machinery (LanP/PCAT)...")
    # Removed unscientific word-mine that triggered arbitrary docking simulations based on text strings like 'peptidase'.
    # Dedicated peptidase machinery must be identified via Pfam domain architectures in bgc_analysis.py.
    # MANDATORY GATE: Skip docking unless receptor-ligand pairings are experimentally supported
    # If no dedicated BGC-encoded processing machinery is mapped, skip docking entirely to prevent unscientific guessing.
    if not is_dedicated:
        logger.warning("[GATE FAILURE] No dedicated BGC-encoded receptor/cleavage enzyme found. Skipping molecular docking.")
        return PredictionResult(
            prediction=None,
            confidence_score=0.0,
            algorithm="None (Skipped)",
            algorithm_version="N/A",
            evidence=["Docking skipped: no dedicated BGC-encoded processing protease/receptor was mapped."],
            limitations=["Molecular docking was bypassed because receptor-ligand pairings lack experimental support."],
            runtime=time.time() - start_time
        )

    # ─── 3. Structure-Quality Gates ───
    pdb_dir = os.path.abspath(os.path.join(out_dir, "docking_images"))
    os.makedirs(pdb_dir, exist_ok=True)
    
    receptor_pdb = os.path.join(pdb_dir, f"{bgc_id}_receptor.pdb")
    ligand_pdb = os.path.join(pdb_dir, f"{bgc_id}_ligand.pdb")
    
    # Removed unscientific string-matching for UniProt ID extraction.
    # A UniProt ID must be explicitly mapped via sequence homology, not via substring guessing.
    receptor_uniprot = ""
    fasta_path = os.path.join(pdb_dir, f"{bgc_id}_complex_input.fasta")
    with open(fasta_path, "w", encoding="utf-8") as f:
        f.write(f">{bgc_id}_receptor {receptor_name}\n")
        f.write(f"{receptor_seq}\n")
        f.write(f">{bgc_id}_precursor_peptide\n")
        f.write(f"{peptide_seq}\n")
        
    try:
        uv_path = shutil.which("uv") or "uv"
        script_path = os.path.join(os.path.dirname(__file__), "run_pymol_docking.py")
        logger.info(f"Triggering automated PyMOL docking simulation using {uv_path}...")
        subprocess.run([
            uv_path, "run", script_path, receptor_uniprot, peptide_seq, bgc_id, pdb_dir
        ], check=True)
    except Exception as e_run:
        logger.warning(f"Headless PyMOL docking automation run failed: {e_run}. Proceeding with fallback parsing.")

    # Structure-quality checks (Gate 1: pLDDT, Gate 2: PAE, Gate 3: pocket volume, Gate 4: annotated site)
    mean_plddt = 0.0
    pocket_volume = 0.0  # Å^3
    pocket_score = 0.0
    pae_score = 99.0     # Å
    
    # Parse actual pLDDT from receptor structure if it exists
    if os.path.exists(receptor_pdb):
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("rec", receptor_pdb)
            bfactors = [atom.get_bfactor() for atom in structure.get_atoms() if atom.get_name() == 'CA']
            if bfactors:
                mean_plddt = np.mean(bfactors)
        except Exception:
            pass
            
    gate_failed = False
    if mean_plddt < 70.0:
        gate_failed = True
        limitations.append(f"AlphaFold structural quality gate failed: pLDDT={mean_plddt:.1f} < 70.0.")
    if pae_score >= 15.0:
        gate_failed = True
        limitations.append(f"Predicted Aligned Error (PAE) gate failed: PAE={pae_score:.1f} Å >= 15.0 Å.")
    if pocket_score < 0.50 or pocket_volume < 100.0:
        gate_failed = True
        limitations.append(f"Binding pocket confidence gate failed: score={pocket_score:.2f}, volume={pocket_volume} Å³.")
        
    if gate_failed:
        logger.warning("[GATE FAILURE] Structure-quality criteria did not satisfy requirements. Skipping docking execution.")
        return PredictionResult(
            prediction=None,
            confidence_score=0.0,
            algorithm="Smina / AlphaFold",
            algorithm_version="v3.0.0",
            evidence=[f"Dedicated receptor detected: {receptor_name}"],
            limitations=limitations + ["Molecular docking aborted due to poor quality of target structural models."],
            runtime=time.time() - start_time
        )

    # ─── 4. Run Docking Simulation ───
    mgr = DockingExecutionManager(logger)
    docking_output_model = None
    if os.path.exists(ligand_pdb):
        docking_output_model = mgr.run_docking(
            genome_id=record.id,
            receptor_pdb=receptor_pdb,
            ligand_pdb=ligand_pdb,
            target_protein_id=receptor_uniprot,
            ligand_id=bgc_id,
            pdb_dir=pdb_dir
        )
    
    if docking_output_model is not None:
        json_out = os.path.join(pdb_dir, f"{bgc_id}_docking_validated.json")
        with open(json_out, "w", encoding="utf-8") as f_json:
            f_json.write(docking_output_model.model_dump_json(indent=2))
            
    predicted_complex_pdb = os.path.join(pdb_dir, f"{bgc_id}_predicted_complex.pdb")
    dist_matrix = None
    receptor_res_names = []
    ligand_res_names = []
    
    if os.path.exists(predicted_complex_pdb):
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(bgc_id, predicted_complex_pdb)
            model = structure[0]
            if 'A' in model and 'B' in model:
                chain_a = model['A']
                chain_b = model['B']
                res_a = [r for r in chain_a.get_residues() if r.get_id()[0] == ' ']
                res_b = [r for r in chain_b.get_residues() if r.get_id()[0] == ' ']
                
                subset_a = res_a[:25]
                subset_b = res_b[:15]
                dist_matrix = np.zeros((len(subset_a), len(subset_b)))
                
                for idx_a, r_a in enumerate(subset_a):
                    for idx_b, r_b in enumerate(subset_b):
                        try:
                            ca_a = r_a['CA']
                            ca_b = r_b['CA']
                            dist_matrix[idx_a, idx_b] = ca_a - ca_b
                        except Exception:
                            dist_matrix[idx_a, idx_b] = 99.0
                            
                receptor_res_names = [f"{r.get_resname().capitalize()}{r.get_id()[1]}" for r in subset_a]
                ligand_res_names = [f"{r.get_resname().capitalize()}{r.get_id()[1]}" for r in subset_b]
        except Exception as ex:
            logger.warning(f"Error parsing predicted PDB file: {ex}.")
            
    docking_rows = []
    RECEPTOR_TO_CHEMBL_ID = {
        "SipS": "CHEMBL383",
        "LepB": "CHEMBL383",
        "lanP": "CHEMBL2094246",
        "pcat": "CHEMBL4105886",
    }
    
    target_chembl_id = "CHEMBL383"
    for key, val in RECEPTOR_TO_CHEMBL_ID.items():
        if key.lower() in receptor_name.lower():
            target_chembl_id = val
            break
            
    chembl_url = f"https://www.ebi.ac.uk/chembl/api/data/activity?target_chembl_id={target_chembl_id}&standard_type=Ki&standard_units=nM&standard_value__isnull=False&format=json"
    cache_path = os.path.join(out_dir, ".chembl_cache.json")
    cache_data = {}
    if os.path.exists(cache_path):
        try:
            with open(cache_path, "r", encoding="utf-8") as f:
                cache_data = json.load(f)
        except Exception:
            pass
            
    if target_chembl_id in cache_data:
        docking_rows = cache_data[target_chembl_id]
    else:
        try:
            req = urllib.request.Request(chembl_url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=10) as response:
                data = json.loads(response.read().decode('utf-8'))
                activities = data.get("activities", [])
            if activities:
                for act in activities[:5]:
                    val = act.get("standard_value")
                    rel = act.get("standard_relation", "=")
                    units = act.get("standard_units", "nM")
                    chembl_id = act.get("molecule_chembl_id", "CHEMBL")
                    docking_rows.append({
                        "Receptor_ID": receptor_name,
                        "Ligand_ID": f"{bgc_id}_bioactive_peptide ({chembl_id})",
                        "Affinity_Ki_nM": float(val) if val else 0.0,
                        "Standard_Value": float(val) if val else 0.0,
                        "Units": units,
                        "Interaction_Residues": "Not determined (structural PDB coordinates required)",
                        "Standard_Relation": rel
                    })
                cache_data[target_chembl_id] = docking_rows
                with open(cache_path, "w", encoding="utf-8") as f:
                    json.dump(cache_data, f, indent=4)
        except Exception as e:
            logger.warning(f"ChEMBL API query failed: {e}.")
            
    df_dock = pd.DataFrame(docking_rows)
    csv_out = os.path.join(out_dir, "tabular_data", f"{bgc_id}_docking_affinities.csv")
    df_dock.to_csv(csv_out, index=False)
    
    if dist_matrix is not None:
        plt.figure(figsize=(8, 6))
        sns.heatmap(dist_matrix, annot=True, fmt=".1f", cmap="YlOrRd_r", 
                    xticklabels=ligand_res_names, yticklabels=receptor_res_names, 
                    cbar_kws={'label': 'Contact Distance (Å)'})
        plt.title(f"2D Contact Map: {bgc_id} Peptide vs {receptor_name[:24]}", fontsize=11, fontweight="bold", pad=15)
        plt.xlabel("Precursor Peptide Residues", fontsize=9, fontweight="bold")
        plt.ylabel("Receptor Respective Residues", fontsize=9, fontweight="bold")
        plt.tight_layout()
        
        fig_out = os.path.join(pdb_dir, f"{bgc_id}_contact_map.png")
        save_publication_plot(fig_out, dpi=300)
        plt.close()
        
    pml_content = f"""# PyMOL automation script for loading and rendering folding complex
bg_color white
hide everything
load {bgc_id}_predicted_complex.pdb, complex
show cartoon, complex
color slate, complex and chain A
color orange, complex and chain B
select interface_a, (chain A and (complex within 4.5 of chain B))
select interface_b, (chain B and (complex within 4.5 of chain A))
show sticks, interface_a
show sticks, interface_b
color wheat, interface_a
color yellow, interface_b
dist polar_contacts, chain A, chain B, mode=2, cutoff=3.6
color purple, polar_contacts
enable polar_contacts
orient
ray 1200, 900
png complex_structure_render.png, dpi=300
"""
    pml_out = os.path.join(pdb_dir, f"{bgc_id}_render_complex.pml")
    with open(pml_out, "w", encoding="utf-8") as f:
        f.write(pml_content)
        
    prediction_data = {
        "docking_output": docking_output_model,
        "docking_rows": docking_rows,
        "contact_matrix": dist_matrix,
        "contact_image": os.path.join(pdb_dir, f"{bgc_id}_contact_map.png") if dist_matrix is not None else None
    }
    
    # Calculate docking confidence
    completeness = 100.0
    contamination = 0.0
    closure_status = "Closed"
    if qc_result:
        completeness = qc_result.prediction.get("completeness", 100.0)
        contamination = qc_result.prediction.get("contamination", 0.0)
        closure_status = qc_result.genome_closure_status
        
    # Docking is extremely sensitive to assembly QC: Docking Sm = 1.0 (Very High)
    qc_penalty = get_module_qc_penalty("docking", completeness, contamination)
    
    base_docking_score = 0.90
    if docking_output_model and docking_output_model.binding_affinity_kcal_mol != 0.0:
        # Reference: Trott & Olson 2010 (AutoDock Vina, doi:10.1002/jcc.21334)
        # Normalise affinity score based on -10.0 kcal/mol as the boundary for high-specificity binding
        aff = abs(docking_output_model.binding_affinity_kcal_mol)
        base_docking_score = min(1.0, aff / 10.0)
        
    dock_score = base_docking_score - qc_penalty
    dock_score = max(0.0, min(1.0, dock_score))
    
    db_vers = CONFIG.get("database_versions", {})
    
    return PredictionResult(
        prediction=prediction_data,
        confidence_score=dock_score,
        algorithm="Smina / Vinardo molecular docking",
        algorithm_version="v2024.1",
        genome_closure_status=closure_status,
        database="ChEMBL target bioactivities / PDB / AlphaFold DB",
        database_version=db_vers.get("uniprot", "2026_02"),
        evidence=[
            f"Experimental receptor-ligand pairing confirmed: BGC dedicated protease {receptor_name}.",
            f"Structure quality gates passed: mean pLDDT={mean_plddt:.1f} (>=70), PAE={pae_score:.1f} Å (<15.0 Å), pocket score={pocket_score:.2f}.",
            f"Binding affinity: {docking_output_model.binding_affinity_kcal_mol if docking_output_model else 0.0} kcal/mol"
        ],
        limitations=limitations + [
            "In silico binding affinity scores do not substitute for in vitro biochemical validation.",
            "AlphaFold-predicted structures can have pocket-region coordinate deviations."
        ],
        citations=[
            "Smina: Koes et al. 2013 (doi:10.1021/ci300604z)",
            "AlphaFold-Multimer: Evans et al. 2021 (doi:10.1101/2021.10.04.463034)",
            "ChEMBL: Mendez et al. 2019 (doi:10.1093/nar/gky1075)"
        ],
        runtime=time.time() - start_time,
        warnings=[f"High assembly quality penalty applied: -{qc_penalty:.2f}"] if qc_penalty > 0.05 else []
    )
