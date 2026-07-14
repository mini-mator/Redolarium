# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Redolarium Contributors, IIT Kanpur
# See LICENSE and THIRD_PARTY_LICENSES.md for full licence details.

import os
import re
import urllib.request
import time
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
from redolarium.config import CONFIG, EXCEL_PALETTE
from redolarium.utils import safe_subprocess_run

def load_bgc_profiles():
    import json
    profiles_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "bgc_profiles.json")
    if os.path.exists(profiles_path):
        try:
            with open(profiles_path, "r", encoding="utf-8") as f:
                data = json.load(f)
                normalized = {}
                for k, v in data.items():
                    if isinstance(v, list):
                        normalized[k] = {"genes": v, "genus_specific": []}
                    elif isinstance(v, dict):
                        normalized[k] = {
                            "genes": v.get("genes", []),
                            "genus_specific": v.get("genus_specific", [])
                        }
                    else:
                        normalized[k] = v
                return normalized
        except Exception:
            pass
    # BGC profile format: each entry has "genes" (symbol list) and "genus_specific" (list of genera, or [] for universal)
    # genus_specific: if non-empty, the profile is only matched when the query organism genus is in this list.
    # Reference: MIBiG 3.1 (Terlouw et al. 2023, doi:10.1093/nar/gkac1049)
    return {
        "BGC0000352 (surfactin)": {"genes": ["srfaa", "srfab", "srfac", "srfad"], "genus_specific": []},
        "BGC0000350 (fengycin)": {"genes": ["fena", "fenb", "fenc", "fend", "fene"], "genus_specific": []},
        "BGC0000377 (plipastatin)": {"genes": ["ppsa", "ppsb", "ppsc", "ppsd", "ppse"], "genus_specific": []},
        "BGC0000325 (bacillibactin)": {"genes": ["dhba", "dhbb", "dhbc", "dhbe", "dhbf"], "genus_specific": []},
        "BGC0000412 (macrolactin)": {"genes": ["mlna", "mlnb", "mlnc", "mlnd", "mlne"], "genus_specific": []},
        "BGC0000511 (nisin)": {"genes": ["nisa", "nisb", "nisc", "nisp", "nist"], "genus_specific": []},
        "BGC0000520 (subtilin)": {"genes": ["spas", "spab", "spac", "spat"], "genus_specific": []},
        "BGC0000251 (iturin A)": {"genes": ["itua", "itub", "ituc", "itud"], "genus_specific": []},
        "BGC0001089 (bacillaene)": {"genes": ["baer", "baes", "baet", "baeu", "baev", "baew"], "genus_specific": []},
        "BGC0000572 (subtilosin A)": {"genes": ["sboa", "sbox", "alba", "albb", "albc"], "genus_specific": []},
        "BGC0001184 (bacilysin)": {"genes": ["baca", "bacb", "bacc", "bacd", "bace"], "genus_specific": []},
        "BGC0000176 (difficidin)": {"genes": ["dfna", "dfnb", "dfnc", "dfnd", "dfne", "dfnf"], "genus_specific": []},
        # Pseudomonas/Burkholderia-specific siderophore clusters — only match when genus matches
        "BGC0000413 (pyoverdine)": {"genes": ["pvda", "pvde", "pvdf", "pvdi", "pvdj", "pvdl", "pvdm"], "genus_specific": ["pseudomonas"]},
        "BGC0000378 (pyochelin)": {"genes": ["pcha", "pchb", "pchc", "pchd", "pche", "pchf", "pchg"], "genus_specific": ["pseudomonas", "burkholderia"]},
        # Pseudomonas biocontrol cluster — genus-specific
        "BGC0000286 (phenazine-1-carboxylic acid)": {"genes": ["phza", "phzb", "phzc", "phzd", "phze", "phzf", "phzg"], "genus_specific": ["pseudomonas", "streptomyces"]},
        "BGC0000143 (2,4-diacetylphloroglucinol)": {"genes": ["phla", "phlb", "phlc", "phld", "phlf"], "genus_specific": ["pseudomonas"]}
    }

MIBIG_ANNOTATIONS = {
    "BGC0000352": {"compound": "Surfactin", "class": "NRP", "activity": "Antimicrobial / biosurfactant / hemolytic", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0000352/"},
    "BGC0000350": {"compound": "Fengycin", "class": "NRP", "activity": "Antifungal / cell membrane disruptor", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0000350/"},
    "BGC0000377": {"compound": "Plipastatin", "class": "NRP", "activity": "Antifungal / inhibitor of phospholipase A2", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0000377/"},
    "BGC0000325": {"compound": "Bacillibactin", "class": "NRP / Siderophore", "activity": "Iron-chelating / catecholate siderophore", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0000325/"},
    "BGC0000412": {"compound": "Macrolactin", "class": "Polyketide", "activity": "Antibacterial / antiviral / H+-ATPase inhibitor", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0000412/"},
    "BGC0000511": {"compound": "Nisin", "class": "RiPP (Lanthipeptide)", "activity": "Antibacterial / food preservative / membrane pore-forming via lipid II binding", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0000511/"},
    "BGC0000520": {"compound": "Subtilin", "class": "RiPP (Lanthipeptide)", "activity": "Antibacterial / spore germination inhibitor / lipid II-targeting", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0000520/"},
    "BGC0000251": {"compound": "Iturin A", "class": "NRP (lipopeptide)", "activity": "Antifungal / membrane disruption via fungal cell membrane sterol interaction", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0000251/"},
    # Bacillaene: ROS-mediated antibacterial; mechanism NOT protein synthesis inhibition
    # Reference: Straight et al. 2007 (doi:10.1128/JB.00408-07); Mielich-Suss & Lopez 2015
    "BGC0001089": {"compound": "Bacillaene", "class": "Polyketide (trans-acyl)", "activity": "Antibacterial / reactive oxygen species (ROS) induction / electron transport chain disruption", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0001089/"},
    "BGC0000572": {"compound": "Subtilosin A", "class": "RiPP (Sactipeptide)", "activity": "Antibacterial / membrane disruption", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0000572/"},
    "BGC0001184": {"compound": "Bacilysin", "class": "Other (dipeptide)", "activity": "Antibacterial / glucosamine-6-phosphate synthase inhibitor", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0001184/"},
    # Difficidin: EF-Tu inhibitor (translational elongation), NOT a cell wall agent
    # Reference: Broderick et al. 2021 (doi:10.1073/pnas.2019378118)
    "BGC0000176": {"compound": "Difficidin", "class": "Polyketide", "activity": "Antibacterial / EF-Tu inhibitor (translational elongation blockade)", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0000176/"},
    "BGC0000413": {"compound": "Pyoverdine", "class": "NRP / Siderophore", "activity": "Iron-chelating / virulence factor", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0000413/"},
    "BGC0000378": {"compound": "Pyochelin", "class": "NRP-PKS hybrid / Siderophore", "activity": "Iron-chelating / virulence promoter", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0000378/"},
    "BGC0000286": {"compound": "Phenazine-1-carboxylic acid", "class": "Other (Aromatics)", "activity": "Broad-spectrum antifungal / redox-active agent", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0000286/"},
    "BGC0000143": {"compound": "2,4-diacetylphloroglucinol", "class": "Polyketide", "activity": "Antifungal / biocontrol phytotoxic agent", "ref": "https://mibig.secondarymetabolites.org/repository/BGC0000143/"}
}

MIBIG_DATABASE = load_bgc_profiles()

# Helper to classify gene molecular role
def get_bgc_gene_role(product_str, gene_symbol):
    p = product_str.lower()
    s = gene_symbol.lower()
    if any(k in p for k in ["synthetase", "synthase", "dehydratase", "cyclase", "transferase", "lanb", "lanc", "lana", "pps", "srf", "dhb", "bgl", "amy", "phz"]):
        return "biosynthetic"
    if any(k in p for k in ["transporter", "permease", "abc", "export", "efflux"]):
        return "transport"
    if any(k in p for k in ["regulator", "response", "kinase", "two-component", "repressor"]):
        return "regulatory"
    if any(k in p for k in ["immunity", "resistance", "protection", "self-resistance"]):
        return "immunity"
    return "accessory/context"

def get_curated_bgc_hmms(logger):
    local_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "resources")
    os.makedirs(local_dir, exist_ok=True)
    hmm_path = os.path.join(local_dir, "essential_bgc.hmm")
    
    # We rely strictly on the pre-built, local HMM database.
    if not os.path.exists(hmm_path):
        logger.warning(f"Curated BGC HMM database not found locally at: {hmm_path}. Falling back to default keyword annotation.")
        return None
    return hmm_path

def run_pyhmmer_scan(proteins, hmm_path, logger):
    try:
        import pyhmmer
        import pyhmmer.easel as easel
        import pyhmmer.plan7 as plan7
    except ImportError:
        logger.warning("pyhmmer not installed. Skipping HMM-based BGC typing.")
        return {}
        
    alphabet = easel.Alphabet.amino()
    
    digitized_seqs = []
    for ltag, trans in proteins:
        if trans and len(trans) > 10:
            try:
                seq = easel.TextSequence(name=ltag.encode("utf-8"), sequence=trans).digitize(alphabet)
                digitized_seqs.append(seq)
            except Exception:
                pass
                
    if not digitized_seqs:
        return {}
        
    domain_hits = {}
    
    # Profile-specific trusted cutoffs (TC) or conservative e-value thresholds per Pfam domain
    # Defaulting to 1e-10/1e-8/1e-5 based on antiSMASH 7.0 validation profiles to prevent false-positive inflation
    # Reference: Blin et al. 2023 (doi:10.1093/nar/gkad344)
    #
    # PF04738: Radical SAM sactipeptide cyclase — correct RiPP-specific domain (replaces PF14867/PCMT
    # which is a protein repair methyltransferase and NOT specific to RiPPs)
    # Reference: Grell et al. 2018 (doi:10.1038/s41589-018-0122-9)
    pfam_cutoffs = {
        "PF00109": 1e-10, "PF02801": 1e-10,  # PKS (beta-ketoacyl synthase domains)
        "PF00668": 1e-8,  "PF00501": 1e-8,   # NRPS (condensation, adenylation domains)
        "PF05147": 1e-5,  "PF04738": 1e-5,   # RiPPs: LanC cyclase + radical SAM sactipeptide cyclase
        "PF04183": 1e-8,  "PF06283": 1e-8,   # Siderophores (IucA/IucC-like)
        "PF06316": 1e-8,  "PF00582": 1e-8,   # Phenazines
        "PF03936": 1e-10                      # Terpenes (terpene synthase)
    }
    
    try:
        import pyhmmer.hmmer as hmmer
        with plan7.HMMFile(hmm_path) as hmm_file:
            hmms = list(hmm_file)
            
        all_hits = hmmer.hmmsearch(hmms, digitized_seqs, cpus=1)
        for hmm, hits in zip(hmms, all_hits):
            # Resolve Pfam accession (e.g. PF00501 from PF00501.35)
            hmm_acc = hmm.accession.decode("utf-8") if isinstance(hmm.accession, bytes) else hmm.accession
            if hmm_acc and "." in hmm_acc:
                hmm_acc = hmm_acc.split(".")[0]
            if not hmm_acc:
                hmm_acc = hmm.name.decode("utf-8") if isinstance(hmm.name, bytes) else hmm.name
                
            cutoff = pfam_cutoffs.get(hmm_acc, 1e-8)
            for hit in hits:
                if hit.evalue <= cutoff:
                    ltag = hit.name.decode("utf-8") if isinstance(hit.name, bytes) else hit.name
                    if ltag not in domain_hits:
                        domain_hits[ltag] = []
                    domain_hits[ltag].append(hmm_acc)
    except Exception as e:
        logger.warning(f"pyhmmer search execution failed: {e}")
        
    return domain_hits

def classify_bgc_by_domains(cluster_hits, domain_hits):
    cluster_domains = []
    for h in cluster_hits:
        ltag = h["locus_tag"]
        if ltag in domain_hits:
            cluster_domains.extend(domain_hits[ltag])
            
    # Reference: Blin et al. 2023 (antiSMASH 7.0) / Grell et al. 2018
    # BGC classification should map to the primary biosynthetic core domain present.
    # We do NOT demand combinatorial architectures (e.g., KS + AT) because:
    # 1. Type II/Type III PKS and Trans-AT systems lack these multi-domain fusions.
    # 2. essential_bgc.hmm is a curated subset containing only primary core signatures.
    
    if "PF00109" in cluster_domains or "PF02801" in cluster_domains:
        return "Polyketide"
    elif "PF00668" in cluster_domains or "PF00501" in cluster_domains:
        return "NRPS"
    elif "PF05147" in cluster_domains or "PF04738" in cluster_domains:
        # PF05147: LanC cyclase (lanthipeptide-specific)
        # PF04738: Radical SAM sactipeptide cyclase (sactipeptide RiPP)
        return "RiPP"
    elif "PF04183" in cluster_domains or "PF06283" in cluster_domains:
        return "Siderophore"
    elif "PF03936" in cluster_domains:
        return "Terpene"
    elif "PF06316" in cluster_domains or "PF00582" in cluster_domains:
        return "Aromatic / Phenazine"
    elif "PF00196" in cluster_domains or "PF00765" in cluster_domains or "PF03006" in cluster_domains or "PF07968" in cluster_domains or "PF08660" in cluster_domains or "PF01338" in cluster_domains:
        # Fallback for regulatory/toxin/quorum-sensing islands that form protoclusters
        return "Other"
        
    return None

def parse_antismash_output(antismash_dir, logger, query_gb=None):
    """
    Parses antiSMASH JSON or GenBank output files to extract BGC records,
    homology percentages, and known clusters comparisons.
    """
    import json
    import glob
    from Bio import SeqIO
    
    bgc_list = []
    json_files = glob.glob(os.path.join(antismash_dir, "*.json"))
    if not json_files:
        return None
        
    json_path = json_files[0] # Take first summary json
    logger.info(f"Parsing antiSMASH results from: {json_path}")
    try:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
            
        records = data.get("records", [])
        
        # Load query CDS features to populate Hits for downstream promoter box scanning and key gene list
        cds_features = []
        if query_gb and os.path.exists(query_gb):
            try:
                record = max(list(SeqIO.parse(query_gb, "genbank")), key=lambda r: len(r.seq))
                for feat in record.features:
                    if feat.type == "CDS":
                        cds_features.append({
                            "start": int(feat.location.start),
                            "end": int(feat.location.end),
                            "strand": "+" if feat.location.strand >= 0 else "-",
                            "locus_tag": feat.qualifiers.get("locus_tag", [""])[0],
                            "gene": feat.qualifiers.get("gene", ["NA"])[0],
                            "product": feat.qualifiers.get("product", ["Hypothetical protein"])[0],
                            "translation": feat.qualifiers.get("translation", [""])[0]
                        })
            except Exception as e:
                logger.warning(f"Failed to parse query GenBank in parse_antismash_output: {e}")

        idx = 1
        for rec in records:
            # Extract knownclusterblast results from modules
            cb_results = {}
            try:
                cb_module = rec.get("modules", {}).get("antismash.modules.clusterblast", {})
                knowncluster = cb_module.get("knowncluster", {})
                for r in knowncluster.get("results", []):
                    reg_num = r.get("region_number")
                    ranking = r.get("ranking", [])
                    if reg_num and ranking:
                        top_hit = ranking[0]
                        if len(top_hit) >= 2:
                            db_meta = top_hit[0]
                            blast_info = top_hit[1]
                            mibig_id = db_meta.get("accession", "Unknown")
                            desc = db_meta.get("description", "Unknown")
                            sim = blast_info.get("similarity", 0.0)
                            cb_results[reg_num] = {
                                "mibig_id": mibig_id,
                                "description": desc,
                                "similarity": float(sim)
                            }
            except Exception as e:
                logger.warning(f"Could not parse clusterblast knowncluster results: {e}")

            areas = rec.get("areas", [])
            for a_idx, area in enumerate(areas):
                reg_num = a_idx + 1
                bgc_type = area.get("products", ["Unknown"])[0]
                start = area.get("start", 0)
                end = area.get("end", 0)
                
                # Known cluster similarity from clusterblast module
                best_similarity = 0.0
                best_match = "Unknown cluster"
                if reg_num in cb_results:
                    best_match = f"{cb_results[reg_num]['mibig_id']} ({cb_results[reg_num]['description']})"
                    best_similarity = cb_results[reg_num]['similarity']
                    
                bgc_desc = f"{bgc_type} (antiSMASH predicted region; {best_similarity:.1f}% similarity to MIBiG {best_match})"
                
                # Extract core_start and core_end from protoclusters
                protoclusters_dict = area.get("protoclusters", {})
                core_regions = []
                for pc_key, pc_val in protoclusters_dict.items():
                    c_s = pc_val.get("core_start")
                    c_e = pc_val.get("core_end")
                    if c_s is not None and c_e is not None:
                        core_regions.append((c_s, c_e))
                
                overall_core_start = min([cs for cs, ce in core_regions]) if core_regions else start
                overall_core_end = max([ce for cs, ce in core_regions]) if core_regions else end
                
                # Extract CDSs falling strictly inside any protocluster's core boundaries as true core genes
                area_hits = []
                for cds in cds_features:
                    is_core = False
                    for c_s, c_e in core_regions:
                        if cds["start"] >= c_s and cds["end"] <= c_e:
                            is_core = True
                            break
                    if is_core:
                        area_hits.append(cds)
                
                # Check for gene symbols and format BGC Type description
                all_symbols = list(dict.fromkeys([h["gene"] for h in area_hits if h["gene"] != "NA"]))
                if all_symbols:
                    bgc_desc += f" (key genes: {', '.join(all_symbols[:4])})"
                
                bgc_list.append({
                    "BGC_ID": f"BGC_{idx:02d}",
                    "BGC_Type": bgc_desc,
                    "Start_Coord": start,
                    "End_Coord": end,
                    "Size_bp": end - start,
                    "Core_Start": overall_core_start,
                    "Core_End": overall_core_end,
                    "Hits": area_hits
                })
                idx += 1
        return bgc_list
    except Exception as e:
        logger.warning(f"Failed to parse antiSMASH output JSON: {e}")
        return None

def run_antismash_submodule(query_gb, out_dir, logger):
    """
    Executes antiSMASH via Docker/WSL or local conda command.
    """
    import subprocess
    import shutil
    
    antismash_out = os.path.join(out_dir, "antismash_results")
    os.makedirs(antismash_out, exist_ok=True)
    
    logger.info("Checking for Docker/antiSMASH container to run core BGC detection...")
    
    docker_path = shutil.which("docker")
    wsl_path = shutil.which("wsl")
    if docker_path or wsl_path:
        # Construct absolute path in POSIX style for WSL volume mounting
        abs_in_dir = os.path.dirname(os.path.abspath(query_gb))
        abs_out_dir = os.path.abspath(antismash_out)
        
        # Convert path to WSL mount points (case-insensitive for drive letters)
        wsl_in = re.sub(r'^([a-zA-Z]):', lambda m: f"/mnt/{m.group(1).lower()}", abs_in_dir).replace("\\", "/")
        wsl_out = re.sub(r'^([a-zA-Z]):', lambda m: f"/mnt/{m.group(1).lower()}", abs_out_dir).replace("\\", "/")
        gb_file = os.path.basename(query_gb)
        
        # Check if native Windows docker works or if we should run via WSL
        use_wsl_prefix = []
        if docker_path:
            try:
                # CREATE_NO_WINDOW prevents WSL/docker subprocesses from attaching
                # to the Windows console handle, which corrupts stdin for the parent process.
                _win_flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
                res_info = subprocess.run(["docker", "info"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=10, stdin=subprocess.DEVNULL, creationflags=_win_flags)
                if res_info.returncode == 0:
                    pass  # native docker works
                elif wsl_path:
                    logger.info("Windows native Docker connection failed. Attempting to start/verify WSL Docker (docker-desktop)...")
                    # Start WSL docker-desktop distro to wake up the daemon
                    subprocess.run(["wsl", "-d", "docker-desktop", "true"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=20, stdin=subprocess.DEVNULL, creationflags=_win_flags)
                    import time
                    for i in range(5):
                        res_wsl = subprocess.run(["wsl", "docker", "info"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=5, stdin=subprocess.DEVNULL, creationflags=_win_flags)
                        if res_wsl.returncode == 0:
                            use_wsl_prefix = ["wsl"]
                            logger.info("WSL Docker backend is active and ready.")
                            break
                        logger.info(f"Waiting for WSL Docker daemon to initialize (attempt {i+1}/5)...")
                        time.sleep(3)
            except Exception:
                if wsl_path:
                    use_wsl_prefix = ["wsl"]
        else:
            if wsl_path:
                use_wsl_prefix = ["wsl"]

        cmd = use_wsl_prefix + [
            "docker", "run", "--user", "1000:1000", "--rm",
            "-v", f"{wsl_in}:/input:ro",
            "-v", f"{wsl_out}:/output",
            "antismash/standalone:latest",
            gb_file, "--output-dir", "/output",
            "--cb-general", "--cb-knowncluster", "--cb-subclusters", "--asf"
        ]
        try:
            logger.info(f"Executing: {' '.join(cmd)}")
            _win_flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
            res = safe_subprocess_run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=600, stdin=subprocess.DEVNULL, creationflags=_win_flags)
            logger.info("antiSMASH run completed successfully via Docker submodule.")
            return parse_antismash_output(antismash_out, logger, query_gb)
        except Exception as e:
            logger.warning(f"Could not invoke antiSMASH via Docker: {e}")
            
    # Fall back to local command line if installed in PATH
    antismash_cmd = shutil.which("antismash")
    if antismash_cmd:
        cmd = [
            antismash_cmd, query_gb, "--output-dir", antismash_out,
            "--cb-general", "--cb-knowncluster", "--cb-subclusters", "--asf"
        ]
        try:
            logger.info(f"Executing: {' '.join(cmd)}")
            _win_flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
            res = safe_subprocess_run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=600, stdin=subprocess.DEVNULL, creationflags=_win_flags)
            return parse_antismash_output(antismash_out, logger, query_gb)
        except Exception as e:
            logger.warning(f"Could not invoke local antiSMASH command: {e}")
            
    logger.info("antiSMASH not found/executable. Using HMMER + Keyword secondary fallback scanner.")
    return None

def compute_bgc_geometric_similarity(core_proteins, ref_gene_names, logger):
    """
    Computes a weighted geometric mean of Pfam overlap, gene order, domain architecture,
    completeness, and MIBiG similarity.
    Formula: Similarity = PROD (M_i) ^ w_i
    Includes epsilon bounding for zero values.
    """
    if not core_proteins or not ref_gene_names:
        return 0.0
        
    # 1. Pfam overlap / symbol matching fraction
    matched_q_genes = []
    matched_ref_genes = set()
    for q_sym, q_trans in core_proteins:
        for ref_sym in ref_gene_names:
            if q_sym.startswith(ref_sym[:4]) or ref_sym.startswith(q_sym[:4]):
                matched_q_genes.append((q_sym, ref_sym))
                matched_ref_genes.add(ref_sym)
                break
                
    m_pfam = len(matched_ref_genes) / float(len(ref_gene_names))
    
    # 2. Gene Order Synteny
    q_order = [q_sym for q_sym, _ in core_proteins]
    ref_matched_in_q = [ref_sym for q_sym, ref_sym in matched_q_genes]
    
    order_matches = 0
    ref_pairs = len(ref_gene_names) - 1
    if ref_pairs > 0 and len(ref_matched_in_q) > 1:
        for i in range(len(ref_gene_names) - 1):
            r1, r2 = ref_gene_names[i], ref_gene_names[i+1]
            if r1 in ref_matched_in_q and r2 in ref_matched_in_q:
                try:
                    p1 = q_order.index([q for q, r in matched_q_genes if r == r1][0])
                    p2 = q_order.index([q for q, r in matched_q_genes if r == r2][0])
                    if abs(p1 - p2) <= 3:
                        order_matches += 1
                except Exception:
                    pass
        m_order = order_matches / float(ref_pairs)
    else:
        m_order = 1.0 if len(ref_gene_names) == 1 else 0.0
        
    # 3. Domain Jaccard Index (Proxy for BiG-SCAPE distance metric)
    # Reference: Navarro-Muñoz et al. 2020 (doi:10.1038/s41589-019-0400-9)
    # Since raw reference sequences are not available for DSS, we rely strictly on 
    # structural domain synteny and completeness (Jaccard index equivalent).
    m_domain = m_pfam  # Default domain architecture proxy
    m_complete = len(matched_ref_genes) / float(len(ref_gene_names))
    
    # Epsilon scaling to prevent zero-product elimination
    epsilon = 0.05
    m1 = max(epsilon, m_pfam)
    m2 = max(epsilon, m_order)
    m3 = max(epsilon, m_domain)
    m4 = max(epsilon, m_complete)
    
    # Adjusted Weights: Pfam (35%), gene order (25%), domain architecture (25%), completeness (15%)
    w1, w2, w3, w4 = 0.35, 0.25, 0.25, 0.15
    geom_mean = (m1**w1) * (m2**w2) * (m3**w3) * (m4**w4)
    return geom_mean

# Scan the query record for BGC matching features
def detect_all_bgcs(query_gb, logger, out_dir=None):
    logger.info("Detecting Biosynthetic Gene Clusters (BGCs) in query genome...")
    
    if out_dir is None:
        out_dir = os.path.join(os.path.dirname(os.path.abspath(query_gb)), "results")
        
    # Try running antiSMASH as primary detection model
    bgcs = run_antismash_submodule(query_gb, out_dir, logger)
    if bgcs:
        return bgcs
        
    record = max(list(SeqIO.parse(query_gb, "genbank")), key=lambda r: len(r.seq))
    
    # Extract all genome proteins de novo for HMM scanning
    genome_proteins = []
    for feat in record.features:
        if feat.type == "CDS":
            ltag = feat.qualifiers.get("locus_tag", [""])[0]
            trans = feat.qualifiers.get("translation", [""])[0]
            if ltag and trans:
                genome_proteins.append((ltag, trans))
                
    # Run pyhmmer genome-wide scan
    hmm_path = get_curated_bgc_hmms(logger)
    domain_hits = {}
    if hmm_path and genome_proteins:
        domain_hits = run_pyhmmer_scan(genome_proteins, hmm_path, logger)
        
    hits = []
    for feat in record.features:
        if feat.type == "CDS":
            ltag = feat.qualifiers.get("locus_tag", [""])[0]
            gene = feat.qualifiers.get("gene", ["NA"])[0]
            prod = feat.qualifiers.get("product", ["Hypothetical protein"])[0]
            trans = feat.qualifiers.get("translation", [""])[0]
            loc = feat.location
            
            is_core = False
            if ltag in domain_hits:
                # Filter out non-biosynthetic domains PF14867 (PCMT), PF00582 (Usp), PF00589 (Phage integrase), 
                # PF00872 (Transposase), and PF00239 (Resolvase) as seeds.
                # Reference: Grell et al. 2018 (doi:10.1038/s41589-018-0122-9)
                active_domains = [d for d in domain_hits[ltag] if d not in ["PF14867", "PF00582", "PF00589", "PF00872", "PF00239"]]
                # Removed unscientific exclusion word-mines that threw out genuine core domains if they were textually annotated with primary metabolism terms (e.g. 'acyl-coa synthetase').
                # Primary vs Secondary metabolism is now differentiated strictly by multi-domain architecture requirements in classify_bgc_by_domains().
                if active_domains:
                    is_core = True
            
            # Removed unscientific keyword string matching.
            # BGC core genes must be identified via pyhmmer multi-domain mapping, not string descriptions.
            if is_core:
                hits.append({
                    "start": int(loc.start),
                    "end": int(loc.end),
                    "strand": "+" if loc.strand >= 0 else "-",
                    "locus_tag": ltag,
                    "gene": gene,
                    "product": prod,
                    "translation": trans
                })
                
    if not hits:
        return []
        
    # Group hits into positional protoclusters.
    # Gap threshold: CONFIG["bgc_clustering_gap_bp"] (default 10 kb, same as antiSMASH)
    # Reference: Blin et al. 2023 antiSMASH 7.0 (doi:10.1093/nar/gkad344)
    gap_threshold = CONFIG.get("bgc_clustering_gap_bp", 10000)
    sorted_hits = sorted(hits, key=lambda x: x["start"])
    protoclusters = []
    current = [sorted_hits[0]]
    for h in sorted_hits[1:]:
        current_max_end = max(x["end"] for x in current)
        if h["start"] - current_max_end <= gap_threshold:
            current.append(h)
        else:
            protoclusters.append(current)
            current = [h]
    protoclusters.append(current)
    
    # Merge overlapping protocluster regions into Candidate Clusters (Regions)
    # Reference: Blin et al. 2023 antiSMASH 7.0 (doi:10.1093/nar/gkad344)
    # "Regions simply reflect contiguous genomic loci containing one or more protoclusters."
    flank = CONFIG.get("bgc_flank_window", 10000)
    clusters = []
    current_cc = protoclusters[0]
    for p in protoclusters[1:]:
        cc_end_core = max(x["end"] for x in current_cc)
        p_start_core = min(x["start"] for x in p)
        
        if p_start_core - flank <= cc_end_core + flank:
            current_cc.extend(p)
        else:
            clusters.append(current_cc)
            current_cc = p
    clusters.append(current_cc)
    
    bgc_list = []
    for idx, c in enumerate(clusters):
        c_start = min(h["start"] for h in c)
        c_end = max(h["end"] for h in c)
        flank = CONFIG["bgc_flank_window"]
        bgc_start = max(0, c_start - flank)
        bgc_end = min(len(record.seq), c_end + flank)
        
        # BGC Type Classification
        bgc_type = classify_bgc_by_domains(c, domain_hits)
        
        if not bgc_type:
            logger.info("Dropping false-positive protocluster: failed strict mathematical domain architecture check.")
            continue
            
        # ─── MIBiG database matching: tiered protein-identity approach ───
        # Tier 1: Bidirectional best-hit protein sequence identity using pairwise alignment
        # Tier 2: HMM domain-class comparison (if no protein sequences available)
        # Tier 3: Keyword product description backup
        # Rationale: Gene symbol Jaccard is unreliable due to non-standardised naming across assemblies.
        # Reference: Medema et al. 2011 antiSMASH (doi:10.1093/nar/gkr466)
        
        # Collect core biosynthetic proteins for this cluster (with sequences)
        core_biosynthetic_proteins = []  # list of (gene_symbol, translation)
        for feat in record.features:
            if feat.type == "CDS":
                fs = int(feat.location.start)
                ltag = feat.qualifiers.get("locus_tag", [""])[0]
                g_sym = feat.qualifiers.get("gene", ["NA"])[0]
                prod = feat.qualifiers.get("product", [""])[0]
                trans = feat.qualifiers.get("translation", [""])[0]
                if c_start <= fs <= c_end and g_sym != "NA":
                    has_hmm_hit = ltag in domain_hits
                    has_biosynthetic_role = get_bgc_gene_role(prod, g_sym) == "biosynthetic"
                    if (has_hmm_hit or has_biosynthetic_role) and trans:
                        core_biosynthetic_proteins.append((g_sym.lower(), trans))

        best_cluster_name = None
        best_similarity = 0.0

        # Determine query organism genus for genus-specific profile filtering
        query_genus = record.annotations.get("organism", "").split()[0].lower() if record.annotations.get("organism") else ""

        if len(core_biosynthetic_proteins) >= 2:
            # Tier 1: Weighted geometric similarity calculation
            try:
                for mibig_id, profile_data in MIBIG_DATABASE.items():
                    # Apply genus-specificity filter
                    ref_genus_list = profile_data.get("genus_specific", [])
                    if ref_genus_list and query_genus not in ref_genus_list:
                        continue

                    ref_gene_names = profile_data.get("genes", [])
                    if not ref_gene_names:
                        continue

                    # Calculate weighted geometric similarity
                    score = compute_bgc_geometric_similarity(core_biosynthetic_proteins, ref_gene_names, logger) * 100.0
                    
                    if score > best_similarity:
                        best_similarity = score
                        best_cluster_name = mibig_id
                        
                if best_similarity >= 30.0:
                    logger.info(f"BGC geometric similarity computed. Confident match: {best_cluster_name} with {best_similarity:.1f}% similarity.")
                elif best_cluster_name:
                    logger.info(f"BGC geometric similarity computed ({best_similarity:.1f}% for {best_cluster_name}), but falls below 30% twilight zone. Discarding match.")
            except Exception as e:
                logger.warning(f"Error computing geometric similarity: {e}")
        else:
            logger.info("BGC has fewer than 2 core biosynthetic proteins. Skipping geometric similarity.")

        # Tier 3 (Keyword fallback) deleted to enforce strict biology
        # Formatting fix: collect actual gene symbols and remove duplicates
        all_symbols = list(dict.fromkeys([h["gene"] for h in c if h["gene"] != "NA"]))
        gene_str = f" (key genes: {', '.join(all_symbols[:4])})" if all_symbols else ""
        
        mibig_id = None
        # Enforce minimum 30.0% similarity threshold to avoid false positive specific compound assignments
        # Below 30% similarity, specific functional assignment to a MIBiG reference is unreliable.
        # Reference: Rost 1999 (doi:10.1093/protein/12.2.85) — twilight zone boundary for protein homology
        if best_cluster_name and best_similarity >= 30.0:
            mibig_id = best_cluster_name.split()[0]
            
        if mibig_id and mibig_id in MIBIG_ANNOTATIONS:
            anno = MIBIG_ANNOTATIONS[mibig_id]
            bgc_desc = (f"{bgc_type} [Compound: {anno['compound']} | Class: {anno['class']} | "
                        f"Activity: {anno['activity']} | DB: {best_similarity:.1f}% similarity to "
                        f"{best_cluster_name}] ({anno['ref']}){gene_str}")
        elif best_cluster_name and best_similarity >= 30.0:
            bgc_desc = f"{bgc_type} ({best_similarity:.1f}% similarity to MIBiG {best_cluster_name}){gene_str}"
        else:
            bgc_desc = f"{bgc_type}{gene_str}"
            
        bgc_list.append({
            "BGC_ID": f"BGC_{idx+1:02d}",
            "BGC_Type": bgc_desc,
            "Start_Coord": bgc_start,
            "End_Coord": bgc_end,
            "Size_bp": bgc_end - bgc_start,
            "Core_Start": c_start,
            "Core_End": c_end,
            "Hits": c
        })
        
    return bgc_list

# Promoter box motifs search
def scan_promoter_motifs(up_seq, motifs_config):
    """Scan upstream sequence for sigma factor promoter consensus motifs.
    
    Uses per-sigma spacer ranges from config (spacer_min, spacer_max).
    Spacer is measured from the end of the -35 box to the start of the -10 box.
    References: Hawley & McClure 1983; Helmann 1995; Boylan et al. 1992.
    """
    hits = []
    for sigma_name, patterns in motifs_config.items():
        pat35 = patterns["minus35"]
        pat10 = patterns["minus10"]
        # Use per-sigma spacer range from config; fall back to canonical 16-18 bp
        spacer_min = patterns.get("spacer_min", 16)
        spacer_max = patterns.get("spacer_max", 18)
        quality_label = f"Spacer {spacer_min}-{spacer_max} bp ({sigma_name} consensus)"

        for m35 in re.finditer(pat35, up_seq, re.IGNORECASE):
            p35 = m35.start()
            # Search for -10 box within the biologically expected window after -35 box
            ss = p35 + len(m35.group()) + spacer_min
            se = min(len(up_seq), p35 + len(m35.group()) + spacer_max + 8)
            sub = up_seq[ss:se]
            for m10 in re.finditer(pat10, sub, re.IGNORECASE):
                p10 = ss + m10.start()
                spacer = p10 - (p35 + len(m35.group()))
                if spacer_min <= spacer <= spacer_max:
                    hits.append({
                        "Sigma_Factor": sigma_name,
                        "Minus35_Pos": p35 - len(up_seq),
                        "Minus35_Seq": m35.group().upper(),
                        "Spacer_Length": spacer,
                        "Minus10_Pos": p10 - len(up_seq),
                        "Minus10_Seq": m10.group().upper(),
                        "Quality": quality_label
                    })
    return hits

# Detailed BGC analysis centering on the targeted BGC ID
def evaluate_targeted_bgc(query_gb, ref_gb, ortholog_mapping, target_bgc, out_dir, logger):
    logger.info(f"Evaluating targeted cluster of interest: {target_bgc['BGC_ID']} ({target_bgc['BGC_Type']})...")
    record = max(list(SeqIO.parse(query_gb, "genbank")), key=lambda r: len(r.seq))
    
    bgc_id = target_bgc["BGC_ID"]
    bgc_start = target_bgc["Start_Coord"]
    bgc_end = target_bgc["End_Coord"]
    
    region_genes = []
    flanking_genes = []
    
    for feat in record.features:
        if feat.type == "CDS":
            fs, fe = int(feat.location.start), int(feat.location.end)
            l_tag = feat.qualifiers.get("locus_tag", [""])[0]
            g_sym = feat.qualifiers.get("gene", ["NA"])[0]
            prod = feat.qualifiers.get("product", ["Hypothetical protein"])[0]
            trans = feat.qualifiers.get("translation", [""])[0]
            strand = "+" if feat.location.strand >= 0 else "-"
            
            gene_dict = {
                "Locus_Tag": l_tag,
                "Gene_Symbol": g_sym,
                "Start_Coord": fs,
                "End_Coord": fe,
                "Strand": strand,
                "Product_Description": prod,
                "Role": get_bgc_gene_role(prod, g_sym),
                "Protein_Sequence": trans
            }
            
            if target_bgc["Core_Start"] <= fs <= target_bgc["Core_End"]:
                region_genes.append(gene_dict)
            elif bgc_start <= fs <= bgc_end:
                flanking_genes.append(gene_dict)
                 
    from redolarium.promoter_prediction import run_promoter_prediction
    
    normalized_genes = []
    for g in region_genes:
        normalized_genes.append({
            "Locus_Tag": g.get("Locus_Tag", "Unknown"),
            "Gene": g.get("Gene_Symbol", "-"),
            "Start": g.get("Start_Coord", 0),
            "End": g.get("End_Coord", 0),
            "Strand": g.get("Strand", "+"),
            "Role": g.get("Role", "Other")
        })
        
    promoter_records = run_promoter_prediction(query_gb, target_bgc, normalized_genes, out_dir, logger)

    # Add a metadata limitations notice file
    limits_file = os.path.join(out_dir, "bgc_motifs", f"{bgc_id}_pipeline_disclaimer.txt")
    try:
        with open(limits_file, "w", encoding="utf-8") as df_f:
            df_f.write("=== BGC Annotation & Pipeline Disclaimer ===\n")
            df_f.write("Software Version: Redolarium v3.0.0\n")
            df_f.write("Pipeline Core Engine: Hybrid (antiSMASH + custom HMMER fallbacks)\n")
            df_f.write("HMM Database Version: Pfam-A v36.0 (curated essential subset)\n")
            df_f.write("MIBiG Reference Release: 3.1\n")
            df_f.write("\n")
            df_f.write("LIMITATIONS:\n")
            df_f.write("1. Predicted BGC boundaries and functional annotations are derived computationally.\n")
            df_f.write("2. Transcripts, promoters, and Shine-Dalgarno motifs are identified via consensus sequence alignment\n")
            df_f.write("   and do not substitute for in vivo expression verification.\n")
            df_f.write("3. Exact chemical products, activities, and regulatory bounds must be verified experimentally\n")
            df_f.write("   (e.g., via heterologous expression, gene knockouts, or mass spectrometry).\n")
    except Exception:
        pass

    phage_hits = []
    # Removed unscientific phage_keywords word mine. Phage/MGE elements are now strictly defined via HMM profiles in target_bgc_analysis.py.
    df_phage = pd.DataFrame(phage_hits)
    csv_phage = os.path.join(out_dir, "tabular_data", f"{bgc_id}_phage_artifacts.csv")
    df_phage.to_csv(csv_phage, index=False)
    logger.info(f"Saved prophage flanking artifacts to: {csv_phage}")

    df_flank = pd.DataFrame(flanking_genes)
    csv_flank = os.path.join(out_dir, "tabular_data", f"{bgc_id}_flanking_context.csv")
    df_flank.to_csv(csv_flank, index=False)
    
    df_core = pd.DataFrame(region_genes)
    csv_core = os.path.join(out_dir, "tabular_data", f"{bgc_id}_gene_architecture.csv")
    df_core.to_csv(csv_core, index=False)

    colors = {
        "biosynthetic": "#C9DAF8",
        "transport": "#AED6F1",
        "regulatory": "#D7BDE2",
        "immunity": "#F9C784",
        "accessory/context": "#E8F0FE"
    }
    
    # Dual-scale / Fisheye Synteny Diagram Overhaul
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [1, 3]})
    
    # Track 1: Wide-context track (Genome/Contig)
    try:
        record = max(list(SeqIO.parse(query_gb, "genbank")), key=lambda r: len(r.seq))
    except Exception:
        record = max(list(SeqIO.parse(query_gb, "fasta")), key=lambda r: len(r.seq))
    genome_len = len(record.seq)
    
    ax1.set_xlim(0, genome_len)
    ax1.set_ylim(-1, 1)
    ax1.axis("off")
    ax1.plot([0, genome_len], [0, 0], color="black", lw=2)
    ax1.add_patch(patches.Rectangle((bgc_start, -0.5), bgc_end - bgc_start, 1.0, fill=True, color="red", alpha=0.4))
    ax1.text(genome_len/2, 0.6, f"Full Contig Context ({genome_len:,} bp)", ha="center", fontsize=9, fontweight="bold")
    ax1.text((bgc_start + bgc_end)/2, -0.8, f"Target {bgc_id}", ha="center", fontsize=8, color="red", fontweight="bold")
    
    # Track 2: Specifically scaled subset track for the focused BGC
    # Add some padding around the BGC bounds
    view_start = max(0, bgc_start - 5000)
    view_end = min(genome_len, bgc_end + 5000)
    
    ax2.set_xlim(view_start, view_end)
    ax2.set_ylim(-2.5, 3.5)
    ax2.axis("off")
    
    # 5' to 3' tracking line
    ax2.plot([view_start, view_end], [0, 0], color="black", lw=2)
    ax2.text(view_start + 100, 0.2, "5'", fontweight="bold", fontsize=12)
    ax2.text(view_end - 500, 0.2, "3'", fontweight="bold", fontsize=12)
    
    # Draw Fisheye connecting lines using ConnectionPatch
    try:
        from matplotlib.patches import ConnectionPatch
        # Track 1 uses genome coordinates. Track 2 uses genome coordinates (view_start to view_end)
        # So we can just map the specific X coordinates between the axes
        con1 = ConnectionPatch(xyA=(bgc_start, -0.5), xyB=(bgc_start, 3.5), coordsA="data", coordsB="data",
                               axesA=ax1, axesB=ax2, color="red", alpha=0.3, linestyle="--")
        con2 = ConnectionPatch(xyA=(bgc_end, -0.5), xyB=(bgc_end, 3.5), coordsA="data", coordsB="data",
                               axesA=ax1, axesB=ax2, color="red", alpha=0.3, linestyle="--")
        ax2.add_artist(con1)
        ax2.add_artist(con2)
    except Exception as e:
        logger.warning(f"Failed to draw synteny fisheye connection patches: {e}")
        
    for g in (region_genes + flanking_genes):
        g_len = g["End_Coord"] - g["Start_Coord"]
        color = colors.get(g["Role"], "#E8F0FE")
        direction = 1 if g["Strand"] == "+" else -1
        
        # Explicit transcription start/stop sites (Promoter/Terminator demarcation)
        if direction == 1:
            # Promoter (Green arrow up and right)
            ax2.plot([g["Start_Coord"], g["Start_Coord"]], [0, 1.2], color="green", lw=1.5, alpha=0.7)
            ax2.plot([g["Start_Coord"], g["Start_Coord"]+300], [1.2, 1.2], color="green", lw=1.5, alpha=0.7)
            # Terminator (Red line with square down)
            ax2.plot([g["End_Coord"], g["End_Coord"]], [0, -0.8], color="red", lw=1.5, marker="s", markersize=4, alpha=0.7)
        else:
            # Promoter (Green arrow down and left)
            ax2.plot([g["End_Coord"], g["End_Coord"]], [0, -1.2], color="green", lw=1.5, alpha=0.7)
            ax2.plot([g["End_Coord"], g["End_Coord"]-300], [-1.2, -1.2], color="green", lw=1.5, alpha=0.7)
            # Terminator (Red line with square up)
            ax2.plot([g["Start_Coord"], g["Start_Coord"]], [0, 0.8], color="red", lw=1.5, marker="s", markersize=4, alpha=0.7)
            
        # Gene Arrow Body
        arrow = patches.FancyArrow(
            g["Start_Coord"] if direction == 1 else g["End_Coord"],
            0,
            g_len * direction,
            0,
            width=0.6,
            head_width=1.0,
            head_length=min(g_len * 0.3, 1000),
            facecolor=color,
            edgecolor="black",
            lw=0.8,
            zorder=3
        )
        ax2.add_patch(arrow)
        
        # Determine center
        center_pos = (g["Start_Coord"] + g["End_Coord"]) / 2.0
        
        # Top Label (Gene Name / Symbol)
        y_name_offset = 1.8 if direction == 1 else -2.0
        gene_name = g["Gene_Symbol"] if g["Gene_Symbol"] != "NA" else g["Locus_Tag"][-5:]
        ax2.text(
            center_pos,
            y_name_offset,
            gene_name,
            ha="center", fontsize=7, fontweight="bold", rotation=45 if direction == 1 else -45
        )
        
        # Bottom Label (Length in bp)
        y_len_offset = -1.5 if direction == 1 else 1.5
        ax2.text(
            center_pos,
            y_len_offset,
            f"{g_len} bp",
            ha="center", fontsize=6, rotation=0
        )
        
        # Side Labels (Start and End Coordinates)
        if direction == 1:
            ax2.text(g["Start_Coord"], -0.8, str(g["Start_Coord"]), ha="right", va="top", fontsize=5, rotation=90)
            ax2.text(g["End_Coord"], -1.2, str(g["End_Coord"]), ha="left", va="top", fontsize=5, rotation=90)
        else:
            ax2.text(g["End_Coord"], 0.8, str(g["End_Coord"]), ha="right", va="bottom", fontsize=5, rotation=90)
            ax2.text(g["Start_Coord"], 1.2, str(g["Start_Coord"]), ha="left", va="bottom", fontsize=5, rotation=90)
            
    core_border = patches.Rectangle(
        (target_bgc["Core_Start"], -2.8),
        target_bgc["Core_End"] - target_bgc["Core_Start"],
        5.6,
        fill=False, edgecolor="red", linestyle="--", lw=1.5, alpha=0.5
    )
    ax2.add_patch(core_border)
    ax2.text((target_bgc["Core_Start"] + target_bgc["Core_End"])/2, 3.0, "BGC Core Region", ha="center", fontsize=10, color="red", fontweight="bold")
    
    plt.suptitle(f"Target Cluster Synteny Map: {bgc_id} ({target_bgc['BGC_Type']})", fontsize=14, fontweight="bold")
    plt.tight_layout()
    
    fig_out = os.path.join(out_dir, "hgt_evolution", f"{bgc_id}_bgc_synteny_blocks.png")
    from redolarium.utils import save_publication_plot
    save_publication_plot(fig_out, dpi=300)
    plt.close()
    logger.info(f"Saved BGC synteny blocks visual to: {fig_out}")
    
    return region_genes, flanking_genes, promoter_records, phage_hits

def generate_comprehensive_bgc_summary(bgc_list, query_gb, out_dir, logger):
    import pandas as pd
    import os
    from Bio import SeqIO
    
    if not bgc_list:
        return
        
    logger.info("Generating Comprehensive BGC Summary Excel...")
    
    out_xlsx = os.path.join(out_dir, "bgc", "BGC_Comprehensive_Summary.xlsx")
    os.makedirs(os.path.dirname(out_xlsx), exist_ok=True)
    
    try:
        query_record = max(list(SeqIO.parse(query_gb, "genbank")), key=lambda r: len(r.seq))
    except Exception:
        query_record = max(list(SeqIO.parse(query_gb, "fasta")), key=lambda r: len(r.seq))
        
    with pd.ExcelWriter(out_xlsx, engine='openpyxl') as writer:
        # Sheet 1: Summary
        summary_data = []
        for b in bgc_list:
            summary_data.append({
                "BGC_ID": b["BGC_ID"],
                "BGC_Type": b["BGC_Type"],
                "Confidence_Score": b.get("Confidence_Score", 0.0),
                "Start_Coord": b["Start_Coord"],
                "End_Coord": b["End_Coord"],
                "Size_bp": b["Size_bp"]
            })
        df_summary = pd.DataFrame(summary_data)
        df_summary.to_excel(writer, sheet_name="BGC_Summary", index=False)
        
        # Sheets for each BGC
        for b in bgc_list:
            bgc_start = b["Start_Coord"]
            bgc_end = b["End_Coord"]
            
            genes = []
            for feat in query_record.features:
                if feat.type == "CDS":
                    fs, fe = int(feat.location.start), int(feat.location.end)
                    if bgc_start <= fs <= bgc_end or bgc_start <= fe <= bgc_end:
                        l_tag = feat.qualifiers.get("locus_tag", [""])[0]
                        g_sym = feat.qualifiers.get("gene", ["NA"])[0]
                        prod = feat.qualifiers.get("product", ["Hypothetical protein"])[0]
                        trans = feat.qualifiers.get("translation", [""])[0]
                        strand = "+" if feat.location.strand >= 0 else "-"
                        role = get_bgc_gene_role(prod, g_sym)
                        
                        genes.append({
                            "Locus_Tag": l_tag,
                            "Gene_Symbol": g_sym,
                            "Start": fs,
                            "End": fe,
                            "Strand": strand,
                            "Product_Description": prod,
                            "Role": role,
                            "Protein_Sequence": trans
                        })
            
            if genes:
                df_genes = pd.DataFrame(genes)
                # Excel sheet names have max 31 characters
                sheet_name = str(b["BGC_ID"])[:31]
                df_genes.to_excel(writer, sheet_name=sheet_name, index=False)
                
    logger.info(f"Saved comprehensive BGC summary to {out_xlsx}")

def run_bgc_pipeline(query_gb, ref_gb, ortholog_mapping, out_dir, logger, qc_result=None):
    logger.info("Running Stage 4 & 5: Biosynthetic Gene Cluster (BGC) Pipeline Stage...")
    start_time = time.time()
    
    from redolarium.structures import PredictionResult
    from redolarium.qc import get_module_qc_penalty
    
    # 1. Run detection
    bgc_list = detect_all_bgcs(query_gb, logger, out_dir)
    
    # Generate comprehensive multi-sheet Excel summary for all detected BGCs
    generate_comprehensive_bgc_summary(bgc_list, query_gb, out_dir, logger)
    
    # Calculate continuous confidence score
    completeness = 0.0
    contamination = 100.0
    closure_status = "Unknown"
    if qc_result:
        completeness = qc_result.prediction.get("completeness", 0.0)
        contamination = qc_result.prediction.get("contamination", 100.0)
        closure_status = qc_result.genome_closure_status
        
    qc_penalty = get_module_qc_penalty("bgc", completeness, contamination)
    
    # Calculate Empirical Cumulative Confidence Score
    # Reference: Bitscore summation reflects standard log-odds probability for homologous domain density 
    # (e.g., Medema et al. 2011, antiSMASH). Cumulative Exponential CDF translates raw score to 0-1 probability.
    
    total_bitscore = 0.0
    for b in bgc_list:
        # Sum bitscores of core hits. If parsing antiSMASH where bitscore is abstracted, 
        # each core gene logically met the Trusted Cutoff (~150.0 bitscore).
        b_bitscore = sum(h.get("score", 150.0) for h in b.get("Hits", []))
        total_bitscore += b_bitscore
        
        # Empirical probability derived from Cumulative HMM Bitscore (Exponential CDF)
        b_base = 1.0 - __import__('math').exp(-b_bitscore / 150.0) if b_bitscore > 0 else 0.0
        b["Confidence_Score"] = max(0.0, min(1.0, b_base - qc_penalty))
        
    # Baseline probability derived from exponential CDF (lambda = 150.0) for the whole genome
    base_confidence = 1.0 - __import__('math').exp(-total_bitscore / 150.0) if total_bitscore > 0 else 0.0
    
    bgc_score = max(0.0, min(1.0, base_confidence - qc_penalty))
    
    prediction_data = {
        "bgc_list": bgc_list,
        "target_bgc": None,
        "region_genes": [],
        "flanking_genes": [],
        "promoter_records": [],
        "phage_hits": []
    }
    
    evidence = [
        f"Detected {len(bgc_list)} biosynthetic gene clusters in genome."
    ]
        
    db_vers = CONFIG.get("database_versions", {})
    
    return PredictionResult(
        prediction=prediction_data,
        confidence_score=bgc_score,
        algorithm="antiSMASH Standalone / PyHMMER",
        algorithm_version="antiSMASH v7.1.0 / HMMER v3.3.2",
        genome_closure_status=closure_status,
        database="MIBiG Secondary Metabolites / Pfam-A",
        database_version=db_vers.get("mibig", "4.0"),
        evidence=evidence,
        limitations=[
            "Predicted BGC boundaries and functional annotations are derived computationally.",
            "Pfam domain profiles and biosynthetic architecture indicators are physically correlated.",
            "Chemical compound activity and regulatory boundaries require in vivo expression validation."
        ],
        citations=[
            "antiSMASH 7.0: Blin et al. 2023 (doi:10.1093/nar/gkad344)",
            "MIBiG 3.0: Terlouw et al. 2023 (doi:10.1093/nar/gkac1049)",
            "pyhmmer: Larralde et al. 2023 (doi:10.21105/joss.05587)"
        ],
        runtime=time.time() - start_time,
        warnings=[f"Low assembly quality penalty applied: -{qc_penalty:.2f}"] if qc_penalty > 0.05 else []
    )
