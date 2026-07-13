# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Redolarium Contributors, IIT Kanpur
# See LICENSE and THIRD_PARTY_LICENSES.md for full licence details.

import os
import re
import time
import json
import csv
import io
import datetime
import urllib.request
import urllib.parse
import shutil
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from redolarium.config import CONFIG, EXCEL_PALETTE
from redolarium.utils import api_retry, save_publication_plot

# 1. Externalized Gene Symbol to EC mapping loading
def load_gene_symbol_to_ec():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    json_path = os.path.join(base_dir, "config", "gene_symbol_to_ec.json")
    if os.path.exists(json_path):
        try:
            with open(json_path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            print(f"Warning: Failed to load externalized gene symbol map: {e}")
    # Fallback default map
    return {
        'alas': '6.1.1.7', 'args': '6.1.1.19', 'asns': '6.1.1.22', 'asps': '6.1.1.12',
        'cyss': '6.1.1.16', 'glns': '6.1.1.18', 'gltx': '6.1.1.17', 'glys': '6.1.1.14',
        'glyq': '6.1.1.14', 'hiss': '6.1.1.21', 'iles': '6.1.1.5', 'leus': '6.1.1.4',
        'lyss': '6.1.1.6', 'mets': '6.1.1.10', 'phes': '6.1.1.20', 'phet': '6.1.1.20',
        'pros': '6.1.1.15', 'sers': '6.1.1.11', 'thrs': '6.1.1.3', 'trps': '6.1.1.2',
        'tyrs': '6.1.1.1', 'vals': '6.1.1.9',
        'pgi': '5.3.1.9', 'pfka': '2.7.1.11', 'pfkb': '2.7.1.11', 'fbaa': '4.1.2.13',
        'fbab': '4.1.2.13', 'tpi': '5.3.1.1', 'tpia': '5.3.1.1', 'gapa': '1.2.1.12',
        'gapb': '1.2.1.12', 'pgk': '2.7.2.3', 'gpma': '5.4.2.11', 'gpmb': '5.4.2.12',
        'eno': '4.2.1.11', 'pyka': '2.7.1.40', 'pykf': '2.7.1.40', 'pdha': '1.2.4.1',
        'pdhb': '1.2.4.1', 'glta': '2.3.3.1', 'citz': '2.3.3.1', 'acna': '4.2.1.3',
        'acnb': '4.2.1.3', 'icd': '1.1.1.42', 'suca': '1.2.4.2', 'succ': '6.2.1.5',
        'sucd': '6.2.1.5', 'sdha': '1.3.5.1', 'sdhb': '1.3.5.1', 'sdhc': '1.3.5.1',
        'sdhd': '1.3.5.1', 'fuma': '4.2.1.2', 'fumb': '4.2.1.2', 'fumc': '4.2.1.2',
        'mdh': '1.1.1.37'
    }

# Fetch wrapper respecting global config (proxy, timeout)
@api_retry(retries=3, backoff_factor=2.0)
def fetch_kegg_url(url, logger=None):
    meta_cfg = CONFIG.get("metabolism", {})
    timeout = meta_cfg.get("kegg_timeout", 30)
    proxy = meta_cfg.get("kegg_proxy", None)

    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    if proxy:
        proxy_handler = urllib.request.ProxyHandler({'http': proxy, 'https': proxy})
        opener = urllib.request.build_opener(proxy_handler)
    else:
        opener = urllib.request.build_opener()
    
    if logger:
        logger.info(f"KEGG REST: Fetching {url}")
    with opener.open(req, timeout=timeout) as response:
        return response.read().decode('utf-8')

def get_kegg_info_version(logger=None):
    try:
        info_text = fetch_kegg_url("https://rest.kegg.jp/info/kegg", logger=logger)
        for line in info_text.splitlines():
            if "Release" in line:
                return line.strip()
    except Exception:
        pass
    return "Unknown KEGG Release"

# Robust Caching Implementation
def load_or_fetch_kegg_data(filename, url, force_refresh, logger):
    meta_cfg = CONFIG.get("metabolism", {})
    cache_dir = meta_cfg.get("kegg_cache_dir", "resources/kegg_cache")
    if not os.path.isabs(cache_dir):
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        cache_dir = os.path.join(project_root, cache_dir)
    expiry_days = meta_cfg.get("kegg_cache_expiry_days", 7)
    
    os.makedirs(cache_dir, exist_ok=True)
    file_path = os.path.join(cache_dir, filename)
    meta_path = os.path.join(cache_dir, "cache_metadata.json")
    
    # Load metadata
    metadata = {}
    if os.path.exists(meta_path):
        try:
            with open(meta_path, "r", encoding="utf-8") as f:
                metadata = json.load(f)
        except Exception:
            pass
            
    is_expired = True
    if filename in metadata:
        fetch_time_str = metadata[filename].get("fetch_date")
        if fetch_time_str:
            try:
                fetch_time = datetime.datetime.fromisoformat(fetch_time_str)
                age = datetime.datetime.now() - fetch_time
                if age.days < expiry_days:
                    is_expired = False
            except Exception:
                pass
                
    if os.path.exists(file_path) and not force_refresh and not is_expired:
        logger.info(f"Reading cached KEGG data for {filename}")
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                return f.read()
        except Exception as e:
            logger.warning(f"Error reading cache for {filename}: {e}. Will refetch.")
            
    # Fetch and Cache
    try:
        data = fetch_kegg_url(url, logger=logger)
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(data)
            
        metadata[filename] = {
            "fetch_date": datetime.datetime.now().isoformat(),
            "kegg_release": get_kegg_info_version(logger)
        }
        with open(meta_path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=4)
        return data
    except Exception as e:
        if os.path.exists(file_path):
            logger.warning(f"Failed to fetch {url} due to {e}. Using expired local cache as fallback.")
            with open(file_path, "r", encoding="utf-8") as f:
                return f.read()
        else:
            logger.error(f"Critical error fetching {url} and no cached local version found.")
            raise e

# BRITE Hierarchy Parser for Pathway Categorization
def load_brite_pathway_categories(force_refresh, logger):
    # Parse br:ko00001 (or br:map00001) to build a mapping of pathway to high-level category
    brite_data = load_or_fetch_kegg_data("brite_ko00001.keg", "https://rest.kegg.jp/get/br:ko00001", force_refresh, logger)
    pathway_to_cat = {}
    
    current_cat_a = "Other Metabolic Pathways"
    current_cat_b = "Other Metabolic Pathways"
    
    for line in brite_data.splitlines():
        if line.startswith("A"):
            # Top-level e.g. "Metabolism"
            current_cat_a = line[1:].strip()
            # Strip HTML tags if present
            current_cat_a = re.sub('<[^<]+?>', '', current_cat_a)
        elif line.startswith("B"):
            # Second-level e.g. "Carbohydrate metabolism"
            current_cat_b = line[1:].strip()
            current_cat_b = re.sub('<[^<]+?>', '', current_cat_b)
        elif line.startswith("C"):
            # Pathway line e.g. "  00010 Glycolysis / Gluconeogenesis [PATH:ko00010]"
            match = re.search(r"(\d{5})", line)
            if match:
                path_id = match.group(1)
                # Map pathway code to Category B under Metabolism
                if "Metabolism" in current_cat_a or "metabolism" in current_cat_b.lower():
                    # Strip leading digits and whitespace (e.g. "  09101 Carbohydrate metabolism" -> "Carbohydrate metabolism")
                    clean_cat_b = re.sub(r'^\s*\d+\s*', '', current_cat_b).strip()
                    pathway_to_cat[path_id] = clean_cat_b
                    pathway_to_cat[f"map{path_id}"] = clean_cat_b
                    pathway_to_cat[f"ko{path_id}"] = clean_cat_b
                else:
                    pathway_to_cat[path_id] = "Other Metabolic Pathways"
                    pathway_to_cat[f"map{path_id}"] = "Other Metabolic Pathways"
                    pathway_to_cat[f"ko{path_id}"] = "Other Metabolic Pathways"
                    
    return pathway_to_cat

def get_category_color(category):
    cat_lower = category.lower()
    if "carbohydrate" in cat_lower:
        return "#8BC34A" # Light Green
    elif "energy" in cat_lower:
        return "#FFC107" # Amber
    elif "lipid" in cat_lower:
        return "#03A9F4" # Light Blue
    elif "nucleotide" in cat_lower:
        return "#FF5722" # Deep Orange
    elif "amino acid" in cat_lower:
        return "#E91E63" # Pink
    elif "cofactors" in cat_lower or "vitamin" in cat_lower:
        return "#FF9800" # Orange
    elif "glycan" in cat_lower:
        return "#CDDC39" # Lime
    elif "secondary" in cat_lower or "biosynthesis" in cat_lower:
        return "#9C27B0" # Purple
    elif "xenobiotics" in cat_lower or "biodegradation" in cat_lower or "transport" in cat_lower:
        return "#00BCD4" # Cyan
    else:
        return "#9E9E9E" # Grey

def export_to_sbml(genes_with_ec, reconstructed_data, filepath, logger):
    """
    Creates a proper SBML Level 3 model using cobra and writes it to filepath.
    Includes gene-protein-reaction (GPR) associations.
    """
    try:
        import cobra
        model = cobra.Model("reconstructed_metabolism")
        
        # Add genes
        added_genes = {}
        for gene in genes_with_ec:
            ltag = gene["Locus_Tag"]
            sym = gene["Gene_Symbol"]
            if ltag not in added_genes:
                c_gene = cobra.Gene(id=ltag, name=sym if sym != "NA" else ltag)
                model.genes.append(c_gene)
                added_genes[ltag] = c_gene
                
        # Add reactions based on EC numbers
        added_reactions = set()
        for row in reconstructed_data:
            ecs = row.get("EC_Numbers", [])
            ltag = row.get("Locus_Tag")
            
            for ec in ecs:
                rxn_id = f"R_EC_{ec.replace('.', '_')}"
                if rxn_id not in added_reactions:
                    # Create generic reaction
                    rxn = cobra.Reaction(id=rxn_id, name=f"Reaction catalyzed by EC {ec}")
                    model.add_reactions([rxn])
                    added_reactions.add(rxn_id)
                    
                # Setup GPR association
                rxn = model.reactions.get_by_id(rxn_id)
                if ltag:
                    if rxn.gene_reaction_rule:
                        rxn.gene_reaction_rule = f"{rxn.gene_reaction_rule} or {ltag}"
                    else:
                        rxn.gene_reaction_rule = ltag
                        
        cobra.io.write_sbml_model(model, filepath)
        logger.info(f"Saved validated SBML Level 3 model to: {filepath}")
    except Exception as e:
        logger.error(f"Failed to generate validated SBML model: {e}")

def run_metabolism_pipeline(query_gb, ortholog_mapping, out_dir, logger, force_refresh=False):
    logger.info("Stage 3: Reconstructing Metabolic Pathways via KEGG REST API...")
    
    meta_cfg = CONFIG.get("metabolism", {})
    mode = meta_cfg.get("mode", "fast")
    
    if mode == "deep":
        logger.info("Deep metabolism mode active. Checking local gapseq and MetaPathPredict installations...")
        gapseq_avail = shutil.which("gapseq") is not None
        metapath_avail = shutil.which("MetaPathPredict") is not None
        if not gapseq_avail:
            logger.warning("gapseq command not found in PATH. Falling back to KEGG fast mode for reconstruction.")
        if not metapath_avail:
            logger.warning("MetaPathPredict command not found in PATH.")
            
    query_record = max(list(SeqIO.parse(query_gb, "genbank")), key=lambda r: len(r.seq))
    query_org = query_record.annotations.get("organism", "Query Species")
    
    # Load gene symbol to EC mapping
    gene_symbol_to_ec = load_gene_symbol_to_ec()
    
    # 1. Collect genes with EC or KO numbers from query genome
    genes_with_ec = []
    orth_dict = {g["Locus_Tag"]: g for g in ortholog_mapping}
    
    stats_counts = {"qualifier_ec": 0, "qualifier_ko": 0, "symbol_map": 0, "regex_product": 0}
    
    for feat in query_record.features:
        if feat.type == "CDS":
            ltag = feat.qualifiers.get("locus_tag", [""])[0]
            if not ltag:
                continue
            ecs = list(feat.qualifiers.get("EC_number", []))
            kos = []
            
            # Extract KO numbers (from /ko, /kegg, /db_xref)
            for qual in ["ko", "kegg", "kegg_orthology"]:
                for val in feat.qualifiers.get(qual, []):
                    m_ko = re.search(r"K\d{5}", val, re.IGNORECASE)
                    if m_ko:
                        kos.append(m_ko.group(0).upper())
                        
            prod = feat.qualifiers.get("product", ["Hypothetical protein"])[0]
            gene_symbol = feat.qualifiers.get("gene", ["NA"])[0]
            
            # Check db_xref for EC numbers or KO numbers
            db_xrefs = feat.qualifiers.get("db_xref", [])
            for ref in db_xrefs:
                if ref.lower().startswith("ec:"):
                    ec_val = ref.split(":")[1].strip()
                    if ec_val not in ecs:
                        ecs.append(ec_val)
                elif ref.upper().startswith("KO:K"):
                    ko_val = ref.split(":")[1].strip().upper()
                    if ko_val not in kos:
                        kos.append(ko_val)
                        
            mapped_method = None
            if ecs:
                mapped_method = "qualifier_ec"
                stats_counts["qualifier_ec"] += 1
                
            # If no EC, try KO mapping
            if not ecs and kos:
                mapped_method = "qualifier_ko"
                stats_counts["qualifier_ko"] += 1
                
            # Map standard gene symbol if no EC/KO was parsed
            if not ecs and not kos and gene_symbol != "NA":
                sym_lower = gene_symbol.lower()
                if sym_lower in gene_symbol_to_ec:
                    ecs.append(gene_symbol_to_ec[sym_lower])
                    mapped_method = "symbol_map"
                    stats_counts["symbol_map"] += 1
                    
            # Fallback to regex search on product description
            if not ecs and not kos:
                m = re.search(r"ec[:\s]+(\d+\.\d+\.\d+\.\d+)", prod, re.IGNORECASE)
                if m:
                    ecs.append(m.group(1))
                    mapped_method = "regex_product"
                    stats_counts["regex_product"] += 1
            
            if ecs or kos:
                genes_with_ec.append({
                    "Locus_Tag": ltag,
                    "Gene_Symbol": gene_symbol,
                    "Product_Description": prod,
                    "EC_Numbers": ecs,
                    "KO_Numbers": kos,
                    "Mapped_Method": mapped_method or "N/A",
                    "Category": orth_dict.get(ltag, {}).get("Category", "Core")
                })
                
    logger.info(f"Annotation mapping summary: EC qualifier: {stats_counts['qualifier_ec']}, KO qualifier: {stats_counts['qualifier_ko']}, Gene Symbol: {stats_counts['symbol_map']}, Product Regex: {stats_counts['regex_product']}")
                    
    if not genes_with_ec:
        logger.warning("[WARNING] No EC or KO numbers detected in query sequence features. Skipping metabolic pathway reconstruction.")
        os.makedirs(os.path.join(out_dir, "metabolic_pathways"), exist_ok=True)
        csv_out = os.path.join(out_dir, "metabolic_pathways", "kegg_enzymes.csv")
        pd.DataFrame(columns=["Locus_Tag", "Gene_Symbol", "Product_Description", "EC_Numbers", "KO_Numbers", "Category"]).to_csv(csv_out, index=False)
        return []

    # 2. Fetch/Load real pathway list, ec-links, ko-links, and categories
    pathway_names = {}
    ec_to_pathway = {}
    ko_to_pathway = {}
    
    logger.info("Fetching / Loading KEGG pathway list...")
    pw_data = load_or_fetch_kegg_data("pathways.tsv", "https://rest.kegg.jp/list/pathway", force_refresh, logger)
    if pw_data:
        for row in csv.reader(io.StringIO(pw_data.strip()), delimiter='\t'):
            if len(row) >= 2:
                pid = row[0].replace('path:', '').strip()
                pathway_names[pid] = row[1].strip()
                
    logger.info("Fetching / Loading KEGG EC-to-pathway links...")
    ec_data = load_or_fetch_kegg_data("ec_links.tsv", "https://rest.kegg.jp/link/pathway/ec", force_refresh, logger)
    if ec_data:
        for row in csv.reader(io.StringIO(ec_data.strip()), delimiter='\t'):
            if len(row) >= 2:
                val1, val2 = row[0].strip(), row[1].strip()
                if val1.startswith('path:'):
                    path = val1.replace('path:', '')
                    ec = val2.replace('ec:', '')
                else:
                    ec = val1.replace('ec:', '')
                    path = val2.replace('path:', '')
                if ec not in ec_to_pathway:
                    ec_to_pathway[ec] = []
                ec_to_pathway[ec].append(path)
                
    logger.info("Fetching / Loading KEGG KO-to-pathway links...")
    ko_data = load_or_fetch_kegg_data("ko_links.tsv", "https://rest.kegg.jp/link/pathway/ko", force_refresh, logger)
    if ko_data:
        for row in csv.reader(io.StringIO(ko_data.strip()), delimiter='\t'):
            if len(row) >= 2:
                val1, val2 = row[0].strip(), row[1].strip()
                if val1.startswith('path:'):
                    path = val1.replace('path:', '')
                    ko = val2.replace('ko:', '')
                else:
                    ko = val1.replace('ko:', '')
                    path = val2.replace('path:', '')
                if ko not in ko_to_pathway:
                    ko_to_pathway[ko] = []
                ko_to_pathway[ko].append(path)
                
    # Load BRITE pathway categories
    logger.info("Parsing KEGG BRITE pathway hierarchy...")
    pathway_to_category = load_brite_pathway_categories(force_refresh, logger)
    
    ec_to_ko = {}
    logger.info("Fetching / Loading KEGG EC-to-KO links...")
    ec_to_ko_data = load_or_fetch_kegg_data("ec_to_ko.tsv", "https://rest.kegg.jp/link/ko/ec", force_refresh, logger)
    if ec_to_ko_data:
        for row in csv.reader(io.StringIO(ec_to_ko_data.strip()), delimiter='\t'):
            if len(row) >= 2:
                val1, val2 = row[0].strip(), row[1].strip()
                if val1.startswith('ko:'):
                    ko = val1.replace('ko:', '')
                    ec = val2.replace('ec:', '')
                else:
                    ec = val1.replace('ec:', '')
                    ko = val2.replace('ko:', '')
                if ec not in ec_to_ko:
                    ec_to_ko[ec] = []
                if ko not in ec_to_ko[ec]:
                    ec_to_ko[ec].append(ko)
                
    if not pathway_names or (not ec_to_pathway and not ko_to_pathway):
        raise Exception("CRITICAL: Failed to load pathway maps or enzyme links from KEGG. Halting execution.")

    # 3. Map query genes to pathways
    reconstructed_data = []
    category_counts = {}
    pathway_completeness_raw = {}  # pid -> set of EC/KO numbers in genome
    
    for g in genes_with_ec:
        ecs = g["EC_Numbers"]
        kos = g["KO_Numbers"]
        
        # Populate KO_Numbers from EC_Numbers using our link map if empty
        if ecs and not kos:
            for ec in ecs:
                if ec in ec_to_ko:
                    for ko in ec_to_ko[ec]:
                        if ko not in kos:
                            kos.append(ko)
                            
        mapped_paths = []
        
        # Map via EC
        for ec in ecs:
            pids = ec_to_pathway.get(ec, [])
            for pid in pids:
                if pid != 'map01100':
                    mapped_paths.append((pid, pathway_names.get(pid, "Unknown Pathway")))
                    if pid not in pathway_completeness_raw:
                        pathway_completeness_raw[pid] = set()
                    pathway_completeness_raw[pid].add(ec)
                    
        # Map via KO
        for ko in kos:
            pids = ko_to_pathway.get(ko, [])
            for pid in pids:
                if pid != 'map01100':
                    mapped_paths.append((pid, pathway_names.get(pid, "Unknown Pathway")))
                    if pid not in pathway_completeness_raw:
                        pathway_completeness_raw[pid] = set()
                    pathway_completeness_raw[pid].add(ko)
                    
        if not mapped_paths:
            cat = "Other Metabolic Pathways"
            color = get_category_color(cat)
            reconstructed_data.append({
                **g,
                "Primary_Category": cat,
                "Primary_Pathway": "Unmapped",
                "Color_Hex": color,
                "All_Mapped_Pathways": "N/A"
            })
            category_counts[cat] = category_counts.get(cat, 0) + 1
            continue
            
        sorted_paths = sorted(mapped_paths, key=lambda x: pathway_to_category.get(x[0], "Other Metabolic Pathways"))
        prim_pid, prim_pname = sorted_paths[0]
        cat = pathway_to_category.get(prim_pid, "Other Metabolic Pathways")
        color = get_category_color(cat)
        all_paths_str = ", ".join([f"{pid} ({name})" for pid, name in mapped_paths])
        
        reconstructed_data.append({
            **g,
            "Primary_Category": cat,
            "Primary_Pathway": f"{prim_pid}: {prim_pname}",
            "Color_Hex": color,
            "All_Mapped_Pathways": all_paths_str
        })
        category_counts[cat] = category_counts.get(cat, 0) + 1
        
    df_metab = pd.DataFrame(reconstructed_data)
    os.makedirs(os.path.join(out_dir, "metabolic_pathways"), exist_ok=True)
    csv_out = os.path.join(out_dir, "metabolic_pathways", "kegg_enzymes.csv")
    df_metab.to_csv(csv_out, index=False)
    logger.info(f"Saved enzyme pathway mapping to: {csv_out}")
    
    # 4. Calculate Pathway Completeness & Topological/Gap-Filling Scores
    pathway_completeness_rows = []
    
    # Invert ec_to_pathway and ko_to_pathway to count total enzymes per pathway
    pathway_total_elements = {}
    for ec, pids in ec_to_pathway.items():
        for pid in pids:
            if pid not in pathway_total_elements:
                pathway_total_elements[pid] = set()
            pathway_total_elements[pid].add(ec)
    for ko, pids in ko_to_pathway.items():
        for pid in pids:
            if pid not in pathway_total_elements:
                pathway_total_elements[pid] = set()
            pathway_total_elements[pid].add(ko)
            
    for pid, name in pathway_names.items():
        if pid == 'map01100':
            continue
        total_set = pathway_total_elements.get(pid, set())
        if not total_set:
            continue
        present_set = pathway_completeness_raw.get(pid, set())
        present_count = len(present_set)
        total_count = len(total_set)
        
        enzyme_presence_score = (present_count / total_count) * 100.0 if total_count > 0 else 0.0
        
        # Pathway completeness is reported as detected enzyme fraction only.
        # Hole-filling heuristics are deliberately excluded: predicting missing enzymes are 'fillable'
        # based solely on pathway position inflates completeness scores without experimental basis.
        # For gap-filling with metabolic network topology, use GapSeq or ModelSEED.
        # Reference: Zimmermann et al. 2021 (GapSeq, doi:10.1186/s13059-021-02295-1)
        if present_count > 0:
            pathway_completeness_rows.append({
                "Pathway_ID": pid,
                "Pathway_Name": name,
                "Category": pathway_to_category.get(pid, "Other Metabolic Pathways"),
                "Total_Pathway_Elements": total_count,
                "Present_Elements": present_count,
                "Missing_Elements": total_count - present_count,
                "Enzyme_Presence_Score_Pct": round(enzyme_presence_score, 2),
                "Topological_Hole_Filling_Score_Pct": 0.0,
                "Completeness_Note": (
                    f"{present_count} of {total_count} pathway elements detected in genome annotation. "
                    "Missing enzymes may reflect annotation gaps, not true metabolic absence. "
                    "No gap-filling prediction applied."
                )
            })
            
    df_comp = pd.DataFrame(pathway_completeness_rows)
    df_comp = df_comp.sort_values(by="Enzyme_Presence_Score_Pct", ascending=False)
    comp_out = os.path.join(out_dir, "metabolic_pathways", "pathway_completeness.csv")
    df_comp.to_csv(comp_out, index=False)
    logger.info(f"Saved pathway completeness & topological scores to: {comp_out}")

    # 5. SBML level 3 export using cobra
    sbml_out = os.path.join(out_dir, "metabolic_pathways", "query_metabolism.xml")
    export_to_sbml(genes_with_ec, reconstructed_data, sbml_out, logger)

    if df_metab.empty:
        logger.warning("[WARNING] Reconstructed metabolism table is empty.")
        return []

    # 6. Pathway Distribution Plotting
    plt.figure(figsize=(10, 6))
    categories = list(category_counts.keys())
    counts = list(category_counts.values())
    
    # Establish a clean custom palette matching the categorization colors
    colors = [get_category_color(cat) for cat in categories]
    
    sns.set_theme(style="whitegrid")
    ax = sns.barplot(x=counts, y=categories, palette=colors, hue=categories, legend=False)
    plt.title(f"Reconstructed Metabolic Pathway Category Distribution: {query_org}", fontsize=12, fontweight="bold", pad=15)
    plt.xlabel("Number of Annotated Genes", fontsize=10)
    plt.ylabel("Pathway Group", fontsize=10)
    plt.tight_layout()
    
    fig_out = os.path.join(out_dir, "metabolic_pathways", "pathway_distribution.png")
    plt.savefig(fig_out, dpi=300)
    plt.close()
    logger.info(f"Saved metabolic pathway distribution bar plot to: {fig_out}")
    
    return reconstructed_data
