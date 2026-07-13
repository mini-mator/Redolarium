#!/usr/bin/env python3
import os
import sys
import json
import tarfile
import urllib.request
import csv
import time
from urllib.error import URLError

TARGET_COUNT = 100
MIBIG_URL = "https://dl.secondarymetabolites.org/mibig/mibig_json_3.1.tar.gz"
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
GENOMES_DIR = os.path.join(DATA_DIR, "genomes")
JSON_DIR = os.path.join(DATA_DIR, "mibig_jsons")

os.makedirs(GENOMES_DIR, exist_ok=True)
os.makedirs(JSON_DIR, exist_ok=True)

def download_mibig():
    tar_path = os.path.join(DATA_DIR, "mibig_json_3.1.tar.gz")
    if not os.path.exists(tar_path):
        print(f"Downloading MIBiG 3.1 database from {MIBIG_URL}...")
        try:
            urllib.request.urlretrieve(MIBIG_URL, tar_path)
        except URLError as e:
            print(f"Network Error: Could not download MIBiG ({e}). Please ensure your network allows access to dl.secondarymetabolites.org.")
            print("To manually bypass, download the JSON tarball from MIBiG and place it at: " + tar_path)
            sys.exit(1)
            
    print("Extracting MIBiG JSONs...")
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(path=JSON_DIR)

def check_accession_length(acc):
    """Query NCBI E-utilities to get the sequence length."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id={acc}&retmode=json"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Redolarium-Benchmark/1.0'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode())
            res = data.get("result", {})
            uids = res.get("uids", [])
            if not uids: return 0
            return int(res[uids[0]].get("slen", 0))
    except Exception:
        return 0

def fetch_genbank(acc, out_path):
    """Download full GenBank record."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={acc}&rettype=gbwithparts&retmode=text"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Redolarium-Benchmark/1.0'})
        with urllib.request.urlopen(req) as response:
            with open(out_path, 'wb') as f:
                f.write(response.read())
        return True
    except Exception as e:
        print(f"Failed to fetch {acc}: {e}")
        return False

def harvest_genomes():
    download_mibig()
    
    # Locate the extracted json dir (sometimes it extracts into a subfolder)
    subdirs = [os.path.join(JSON_DIR, d) for d in os.listdir(JSON_DIR) if os.path.isdir(os.path.join(JSON_DIR, d))]
    target_json_dir = subdirs[0] if subdirs else JSON_DIR
    
    manifest = []
    collected_accs = set()
    
    print("Scanning MIBiG JSONs for suitable complete genomes...")
    for filename in os.listdir(target_json_dir):
        if not filename.endswith(".json"): continue
        if len(collected_accs) >= TARGET_COUNT: break
        
        filepath = os.path.join(target_json_dir, filename)
        with open(filepath, 'r') as f:
            try:
                data = json.load(f)
            except json.JSONDecodeError:
                continue
                
        cluster = data.get("cluster", {})
        mibig_id = cluster.get("mibig_accession", filename.split(".")[0])
        loci = cluster.get("loci", {})
        if not loci: continue
        
        acc = loci.get("accession")
        start = loci.get("start_coord", -1)
        end = loci.get("end_coord", -1)
        classes = cluster.get("biosyn_class", [])
        
        if not acc or start == -1 or end == -1: continue
        
        # We want unique genomes that are actually full chromosomes
        # MIBiG often contains fragments. We filter for >1Mbp.
        if acc in collected_accs:
            # We already have this genome, just add the BGC to manifest
            manifest.append({
                "Genome_Accession": acc,
                "MIBiG_ID": mibig_id,
                "Start_Coord": start,
                "End_Coord": end,
                "BGC_Class": "|".join(classes)
            })
            continue
            
        # Check length (rate limited to 3/sec without API key)
        length = check_accession_length(acc)
        time.sleep(0.35)
        
        if length > 1000000:
            print(f"Found Golden Genome: {acc} ({length:,} bp) containing {mibig_id}")
            out_gbk = os.path.join(GENOMES_DIR, f"{acc}.gbk")
            if not os.path.exists(out_gbk):
                success = fetch_genbank(acc, out_gbk)
                if not success: continue
                time.sleep(0.35)
                
            collected_accs.add(acc)
            manifest.append({
                "Genome_Accession": acc,
                "MIBiG_ID": mibig_id,
                "Start_Coord": start,
                "End_Coord": end,
                "BGC_Class": "|".join(classes)
            })
            
    print(f"\nSuccessfully harvested {len(collected_accs)} genomes containing {len(manifest)} verified BGCs.")
    
    manifest_path = os.path.join(DATA_DIR, "ground_truth_manifest.csv")
    ledger_path = os.path.join(DATA_DIR, "accession_ledger.csv")
    
    with open(manifest_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Genome_Accession", "MIBiG_ID", "Start_Coord", "End_Coord", "BGC_Class"])
        writer.writeheader()
        writer.writerows(manifest)
        
    with open(ledger_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Genome_Accession"])
        for acc in collected_accs:
            writer.writerow([acc])
            
    print(f"Manifest written to: {manifest_path}")
    print(f"Ledger written to: {ledger_path}")
    print("Phase 1 Complete. Ready for run_validation_suite.py.")

if __name__ == "__main__":
    harvest_genomes()
