import os
import gzip
import shutil
import random
import ftplib
import urllib.request
import pandas as pd
from Bio import Entrez
import time

Entrez.email = "researcher@redolarium.edu"

BENCHMARK_DIR = "F:/IITk/Redolarium/redolarium_github/benchmarks/data"
GENOMES_DIR = os.path.join(BENCHMARK_DIR, "genomes")
os.makedirs(GENOMES_DIR, exist_ok=True)

def scrape_diverse_genomes(total_genomes=150, seed=42):
    print("Stage 1: Establishing Reproducible Ground Truth...")
    random.seed(seed)
    
    # Query for completed, reference prokaryotic genomes
    query = '("Bacteria"[Organism] OR "Archaea"[Organism]) AND "complete genome"[Assembly Level] AND "reference genome"[filter]'
    
    try:
        handle = Entrez.esearch(db="assembly", term=query, retmax=10000)
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"Failed to reach NCBI: {e}")
        return

    id_list = record["IdList"]
    print(f"Found {len(id_list)} highly curated representative genomes. Randomly sampling {total_genomes}...")
    
    sampled_ids = random.sample(id_list, min(total_genomes, len(id_list)))
    
    manifest = []
    
    # Download them
    for i, asm_id in enumerate(sampled_ids):
        try:
            sum_handle = Entrez.esummary(db="assembly", id=asm_id)
            summary = Entrez.read(sum_handle)
            sum_handle.close()
            
            doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
            acc = doc["AssemblyAccession"]
            org = doc["SpeciesName"]
            ftp_path = doc["FtpPath_RefSeq"]
            taxid = doc["Taxid"]
            
            if not ftp_path:
                continue
                
            file_prefix = ftp_path.split('/')[-1]
            download_url = f"{ftp_path}/{file_prefix}_genomic.gbff.gz"
            
            gz_path = os.path.join(GENOMES_DIR, f"{acc}.gbff.gz")
            gbk_path = os.path.join(GENOMES_DIR, f"{acc}.gbk")
            
            if not os.path.exists(gbk_path) or os.path.getsize(gbk_path) < 1000:
                if os.path.exists(gbk_path):
                    print(f"[{i+1}/{total_genomes}] {acc} appears truncated or corrupted. Re-downloading...")
                else:
                    print(f"[{i+1}/{total_genomes}] Downloading {acc} ({org})...")
                    
                # Use Request to provide a User-Agent
                req = urllib.request.Request(download_url, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(req) as response, open(gz_path, 'wb') as out_file:
                    shutil.copyfileobj(response, out_file)
                
                # Extract and delete to save space (9GB limit rule)
                with gzip.open(gz_path, 'rb') as f_in:
                    with open(gbk_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                        
                os.remove(gz_path)
            else:
                print(f"[{i+1}/{total_genomes}] {acc} already exists. Skipping download.")
                
            manifest.append({
                "Accession": acc,
                "Organism": org,
                "TaxID": taxid,
                "Assembly_ID": asm_id,
                "File_Path": gbk_path
            })
            
            time.sleep(0.5) # Be nice to NCBI
            
        except Exception as e:
            print(f"Error fetching {asm_id}: {e}")
            
    df = pd.DataFrame(manifest)
    manifest_path = os.path.join(BENCHMARK_DIR, "benchmark_manifest.csv")
    df.to_csv(manifest_path, index=False)
    print(f"\nGround Truth Manifest saved to: {manifest_path}")
    print("Scraping complete!")

if __name__ == "__main__":
    scrape_diverse_genomes(total_genomes=100, seed=42)
