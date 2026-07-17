#!/usr/bin/env python3
import os
import sys
import pandas as pd
import zipfile
import shutil

BASE_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(BASE_DIR, "data")
MANIFEST = os.path.join(DATA_DIR, "ground_truth_manifest.csv")
RESULTS_DIR = os.path.join(DATA_DIR, "results")
BUNDLE_DIR = os.path.join(BASE_DIR, "Reproducibility_Bundle")
ANALYTICS_DIR = os.path.join(BUNDLE_DIR, "analytics")

def evaluate_predictions():
    print("Evaluating Redolarium Predictions vs MIBiG Ground Truth...")
    
    if not os.path.exists(MANIFEST):
        print("Manifest not found.")
        return
        
    df_manifest = pd.read_csv(MANIFEST)
    
    total_mibig_clusters = len(df_manifest)
    true_positives = 0
    false_negatives = 0
    
    validation_records = []
    
    for _, row in df_manifest.iterrows():
        acc = row["Genome_Accession"]
        mibig_id = row["MIBiG_ID"]
        category = row.get("Category", "Unknown")
        gt_start = int(row["Start_Coord"])
        gt_end = int(row["End_Coord"])
        gt_class = row["BGC_Class"]
        
        summary_file = os.path.join(RESULTS_DIR, category, acc, "bgc", "BGC_Comprehensive_Summary.xlsx")
        
        detected = False
        pred_id = "None"
        overlap_pct = 0.0
        
        if os.path.exists(summary_file):
            try:
                df_pred = pd.read_excel(summary_file, sheet_name="BGC_Summary")
                for _, p_row in df_pred.iterrows():
                    p_start = int(p_row.get("Start_Coord", -1))
                    p_end = int(p_row.get("End_Coord", -1))
                    
                    if p_start == -1 or p_end == -1: continue
                    
                    # Calculate spatial overlap
                    overlap_start = max(gt_start, p_start)
                    overlap_end = min(gt_end, p_end)
                    
                    if overlap_end > overlap_start:
                        overlap_len = overlap_end - overlap_start
                        gt_len = gt_end - gt_start
                        overlap_pct = (overlap_len / gt_len) * 100
                        
                        # If Redolarium's cluster overlaps the MIBiG cluster by >20%, it's a hit.
                        if overlap_pct > 20.0:
                            detected = True
                            pred_id = p_row.get("BGC_ID", "Unknown")
                            break
            except Exception as e:
                pass
                
        if detected:
            true_positives += 1
        else:
            false_negatives += 1
            
        validation_records.append({
            "Genome": acc,
            "MIBiG_ID": mibig_id,
            "Ground_Truth_Class": gt_class,
            "Category": category,
            "Status": "True Positive" if detected else "False Negative",
            "Predicted_BGC_ID": pred_id,
            "Overlap_Percent": f"{overlap_pct:.1f}%"
        })
        
    os.makedirs(ANALYTICS_DIR, exist_ok=True)
    df_val = pd.DataFrame(validation_records)
    df_val.to_csv(os.path.join(ANALYTICS_DIR, "Validation_Matrix.csv"), index=False)
    
    print("\n" + "="*50)
    print(" GROUND TRUTH VALIDATION RESULTS (Phase 3)")
    print("="*50)
    
    if not df_val.empty:
        for cat in df_val["Category"].unique():
            cat_df = df_val[df_val["Category"] == cat]
            tot = len(cat_df)
            tp = len(cat_df[cat_df["Status"] == "True Positive"])
            fn = tot - tp
            recall = (tp / tot) * 100 if tot > 0 else 0
            print(f"[{cat}] Total: {tot} | True Positives: {tp} | Missed: {fn} | Recall: {recall:.2f}%")
            
    print("="*50)
    print("NOTE: Precision requires manual curation of 'novel' predictions which MIBiG lacks.")
    print("="*50 + "\n")

def build_reproducibility_bundle():
    print("Assembling Journal-Compliant Reproducibility Bundle (Phase 4)...")
    
    os.makedirs(BUNDLE_DIR, exist_ok=True)
    
    # 1. Accession Ledger
    ledger_src = os.path.join(DATA_DIR, "accession_ledger.csv")
    if os.path.exists(ledger_src):
        shutil.copy(ledger_src, os.path.join(BUNDLE_DIR, "Accession_Ledger.csv"))
        
    # 2. Version Provenance
    provenance_path = os.path.join(BUNDLE_DIR, "methods_provenance.txt")
    with open(provenance_path, "w") as f:
        f.write("REDOLARIUM PIPELINE PROVENANCE LOG\n")
        f.write("==================================\n")
        f.write("Redolarium Pipeline: v1.1.0 (Strict HMM Branch)\n")
        f.write("antiSMASH Version: 7.0.0 (via Docker)\n")
        f.write("PyHMMER Version: 0.10.4\n")
        f.write("Ground Truth Database: MIBiG JSON 3.1\n")
        
    # 3. Zip it all up
    print("Compressing bundle into Reproducibility_Bundle.zip...")
    zip_path = os.path.join(BASE_DIR, "Reproducibility_Bundle.zip")
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(BUNDLE_DIR):
            for file in files:
                abs_path = os.path.join(root, file)
                rel_path = os.path.relpath(abs_path, BASE_DIR)
                zipf.write(abs_path, rel_path)
                
    print(f"Done! Submittable archive ready at: {zip_path}")

if __name__ == "__main__":
    evaluate_predictions()
    build_reproducibility_bundle()
