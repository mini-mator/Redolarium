import os
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import f1_score

# =====================================================================
# Redolarium Empirical Weight Optimization Benchmark
# =====================================================================

def main():
    os.makedirs("benchmarks/results", exist_ok=True)
    print("Initiating Grid Search for Scoring Weights against validation dataset...")

    dataset_path = "benchmarks/data/validation_alignments.csv"
    if not os.path.exists(dataset_path):
        print(f"Error: Validation dataset not found at {dataset_path}.")
        print("Cannot perform weight optimization without empirical data.")
        print("Please provide a CSV with columns: identity, coverage, bitscore, label (1 for TP, 0 for FP).")
        return

    df = pd.read_csv(dataset_path)
    if not all(col in df.columns for col in ['identity', 'coverage', 'bitscore', 'label']):
        print("Error: Dataset missing required columns.")
        return

    X_id = df['identity'].values
    X_cov = df['coverage'].values
    X_bit = df['bitscore'].values
    y_true = df['label'].values

    # Grid search combinations (increments of 0.1)
    weights = np.arange(0.0, 1.1, 0.1)
    combinations = [
        (round(w1, 1), round(w2, 1), round(w3, 1))
        for w1 in weights for w2 in weights for w3 in weights
        if np.isclose(w1 + w2 + w3, 1.0)
    ]

    best_f1 = 0
    best_w = None
    results = []

    for (w1, w2, w3) in combinations:
        # Composite score
        scores = (w1 * X_id) + (w2 * X_cov) + (w3 * X_bit)
        
        # Optimal decision boundary for this distribution
        y_pred = (scores > 0.65).astype(int)
        f1 = f1_score(y_true, y_pred, zero_division=0)
        
        results.append((w1, w2, w3, f1))
        
        if f1 > best_f1:
            best_f1 = f1
            best_w = (w1, w2, w3)

    if best_w:
        print(f"Optimal Empirical Weights -> Identity: {best_w[0]}, Coverage: {best_w[1]}, Bitscore: {best_w[2]}")
        print(f"Maximum F1-Score: {best_f1:.4f}")

        # Plotting
        f1_scores = [r[3] for r in results]
        ids = [r[0] for r in results]
        covs = [r[1] for r in results]
        
        plt.figure(figsize=(10, 7))
        sc = plt.scatter(ids, covs, c=f1_scores, cmap='viridis', s=100, edgecolor='k', alpha=0.8)
        plt.colorbar(sc, label="F1-Score (Empirical Validation)")
        plt.xlabel("Weight: Identity ($w_1$)")
        plt.ylabel("Weight: Coverage ($w_2$)")
        plt.title("Reference Dataset: Grid-Search Weight Optimization\n($w_1 + w_2 + w_3 = 1$)")
        
        plt.scatter([best_w[0]], [best_w[1]], c='red', marker='*', s=300, edgecolor='black', label=f'Optimal {best_w}')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.5)
        
        out_png = "benchmarks/results/weight_optimization_f1.png"
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        print(f"Scatter plot saved to: {out_png}")

if __name__ == "__main__":
    main()
