#!/usr/bin/env python3
"""
02_mofa_integration.py
========================
Run MOFA+ (Multi-Omics Factor Analysis v2) on the simulated TCGA-GBM data.

MOFA+ learns a latent factor model:
    X[m]  ≈  Z · W[m]ᵀ    (for each view m)

where  Z  (N × K) are sample-level factor scores shared across all views
and   W[m] (D_m × K) are view-specific feature weights.

The model is fitted via variational Bayes with:
  • ARD (Automatic Relevance Determination) priors on W → per-view sparsity
  • Spike-and-slab priors on W → feature-level sparsity
  • Gaussian likelihood for expression and methylation
  • Bernoulli likelihood for binary mutation data

Outputs (saved to results/tables/):
  factor_scores.csv          – (N × K) MOFA+ factor scores per sample
  factor_weights_<view>.csv  – (D_m × K) feature weights per view
  variance_explained.csv     – (K × 3) R² per factor per view
"""

import numpy as np
import pandas as pd
import h5py
import os
import sys
import warnings
warnings.filterwarnings("ignore")

from mofapy2.run.entry_point import entry_point

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA  = os.path.join(BASE, "data", "simulated")
TABS  = os.path.join(BASE, "results", "tables")
MODEL = os.path.join(BASE, "results", "mofa_model.hdf5")
os.makedirs(TABS, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════════
#  MOFA+ TRAINING
# ═══════════════════════════════════════════════════════════════════════════════

def load_data():
    """Load simulated multi-omics data matrices."""
    expr = pd.read_csv(os.path.join(DATA, "expression_matrix.csv"), index_col=0)
    meth = pd.read_csv(os.path.join(DATA, "methylation_matrix.csv"), index_col=0)
    muts = pd.read_csv(os.path.join(DATA, "mutation_matrix.csv"),    index_col=0)
    clin = pd.read_csv(os.path.join(DATA, "clinical_data.csv"))
    print(f"  Loaded: expr {expr.shape}, meth {meth.shape}, mut {muts.shape}")
    return expr, meth, muts, clin


def run_mofa(expr: pd.DataFrame, meth: pd.DataFrame,
             muts: pd.DataFrame, n_factors: int = 10) -> str:
    """
    Train MOFA+ and save model to HDF5.

    Parameters
    ----------
    n_factors : int
        Maximum number of latent factors (inactive ones are pruned).

    Returns
    -------
    str  – path to saved HDF5 model file
    """
    print(f"\n  Initialising MOFA+ with K={n_factors} factors…")

    # MOFA+ expects: data[views][groups], each array (n_samples × n_features)
    data_mats = [
        [expr.values.astype(np.float64)],   # view 0: expression
        [meth.values.astype(np.float64)],   # view 1: methylation
        [muts.values.astype(np.float64)],   # view 2: mutations
    ]

    ent = entry_point()

    ent.set_data_matrix(
        data_mats,
        likelihoods=["gaussian", "gaussian", "bernoulli"],
        views_names=["expression", "methylation", "mutations"],
        groups_names=["TCGA-GBM"],
    )

    ent.set_data_options(scale_views=True)

    ent.set_model_options(
        factors=n_factors,
        spikeslab_weights=True,     # feature-level sparsity
        ard_weights=True,           # view+factor-level ARD (key for modality attribution)
        spikeslab_factors=False,
        ard_factors=False,
    )

    ent.set_train_options(
        iter=1500,
        convergence_mode="medium",  # good balance of speed vs quality
        dropR2=0.001,               # drop factors explaining <0.1% total R²
        seed=42,
        verbose=False,
    )

    ent.build()
    ent.run()
    ent.save(MODEL)

    print(f"\n  ✓  Model saved → {MODEL}")
    return MODEL


# ═══════════════════════════════════════════════════════════════════════════════
#  EXTRACT RESULTS
# ═══════════════════════════════════════════════════════════════════════════════

def extract_results(model_path: str,
                    expr: pd.DataFrame,
                    meth: pd.DataFrame,
                    muts: pd.DataFrame,
                    clin: pd.DataFrame) -> dict:
    """
    Parse the MOFA+ HDF5 model and return tidy DataFrames.

    H5 structure (confirmed for mofapy2 v0.7.4):
      expectations/Z/TCGA-GBM     (K_active, N)
      expectations/W/expression   (K_active, D_expr)
      variance_explained/r2_per_factor/TCGA-GBM  (3_views, K_active)
      variance_explained/r2_total/TCGA-GBM       (3_views,)

    Returns
    -------
    dict with keys: Z, W_expr, W_meth, W_muts, r2, r2_tot
    """
    GROUP = "TCGA-GBM"
    f     = h5py.File(model_path, "r")

    # ── Factor scores Z : (K, N) → (N, K) ────────────────────────────────
    Z_raw = f["expectations"]["Z"][GROUP][:]       # (K, N)
    K, N  = Z_raw.shape
    Z     = Z_raw.T                                # (N, K)

    factor_names = [f"Factor{k+1}" for k in range(K)]
    sample_ids   = expr.index.tolist()
    Z_df = pd.DataFrame(Z, index=sample_ids, columns=factor_names)

    # ── Feature weights W : (K, D) → (D, K) ──────────────────────────────
    def get_W(view_key, feature_names):
        W_raw = f["expectations"]["W"][view_key][:]   # (K, D)
        return pd.DataFrame(W_raw.T, index=feature_names, columns=factor_names)

    W_expr = get_W("expression", expr.columns.tolist())
    W_meth = get_W("methylation", meth.columns.tolist())
    W_muts = get_W("mutations",   muts.columns.tolist())

    # ── Variance explained ─────────────────────────────────────────────────
    # r2_per_factor[GROUP] : (3 views, K)  rows=[expr,meth,mut]
    r2_raw = f["variance_explained"]["r2_per_factor"][GROUP][:]   # (3, K)
    r2_df  = pd.DataFrame(
        r2_raw.T,                                   # (K, 3)
        index=factor_names,
        columns=["Expression", "Methylation", "Mutations"],
    )

    # r2_total[GROUP] : (3,) total variance explained per view
    r2_total = f["variance_explained"]["r2_total"][GROUP][:]      # (3,)
    r2_tot   = pd.Series(r2_total,
                         index=["Expression", "Methylation", "Mutations"])

    f.close()

    return {
        "Z":      Z_df,
        "W_expr": W_expr,
        "W_meth": W_meth,
        "W_muts": W_muts,
        "r2":     r2_df,
        "r2_tot": r2_tot,
    }


def save_results(results: dict) -> None:
    """Save extracted results to results/tables/."""
    results["Z"].to_csv(os.path.join(TABS, "factor_scores.csv"))
    results["W_expr"].to_csv(os.path.join(TABS, "weights_expression.csv"))
    results["W_meth"].to_csv(os.path.join(TABS, "weights_methylation.csv"))
    results["W_muts"].to_csv(os.path.join(TABS, "weights_mutations.csv"))
    results["r2"].to_csv(os.path.join(TABS, "variance_explained.csv"))
    results["r2_tot"].to_csv(os.path.join(TABS, "variance_explained_total.csv"))
    print("  ✓  Results saved to results/tables/")


def print_summary(results: dict) -> None:
    """Print MOFA+ summary to stdout."""
    Z   = results["Z"]
    r2  = results["r2"]
    K   = len(Z.columns)

    print(f"\n  ── MOFA+ Summary ─────────────────────────────────────")
    print(f"  Active factors   : {K}")
    print(f"  Samples          : {len(Z)}")
    print(f"\n  Marginal variance explained per factor (%):")
    header = f"  {'Factor':<12}  {'Expression':>12}  {'Methylation':>12}  {'Mutations':>10}"
    print(header)
    print("  " + "─" * 52)
    for fac in Z.columns:
        e = r2.loc[fac, "Expression"]      # already in percent
        m = r2.loc[fac, "Methylation"]
        u = r2.loc[fac, "Mutations"]
        print(f"  {fac:<12}  {e:>11.2f}%  {m:>11.2f}%  {u:>9.2f}%")

    print(f"\n  Total R² explained by all factors (%):")
    for view, r2v in results["r2_tot"].items():
        print(f"    {view:<14}: {r2v:.1f}%")


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("  MOFA+ Multi-Omics Integration  –  TCGA-GBM")
    print("=" * 60)

    expr, meth, muts, clin = load_data()

    model_path = run_mofa(expr, meth, muts, n_factors=10)

    print("\n  Extracting factor scores and weights…")
    results = extract_results(model_path, expr, meth, muts, clin)

    print_summary(results)
    save_results(results)

    print("\n  ✓  MOFA+ integration complete.")
    return results


if __name__ == "__main__":
    main()
