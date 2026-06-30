#!/usr/bin/env python3
"""
patch_scripts.py
================
One-time patch that updates hardcoded labels and cosmetic strings in
scripts/03_visualize.py and scripts/04_survival_analysis.py to reflect
the real TCGA-GBM MOFA+ factor structure:

  Factor1 → Mutation landscape      (39% mutation variance)
  Factor2 → Transcriptional subtype (22% expression variance)
  Factor3 → Expression module A     (13% expression variance)
  Factor4 → Expression module B     (9%  expression variance)
  Factor5 → Methylation axis        (9%  methylation variance)

Run ONCE from the repo root, then run 03, 04, 05 as normal.

Usage:
    cd "C:\\Users\\Szonja\\Desktop\\GitHub\\8. Multi-omics-integration"
    python gbm-real-data\\patch_scripts.py
"""

import os
import re

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
# patch_scripts.py lives in gbm-real-data/, so repo root is one level up
REPO_ROOT = os.path.dirname(REPO_ROOT)

VIZ  = os.path.join(REPO_ROOT, "scripts", "03_visualize.py")
SURV = os.path.join(REPO_ROOT, "scripts", "04_survival_analysis.py")


def patch(filepath, replacements, label):
    with open(filepath, "r", encoding="utf-8") as f:
        src = f.read()

    original = src
    for old, new in replacements:
        if old not in src:
            print(f"  !! [{label}] Could not find: {repr(old[:80])}")
        else:
            src = src.replace(old, new, 1)  # replace only the first occurrence
            print(f"  ✓  [{label}] Patched: {repr(old[:60])} …")

    if src != original:
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(src)
        print(f"  → Saved {filepath}\n")
    else:
        print(f"  (no changes written to {filepath})\n")


# ── 03_visualize.py patches ───────────────────────────────────────────────────

viz_patches = [

    # 1. Figure 1 title: hardcoded N=150 → real N from data shape
    (
        '"TCGA-GBM Multi-Omics Data Overview  (N=150)", y=1.07,',
        'f"TCGA-GBM Multi-Omics Data Overview  (N={expr.shape[0]})", y=1.07,',
    ),

    # 2. Expression xlabel: log2 TPM → log2 RSEM+1 (real RNA-seq RSEM data)
    (
        '"Gene Expression\\n(log\\u2082 TPM)"',
        '"Gene Expression\\n(log\\u2082 RSEM+1)"',
    ),

    # 3. Methylation xlabel: CpG sites → gene-level (API returns gene-level means)
    (
        'f"{meth.shape[1]} CpG sites  (80 shown)"',
        'f"{meth.shape[1]} genes (HM450, gene-level)"',
    ),

    # 4. Figure 2: variance-explained factor annotations
    (
        '''    factor_annot = {
        "Factor1": "Proneural/Mes\\naxis",
        "Factor2": "Mutation\\nsignature",
        "Factor3": "G-CIMP /\\nMGMT",
        "Factor4": "Expression\\nmodule",
        "Factor5": "Minor\\nsignal",
    }''',
        '''    factor_annot = {
        "Factor1": "Mutation\\nlandscape",
        "Factor2": "Transcriptional\\nsubtype",
        "Factor3": "Expression\\nmodule A",
        "Factor4": "Expression\\nmodule B",
        "Factor5": "Methylation\\naxis",
    }''',
    ),

    # 5. Figure 3: UMAP right-panel title (Factor1 label)
    (
        '"MOFA+ Factor Space \\u2013 Factor 1 (Proneural \\u2194 Mesenchymal)"',
        '"MOFA+ Factor Space \\u2013 Factor 1 (Mutation landscape)"',
    ),

    # 6. Figure 4: violin-plot factor titles
    (
        '''    factor_titles = [
        "Factor 1\\nProneural \\u2194 Mes",
        "Factor 2\\nMutation signature",
        "Factor 3\\nG-CIMP / MGMT",
        "Factor 4\\nExpression module",
        "Factor 5\\nMinor signal",
    ]''',
        '''    factor_titles = [
        "Factor 1\\nMutation landscape",
        "Factor 2\\nTranscriptional\\nsubtype",
        "Factor 3\\nExpression\\nmodule A",
        "Factor 4\\nExpression\\nmodule B",
        "Factor 5\\nMethylation\\naxis",
    ]''',
    ),

    # 7. Figure 5: top-weights configs — swap Factor3 row for Factor2
    (
        '        ("Factor3", W_e, "Expression",   VIEW_COL["Expression"],  16, 1, 0),\n'
        '        ("Factor3", W_m, "Methylation",  VIEW_COL["Methylation"], 16, 1, 1),\n'
        '        ("Factor3", W_u, "Mutations",    VIEW_COL["Mutations"],   15, 1, 2),',
        '        ("Factor2", W_e, "Expression",   VIEW_COL["Expression"],  16, 1, 0),\n'
        '        ("Factor2", W_m, "Methylation",  VIEW_COL["Methylation"], 16, 1, 1),\n'
        '        ("Factor2", W_u, "Mutations",    VIEW_COL["Mutations"],   15, 1, 2),',
    ),

    # 8. Figure 5: factor_labels dict in top-weights (Factor3 → Factor2, new descriptions)
    (
        '''    factor_labels = {
        "Factor1": "Factor 1  (Proneural \\u2194 Mesenchymal axis)",
        "Factor3": "Factor 3  (G-CIMP / MGMT methylation)",
    }''',
        '''    factor_labels = {
        "Factor1": "Factor 1  (Mutation landscape)",
        "Factor2": "Factor 2  (Transcriptional subtype axis)",
    }''',
    ),
]


# ── 04_survival_analysis.py patches ──────────────────────────────────────────

surv_patches = [

    # 9. Cox forest plot: factor_labels dict
    (
        '''    factor_labels = {
        "Factor1": "Factor 1\\n(Proneural/Mes axis)",
        "Factor2": "Factor 2\\n(Mutation signature)",
        "Factor3": "Factor 3\\n(G-CIMP/MGMT)",
        "Factor4": "Factor 4\\n(Expression module)",
        "Factor5": "Factor 5\\n(Minor signal)",
    }''',
        '''    factor_labels = {
        "Factor1": "Factor 1\\n(Mutation landscape)",
        "Factor2": "Factor 2\\n(Transcriptional subtype)",
        "Factor3": "Factor 3\\n(Expression module A)",
        "Factor4": "Factor 4\\n(Expression module B)",
        "Factor5": "Factor 5\\n(Methylation axis)",
    }''',
    ),
]


if __name__ == "__main__":
    print("=" * 60)
    print("  Patching scripts for real TCGA-GBM data")
    print("=" * 60)
    print(f"\nRepo root: {REPO_ROOT}\n")

    print("── 03_visualize.py ──────────────────────────────────")
    patch(VIZ, viz_patches, "03")

    print("── 04_survival_analysis.py ──────────────────────────")
    patch(SURV, surv_patches, "04")

    print("=" * 60)
    print("Done. Now run:")
    print("  python scripts\\03_visualize.py")
    print("  python scripts\\04_survival_analysis.py")
    print("  python scripts\\05_fusion_comparison.py")
    print("=" * 60)
