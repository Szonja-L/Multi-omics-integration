#!/usr/bin/env python3
"""
run_pipeline.py
═══════════════════════════════════════════════════════════════════════════════
Master runner for the TCGA-GBM Multi-Omics Integration Pipeline.

Executes all five analysis scripts in order:

  01_simulate_data.py      — Simulate multi-omics data (GBM biology)
  02_mofa_integration.py   — Fit MOFA+ model, extract latent factors
  03_visualize.py          — Core MOFA+ figures (data overview, weights, UMAP)
  04_survival_analysis.py  — Kaplan-Meier & Cox regression with factor scores
  05_fusion_comparison.py  — Early / late / MOFA+ C-index comparison

Usage:
  python run_pipeline.py            # run all steps
  python run_pipeline.py --step 3   # re-run from step 3

Outputs:
  results/figures/   — 9 publication-quality PNGs
  results/tables/    — CSVs with factor scores, weights, survival stats
  results/mofa_model.hdf5
═══════════════════════════════════════════════════════════════════════════════
"""

import argparse
import pathlib
import subprocess
import sys
import time

BASE    = pathlib.Path(__file__).resolve().parent
SCRIPTS = BASE / "scripts"

STEPS = [
    (1, "01_simulate_data.py",     "Simulate biologically-realistic GBM multi-omics data"),
    (2, "02_mofa_integration.py",  "Fit MOFA+ model (mofapy2)"),
    (3, "03_visualize.py",         "Generate MOFA+ figures (Figs 1–5)"),
    (4, "04_survival_analysis.py", "Kaplan-Meier & Cox regression (Figs 6–7)"),
    (5, "05_fusion_comparison.py", "Early / late / MOFA+ fusion comparison (Figs 8–9)"),
]


def run_step(step_num: int, script_name: str, description: str) -> bool:
    script = SCRIPTS / script_name
    print(f"\n{'='*64}")
    print(f"  Step {step_num}/5 — {description}")
    print(f"  Script: {script.name}")
    print(f"{'='*64}")
    t0 = time.time()
    result = subprocess.run([sys.executable, str(script)],
                            capture_output=False)
    elapsed = time.time() - t0
    if result.returncode == 0:
        print(f"\n  ✓  Step {step_num} completed in {elapsed:.1f}s")
        return True
    else:
        print(f"\n  ✗  Step {step_num} FAILED (exit {result.returncode})")
        return False


def main():
    parser = argparse.ArgumentParser(description="Run TCGA-GBM multi-omics pipeline")
    parser.add_argument("--step", type=int, default=1, metavar="N",
                        help="Start from step N (1–5, default: 1)")
    parser.add_argument("--only", type=int, default=None, metavar="N",
                        help="Run only step N")
    args = parser.parse_args()

    print("\n" + "█"*64)
    print("  TCGA-GBM Multi-Omics Integration Pipeline")
    print("  MOFA+ · lifelines · scikit-learn")
    print("█"*64)

    start = time.time()
    failed = []

    for num, script, desc in STEPS:
        if args.only is not None:
            if num != args.only:
                continue
        elif num < args.step:
            print(f"  ⊘  Skipping step {num} (--step={args.step})")
            continue

        ok = run_step(num, script, desc)
        if not ok:
            failed.append(num)
            ans = input(f"\n  Continue despite failure in step {num}? [y/N] ").strip().lower()
            if ans != "y":
                break

    total = time.time() - start
    print(f"\n{'='*64}")
    if not failed:
        print(f"  ✓  Pipeline complete in {total:.0f}s")
        print(f"  Results: {BASE / 'results'}")
    else:
        print(f"  ✗  Pipeline finished with failures in steps: {failed}")
    print("=" * 64 + "\n")


if __name__ == "__main__":
    main()
