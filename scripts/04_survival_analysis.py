#!/usr/bin/env python3
"""
04_survival_analysis.py
========================
Link MOFA+ factors to patient survival using lifelines.

Analyses:
  1. Kaplan-Meier curves – patients stratified by Factor1 and Factor2 scores
  2. Log-rank tests – significance of factor-based stratification
  3. Cox proportional-hazards regression – MOFA+ factors as covariates
  4. Forest plot of hazard ratios

Figures produced:
  06_km_curves.png       – KM curves for top factors
  07_cox_forest.png      – Cox regression forest plot
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings, os
warnings.filterwarnings("ignore")

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from scipy import stats

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "savefig.facecolor": "white",
})

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data", "simulated")
TABS = os.path.join(BASE, "results", "tables")
FIGS = os.path.join(BASE, "results", "figures")

# ── Palette ───────────────────────────────────────────────────────────────────
HI_COL = {"high": "#E15759", "low": "#4E79A7"}
SUB_COL = {"PN": "#4E79A7", "CL": "#F28E2B", "MES": "#E15759"}
SUB_FULL = {"PN": "Proneural", "CL": "Classical", "MES": "Mesenchymal"}


def load_data():
    clin = pd.read_csv(os.path.join(DATA, "clinical_data.csv"))
    Z    = pd.read_csv(os.path.join(TABS, "factor_scores.csv"), index_col=0)
    # Align on sample order
    clin.index = Z.index
    merged = pd.concat([clin.reset_index(drop=True), Z.reset_index(drop=True)], axis=1)
    return merged


# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE 6 – KAPLAN-MEIER CURVES
# ═══════════════════════════════════════════════════════════════════════════════

def fig_km_curves(df: pd.DataFrame) -> dict:
    """
    KM curves for:
      Panel A – By molecular subtype (reference)
      Panel B – Factor1 high vs low (median split)
      Panel C – Factor2 high vs low
      Panel D – Factor3 high vs low (MGMT axis)
    """
    print("  → Figure 6: Kaplan-Meier survival curves…")

    fig, axes = plt.subplots(1, 4, figsize=(18, 5.5))
    lr_results = {}

    # ── Panel A: by subtype ────────────────────────────────────────────────
    ax = axes[0]
    kmf = KaplanMeierFitter()
    for st in ["PN", "CL", "MES"]:
        mask = df["subtype"] == st
        kmf.fit(df.loc[mask, "os_months"],
                event_observed=df.loc[mask, "os_event"],
                label=SUB_FULL[st])
        kmf.plot_survival_function(ax=ax, ci_show=True,
                                   color=SUB_COL[st], linewidth=2)

    # Multivariate log-rank
    res_st = multivariate_logrank_test(
        df["os_months"], df["subtype"], df["os_event"])
    ax.set_title(f"By Molecular Subtype\nLog-rank p = {res_st.p_value:.2e}",
                 fontsize=9.5, fontweight="bold")
    ax.set_xlabel("Time (months)", fontsize=9)
    ax.set_ylabel("Survival probability", fontsize=9)
    ax.set_ylim(-0.05, 1.10)
    ax.legend(loc="upper right", fontsize=7.5, frameon=False)
    ax.text(0.05, 0.06, f"N = {len(df)}", transform=ax.transAxes,
            fontsize=8, color="#666666")
    lr_results["subtype"] = res_st.p_value

    # ── Panels B–D: Factor high vs low ────────────────────────────────────
    factor_panels = [
        ("Factor1", "Factor 1  (Proneural ↔ Mes)",
         "↑ Proneural (high)", "↑ Mesenchymal (low)"),
        ("Factor2", "Factor 2  (Mutation signature)",
         "High Factor 2",      "Low Factor 2"),
        ("Factor3", "Factor 3  (G-CIMP / MGMT)",
         "MGMT-methylated",    "MGMT-unmethylated"),
    ]

    for i, (fac, title, lab_hi, lab_lo) in enumerate(factor_panels):
        ax = axes[i + 1]
        if fac not in df.columns:
            ax.set_visible(False)
            continue

        median  = df[fac].median()
        hi_mask = df[fac] >= median
        lo_mask = ~hi_mask

        kmf_hi = KaplanMeierFitter()
        kmf_lo = KaplanMeierFitter()
        kmf_hi.fit(df.loc[hi_mask, "os_months"], df.loc[hi_mask, "os_event"],
                   label=lab_hi)
        kmf_lo.fit(df.loc[lo_mask, "os_months"], df.loc[lo_mask, "os_event"],
                   label=lab_lo)
        kmf_hi.plot_survival_function(ax=ax, ci_show=True,
                                      color=HI_COL["high"], linewidth=2)
        kmf_lo.plot_survival_function(ax=ax, ci_show=True,
                                      color=HI_COL["low"],  linewidth=2)

        lr = logrank_test(
            df.loc[hi_mask, "os_months"], df.loc[lo_mask, "os_months"],
            df.loc[hi_mask, "os_event"],  df.loc[lo_mask, "os_event"],
        )
        p = lr.p_value
        lr_results[fac] = p

        # Significance stars
        stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        ax.set_title(f"{title}\nLog-rank p = {p:.2e}  {stars}",
                     fontsize=9.5, fontweight="bold")
        ax.set_xlabel("Time (months)", fontsize=9)
        ax.set_ylabel("", fontsize=9)
        ax.set_ylim(-0.05, 1.10)
        ax.legend(loc="upper right", fontsize=7.5, frameon=False)

        # Median survival lines
        m_hi = kmf_hi.median_survival_time_
        m_lo = kmf_lo.median_survival_time_
        ax.axhline(0.5, ls=":", lw=0.8, c="#aaaaaa")
        ax.text(0.97, 0.45, f"Med. Hi={m_hi:.1f}m / Lo={m_lo:.1f}m",
                transform=ax.transAxes, fontsize=7.5, ha="right",
                color="#555555")

    fig.suptitle("Survival Analysis: MOFA+ Factors Predict GBM Patient Outcomes",
                 fontsize=12, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGS, "06_km_curves.png"))
    plt.close(fig)
    return lr_results


# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE 7 – COX REGRESSION FOREST PLOT
# ═══════════════════════════════════════════════════════════════════════════════

def fig_cox_forest(df: pd.DataFrame) -> pd.DataFrame:
    """Fit univariate and multivariate Cox PH and plot forest of HRs."""
    print("  → Figure 7: Cox regression & forest plot…")

    # ── Univariate Cox ─────────────────────────────────────────────────────
    factor_cols = [c for c in df.columns if c.startswith("Factor")]
    rows        = []
    for fac in factor_cols:
        tmp = df[["os_months", "os_event", fac]].copy()
        tmp = tmp.rename(columns={fac: "feature"})
        cph = CoxPHFitter()
        cph.fit(tmp, duration_col="os_months", event_col="os_event",
                formula="feature", show_progress=False)
        s   = cph.summary
        hr  = float(s["exp(coef)"])
        lo  = float(s["exp(coef) lower 95%"])
        hi  = float(s["exp(coef) upper 95%"])
        pv  = float(s["p"])
        rows.append({"Factor": fac, "HR": hr, "CI_lo": lo,
                     "CI_hi": hi, "p": pv, "type": "Univariate"})

    # ── Multivariate Cox ───────────────────────────────────────────────────
    cov_cols = factor_cols + ["age_at_dx"]
    tmp = df[["os_months", "os_event"] + cov_cols].copy()
    cph_multi = CoxPHFitter()
    cph_multi.fit(tmp, duration_col="os_months", event_col="os_event",
                  show_progress=False)
    s = cph_multi.summary
    for fac in factor_cols:
        if fac in s.index:
            hr  = float(s.loc[fac, "exp(coef)"])
            lo  = float(s.loc[fac, "exp(coef) lower 95%"])
            hi  = float(s.loc[fac, "exp(coef) upper 95%"])
            pv  = float(s.loc[fac, "p"])
            rows.append({"Factor": fac, "HR": hr, "CI_lo": lo,
                         "CI_hi": hi, "p": pv, "type": "Multivariate"})

    cox_df = pd.DataFrame(rows)
    cox_df.to_csv(os.path.join(TABS, "cox_results.csv"), index=False)

    # ── Forest plot ────────────────────────────────────────────────────────
    factor_labels = {
        "Factor1": "Factor 1\n(Proneural/Mes axis)",
        "Factor2": "Factor 2\n(Mutation signature)",
        "Factor3": "Factor 3\n(G-CIMP/MGMT)",
        "Factor4": "Factor 4\n(Expression module)",
        "Factor5": "Factor 5\n(Minor signal)",
    }

    fig, axes = plt.subplots(1, 2, figsize=(14, 6),
                             gridspec_kw={"wspace": 0.45})

    for ax, typ, color in [(axes[0], "Univariate", "#4E79A7"),
                           (axes[1], "Multivariate", "#E15759")]:
        sub = cox_df[cox_df["type"] == typ].copy().reset_index(drop=True)
        if sub.empty:
            ax.set_visible(False)
            continue

        y_pos = list(range(len(sub)))[::-1]

        ax.axvline(1.0, ls="--", lw=1.0, c="#888888")

        for i, (_, row) in enumerate(sub.iterrows()):
            yp  = y_pos[i]
            lab = factor_labels.get(row["Factor"], row["Factor"])
            pv  = row["p"]
            sig = "***" if pv < 0.001 else "**" if pv < 0.01 else \
                  "*" if pv < 0.05 else ""

            c_dot = "#cc0000" if row["HR"] > 1 else "#0055aa"
            ax.errorbar(row["HR"], yp,
                        xerr=[[row["HR"] - row["CI_lo"]],
                               [row["CI_hi"] - row["HR"]]],
                        fmt="o", color=c_dot, ecolor=c_dot,
                        capsize=4, markersize=9, linewidth=1.8, zorder=3)

            ax.text(-0.02, yp, lab, ha="right", va="center",
                    fontsize=8.5, transform=ax.get_yaxis_transform())
            ax.text(1.02, yp,
                    f"HR={row['HR']:.2f} ({row['CI_lo']:.2f}–{row['CI_hi']:.2f})  "
                    f"p={pv:.2e} {sig}",
                    ha="left", va="center", fontsize=7.5,
                    transform=ax.get_yaxis_transform(),
                    color="#cc0000" if pv < 0.05 else "#555555")

        ax.set_xlim(-0.5, max(sub["CI_hi"].max() + 0.5, 3))
        ax.set_yticks([])
        ax.set_xlabel("Hazard Ratio (HR)  [95% CI]", fontsize=9.5)
        ax.set_title(f"{typ} Cox PH Model\n(HR > 1  →  worse survival)",
                     fontsize=10, fontweight="bold")

        # Shade regions
        xl = ax.get_xlim()
        ax.axvspan(xl[0], 1, alpha=0.04, color="#4E79A7", zorder=0)
        ax.axvspan(1, xl[1], alpha=0.04, color="#E15759", zorder=0)
        ax.text(0.25, 0.02, "← Protective", transform=ax.transAxes,
                fontsize=8, color="#4E79A7", style="italic")
        ax.text(0.62, 0.02, "Risk factor →", transform=ax.transAxes,
                fontsize=8, color="#E15759", style="italic")

    fig.suptitle("Cox Proportional-Hazards Regression: MOFA+ Factors and Overall Survival",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    fig.savefig(os.path.join(FIGS, "07_cox_forest.png"))
    plt.close(fig)
    return cox_df


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("  Survival Analysis  –  TCGA-GBM")
    print("=" * 60)

    df = load_data()
    print(f"  Patients: {len(df)} | Events: {df['os_event'].sum()}")
    print(f"  Subtype OS (median):")
    for st in ["PN","CL","MES"]:
        med = df.loc[df["subtype"]==st, "os_months"].median()
        print(f"    {st}: {med:.1f} months")

    lr_results = fig_km_curves(df)
    cox_df     = fig_cox_forest(df)

    print("\n  ── Log-rank test p-values ──────────────────")
    for k, v in lr_results.items():
        stars = "***" if v < 0.001 else "**" if v < 0.01 else "*" if v < 0.05 else "ns"
        print(f"    {k:<12}: p = {v:.2e}  {stars}")

    print("\n  ── Significant Cox factors (p < 0.05) ─────")
    sig = cox_df[cox_df["p"] < 0.05]
    for _, r in sig.iterrows():
        print(f"    {r['Factor']:<12} [{r['type']:<12}]  "
              f"HR={r['HR']:.3f}  p={r['p']:.2e}")

    print("\n  ✓  Survival analysis complete.")


if __name__ == "__main__":
    main()
