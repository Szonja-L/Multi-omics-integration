"""
05_fusion_comparison.py
═══════════════════════════════════════════════════════════════════════════════
Compare integration strategies for TCGA-GBM multi-omics survival prediction.

Strategies:
  • Single-view baselines  — Expression-only, Methylation-only, Mutation-only
  • Early fusion           — Concatenate all features → PCA(50) → CPH
  • MOFA+ integration      — Latent factor scores → CPH
  • Late fusion            — Per-view CPH scores → ensemble (mean) → C-index

Metric: Harrell's C-index (5-fold CV, mean ± SD)
Output: figures/08_fusion_comparison.png
        figures/09_modality_contribution.png
═══════════════════════════════════════════════════════════════════════════════
"""

import warnings, pathlib
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold

# ──────────────────────────────────────────────────────────────────────────────
#  Paths
# ──────────────────────────────────────────────────────────────────────────────
BASE   = pathlib.Path(__file__).resolve().parents[1]
DATA   = BASE / "data" / "simulated"
RES    = BASE / "results"
FIG    = RES / "figures"
TABLES = RES / "tables"
FIG.mkdir(parents=True, exist_ok=True)

# ──────────────────────────────────────────────────────────────────────────────
#  Palette (Tableau-10 consistent with earlier figures)
# ──────────────────────────────────────────────────────────────────────────────
COLORS = {
    "Expression":   "#59A14F",
    "Methylation":  "#B07AA1",
    "Mutations":    "#E15759",
    "Early Fusion": "#4E79A7",
    "MOFA+":        "#F28E2B",
    "Late Fusion":  "#76B7B2",
}
METHOD_ORDER = ["Expression", "Methylation", "Mutations",
                "Early Fusion", "MOFA+", "Late Fusion"]

print("=" * 60)
print("  Fusion Strategy Comparison  –  TCGA-GBM")
print("=" * 60)

# ──────────────────────────────────────────────────────────────────────────────
#  Load data
# ──────────────────────────────────────────────────────────────────────────────
expr   = pd.read_csv(DATA / "expression_matrix.csv",   index_col=0)
meth   = pd.read_csv(DATA / "methylation_matrix.csv",  index_col=0)
mut    = pd.read_csv(DATA / "mutation_matrix.csv",     index_col=0)
clin   = pd.read_csv(DATA / "clinical_data.csv",       index_col=0)
factors= pd.read_csv(TABLES / "factor_scores.csv",     index_col=0)

# Align indices
common = clin.index.intersection(expr.index).intersection(
         meth.index).intersection(mut.index).intersection(factors.index)
clin    = clin.loc[common]
expr    = expr.loc[common]
meth    = meth.loc[common]
mut     = mut.loc[common]
factors = factors.loc[common]

os_time  = clin["os_months"].values
os_event = clin["os_event"].values.astype(int)
subtypes = clin["subtype"].values

print(f"  Patients: {len(common)} | Events: {os_event.sum()}")

# ──────────────────────────────────────────────────────────────────────────────
#  Helper: C-index via 5-fold CV using CoxPH on reduced features
# ──────────────────────────────────────────────────────────────────────────────
N_SPLITS = 5
rng      = np.random.RandomState(42)

def cv_cindex(X: np.ndarray, T: np.ndarray, E: np.ndarray,
              n_splits: int = N_SPLITS) -> tuple[float, float]:
    """5-fold CV Harrell C-index for a CoxPH model fit on X."""
    scaler = StandardScaler()
    skf    = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    scores = []
    for train_idx, test_idx in skf.split(X, E):
        Xtr = scaler.fit_transform(X[train_idx])
        Xte = scaler.transform(X[test_idx])
        Ttr, Etr = T[train_idx], E[train_idx]
        Tte, Ete = T[test_idx],  E[test_idx]

        df_tr = pd.DataFrame(Xtr, columns=[f"f{i}" for i in range(Xtr.shape[1])])
        df_tr["T"] = Ttr; df_tr["E"] = Etr

        df_te = pd.DataFrame(Xte, columns=[f"f{i}" for i in range(Xte.shape[1])])

        cph = CoxPHFitter(penalizer=0.1)
        cph.fit(df_tr, duration_col="T", event_col="E")
        risk = cph.predict_log_partial_hazard(df_te).values
        ci   = concordance_index(Tte, -risk, Ete)
        scores.append(ci)
    return float(np.mean(scores)), float(np.std(scores))


def pca_reduce(X: np.ndarray, n_components: int = 20) -> np.ndarray:
    """PCA to n_components (or min(X.shape)-1)."""
    nc = min(n_components, X.shape[0]-1, X.shape[1])
    return PCA(n_components=nc, random_state=42).fit_transform(X)


# ──────────────────────────────────────────────────────────────────────────────
#  1.  Single-view baselines
# ──────────────────────────────────────────────────────────────────────────────
print("\n  ── Single-view baselines ─────────────────────────────")
results = {}

for name, mat in [("Expression", expr.values),
                  ("Methylation", meth.values),
                  ("Mutations",   mut.values)]:
    X_pca = pca_reduce(mat)
    mu, sd = cv_cindex(X_pca, os_time, os_event)
    results[name] = (mu, sd)
    print(f"    {name:<14} C-index = {mu:.3f} ± {sd:.3f}")

# ──────────────────────────────────────────────────────────────────────────────
#  2.  Early fusion
# ──────────────────────────────────────────────────────────────────────────────
print("\n  ── Early fusion (concat → PCA) ───────────────────────")
X_concat = np.hstack([expr.values, meth.values, mut.values.astype(float)])
X_early  = pca_reduce(X_concat, n_components=50)
mu, sd   = cv_cindex(X_early, os_time, os_event)
results["Early Fusion"] = (mu, sd)
print(f"    Early Fusion   C-index = {mu:.3f} ± {sd:.3f}")

# ──────────────────────────────────────────────────────────────────────────────
#  3.  MOFA+
# ──────────────────────────────────────────────────────────────────────────────
print("\n  ── MOFA+ (latent factors → CPH) ─────────────────────")
X_mofa  = factors.values
mu, sd  = cv_cindex(X_mofa, os_time, os_event)
results["MOFA+"] = (mu, sd)
print(f"    MOFA+          C-index = {mu:.3f} ± {sd:.3f}")

# ──────────────────────────────────────────────────────────────────────────────
#  4.  Late fusion  (per-view risk scores → ensemble)
# ──────────────────────────────────────────────────────────────────────────────
print("\n  ── Late fusion (per-view → ensemble) ────────────────")

def cv_risk_scores(X: np.ndarray, T: np.ndarray, E: np.ndarray,
                   n_splits: int = N_SPLITS):
    """Return out-of-fold risk score array and corresponding T, E."""
    scaler = StandardScaler()
    skf    = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    oof_risk = np.zeros(len(X))
    oof_T    = np.zeros(len(X))
    oof_E    = np.zeros(len(X), dtype=int)
    for train_idx, test_idx in skf.split(X, E):
        Xtr = scaler.fit_transform(X[train_idx])
        Xte = scaler.transform(X[test_idx])
        Ttr, Etr = T[train_idx], E[train_idx]
        df_tr = pd.DataFrame(Xtr, columns=[f"f{i}" for i in range(Xtr.shape[1])])
        df_tr["T"] = Ttr; df_tr["E"] = Etr
        df_te = pd.DataFrame(Xte, columns=[f"f{i}" for i in range(Xte.shape[1])])
        cph = CoxPHFitter(penalizer=0.1)
        cph.fit(df_tr, duration_col="T", event_col="E")
        oof_risk[test_idx] = cph.predict_log_partial_hazard(df_te).values
        oof_T[test_idx]    = T[test_idx]
        oof_E[test_idx]    = E[test_idx]
    return oof_risk, oof_T, oof_E

late_scores_list = []
for mat in [pca_reduce(expr.values), pca_reduce(meth.values),
            pca_reduce(mut.values.astype(float))]:
    r, _, _ = cv_risk_scores(mat, os_time, os_event)
    late_scores_list.append(r)

ensemble_risk = np.mean(late_scores_list, axis=0)
ci_late = concordance_index(os_time, -ensemble_risk, os_event)
results["Late Fusion"] = (ci_late, 0.0)   # single pass, no SD
print(f"    Late Fusion    C-index = {ci_late:.3f}  (ensemble OOF)")

# ──────────────────────────────────────────────────────────────────────────────
#  Save C-index table
# ──────────────────────────────────────────────────────────────────────────────
df_res = pd.DataFrame(
    {m: {"c_index_mean": v[0], "c_index_std": v[1]} for m, v in results.items()}
).T
df_res.index.name = "method"
df_res.to_csv(TABLES / "fusion_comparison.csv")
print(f"\n  Saved fusion_comparison.csv")

# ──────────────────────────────────────────────────────────────────────────────
#  Figure 8 — Fusion comparison bar chart
# ──────────────────────────────────────────────────────────────────────────────
print("\n  → Figure 8: Fusion strategy comparison…")

fig, ax = plt.subplots(figsize=(9, 5.5))

means = [results[m][0] for m in METHOD_ORDER]
stds  = [results[m][1] for m in METHOD_ORDER]
colors= [COLORS[m]     for m in METHOD_ORDER]
x     = np.arange(len(METHOD_ORDER))

bars = ax.bar(x, means, width=0.62, color=colors, alpha=0.88,
              linewidth=0, zorder=3)
ax.errorbar(x, means, yerr=stds, fmt="none", color="#333333",
            capsize=5, capthick=1.5, linewidth=1.5, zorder=4)

# Annotate values
for xi, (mu, sd) in zip(x, zip(means, stds)):
    label = f"{mu:.3f}"
    ax.text(xi, mu + sd + 0.008, label, ha="center", va="bottom",
            fontsize=10, fontweight="bold", color="#222222")

# Reference line at 0.5
ax.axhline(0.5, color="#888888", linewidth=1.2, linestyle="--",
           label="Random chance (C=0.5)", zorder=2)

# Category shading
ax.axvspan(-0.5, 2.5, color="#F0F0F0", alpha=0.55, zorder=1)
ax.axvspan(2.5,  5.5, color="#E8F4FD", alpha=0.55, zorder=1)
ax.text(1.0,  0.525, "Single-view\nbaselines", ha="center", va="bottom",
        fontsize=9, color="#666666", style="italic")
ax.text(4.0,  0.525, "Integration\nstrategies", ha="center", va="bottom",
        fontsize=9, color="#666666", style="italic")

ax.set_xticks(x)
ax.set_xticklabels(METHOD_ORDER, fontsize=11)
ax.set_ylabel("Harrell's C-index (5-fold CV)", fontsize=12)
ax.set_title("Multi-Omics Integration Strategy Comparison\nTCGA-GBM Survival Prediction",
             fontsize=13, fontweight="bold", pad=12)
ax.set_ylim(0.48, max(means) + max(stds) + 0.08)
ax.yaxis.grid(True, alpha=0.4, zorder=0)
ax.set_axisbelow(True)
ax.spines[["top","right"]].set_visible(False)
ax.legend(fontsize=9, loc="upper left")

# Highlight MOFA+ bar
mofa_idx = METHOD_ORDER.index("MOFA+")
bars[mofa_idx].set_edgecolor("#222222")
bars[mofa_idx].set_linewidth(2.0)

plt.tight_layout()
fig.savefig(FIG / "08_fusion_comparison.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print(f"    Saved 08_fusion_comparison.png")

# ──────────────────────────────────────────────────────────────────────────────
#  Figure 9 — Modality contribution (R² per factor × view)
# ──────────────────────────────────────────────────────────────────────────────
print("\n  → Figure 9: Modality contribution per factor…")

var_df = pd.read_csv(TABLES / "variance_explained.csv", index_col=0)
# rows = factors, cols = views
var_df = var_df.loc[var_df.index.notna()]

# Clip to active factors only (non-zero rows)
row_sums = var_df.sum(axis=1)
active   = var_df.loc[row_sums > 0].copy()

n_factors = len(active)
n_views   = len(active.columns)
bar_width  = 0.22
x2         = np.arange(n_factors)

view_colors = {"expression": "#59A14F", "methylation": "#B07AA1", "mutations": "#E15759"}

fig2, ax2 = plt.subplots(figsize=(10, 5.5))

offsets = np.linspace(-(n_views-1)/2, (n_views-1)/2, n_views) * bar_width

for i, view in enumerate(active.columns):
    vals   = active[view].values
    color  = view_colors.get(view.lower(), "#999999")
    label  = view.capitalize()
    ax2.bar(x2 + offsets[i], vals, bar_width, label=label,
            color=color, alpha=0.85, zorder=3)

ax2.set_xticks(x2)
ax2.set_xticklabels(active.index, fontsize=11)
ax2.set_ylabel("Variance Explained (%)", fontsize=12)
ax2.set_title("MOFA+ Modality Contribution per Latent Factor\nTCGA-GBM",
              fontsize=13, fontweight="bold", pad=12)
ax2.yaxis.grid(True, alpha=0.4, zorder=0)
ax2.set_axisbelow(True)
ax2.spines[["top","right"]].set_visible(False)
ax2.legend(title="Modality", fontsize=10, title_fontsize=10, loc="upper right")

plt.tight_layout()
fig2.savefig(FIG / "09_modality_contribution.png", dpi=160, bbox_inches="tight")
plt.close(fig2)
print(f"    Saved 09_modality_contribution.png")

print("\n  ✓  Fusion comparison complete.")
print("=" * 60)
