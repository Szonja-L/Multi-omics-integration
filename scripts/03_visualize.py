#!/usr/bin/env python3
"""
03_visualize.py
================
Generate all publication-quality figures for the Multi-Omics Integration project.

Figures produced:
  01_data_overview.png          – Heatmaps of all three modalities
  02_variance_explained.png     – MOFA+ variance decomposition
  03_factor_umap.png            – UMAP of factor space
  04_factor_distributions.png  – Factor scores by subtype (violin)
  05_top_weights.png            – Top feature weights per factor/view
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, Normalize
import seaborn as sns
import umap
import warnings, os

warnings.filterwarnings("ignore")
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
os.makedirs(FIGS, exist_ok=True)

# ── Palette ───────────────────────────────────────────────────────────────────
SUB_COL = {
    "PN":  "#4E79A7",
    "CL":  "#F28E2B",
    "MES": "#E15759",
}
SUB_FULL = {"PN": "Proneural", "CL": "Classical", "MES": "Mesenchymal"}
VIEW_COL = {
    "Expression":   "#59A14F",
    "Methylation":  "#B07AA1",
    "Mutations":    "#E15759",
}


# ═══════════════════════════════════════════════════════════════════════════════
#  HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def load_all():
    expr = pd.read_csv(os.path.join(DATA, "expression_matrix.csv"),  index_col=0)
    meth = pd.read_csv(os.path.join(DATA, "methylation_matrix.csv"), index_col=0)
    muts = pd.read_csv(os.path.join(DATA, "mutation_matrix.csv"),    index_col=0)
    clin = pd.read_csv(os.path.join(DATA, "clinical_data.csv"))
    Z    = pd.read_csv(os.path.join(TABS, "factor_scores.csv"),      index_col=0)
    r2   = pd.read_csv(os.path.join(TABS, "variance_explained.csv"), index_col=0)
    r2t  = pd.read_csv(os.path.join(TABS, "variance_explained_total.csv"),
                       index_col=0, header=None).squeeze()
    r2t  = r2t[r2t.index.notna()]   # drop NaN-indexed header artefact row
    W_e  = pd.read_csv(os.path.join(TABS, "weights_expression.csv"), index_col=0)
    W_m  = pd.read_csv(os.path.join(TABS, "weights_methylation.csv"),index_col=0)
    W_u  = pd.read_csv(os.path.join(TABS, "weights_mutations.csv"),  index_col=0)
    return expr, meth, muts, clin, Z, r2, r2t, W_e, W_m, W_u


def subtype_colors(clin: pd.DataFrame) -> list:
    return [SUB_COL[s] for s in clin["subtype"].tolist()]


def figpath(name: str) -> str:
    return os.path.join(FIGS, name)


# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE 1 – DATA OVERVIEW
# ═══════════════════════════════════════════════════════════════════════════════

def fig_data_overview(expr, meth, muts, clin):
    """Three annotated heatmaps of the input data matrices."""
    print("  → Figure 1: Data overview heatmaps…")

    subtypes = clin["subtype"].values
    order    = np.argsort(subtypes)           # sort samples by subtype

    fig = plt.figure(figsize=(18, 6))
    gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.35)

    def _draw_heatmap(ax, mat, title, cmap, vmin, vmax, xlabel):
        n_show = min(80, mat.shape[1])
        idx    = np.random.choice(mat.shape[1], n_show, replace=False)
        data   = mat.iloc[order, :].iloc[:, idx]

        im = ax.imshow(data.values, aspect="auto", cmap=cmap,
                       vmin=vmin, vmax=vmax, interpolation="nearest")
        plt.colorbar(im, ax=ax, shrink=0.65, pad=0.02)

        # Subtype bars on the left
        tick_pos, tick_lab = [], []
        for st in ["PN", "CL", "MES"]:
            idx_st = np.where(subtypes[order] == st)[0]
            ax.axhspan(idx_st.min()-0.5, idx_st.max()+0.5,
                       xmin=0, xmax=0.04, color=SUB_COL[st], clip_on=False, lw=0)
            tick_pos.append(idx_st.mean())
            tick_lab.append(SUB_FULL[st])

        ax.set_yticks(tick_pos)
        ax.set_yticklabels(tick_lab, fontsize=9)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_xticks([])
        ax.set_title(title, fontweight="bold", fontsize=11, pad=8)

    _draw_heatmap(fig.add_subplot(gs[0]),
                  expr, "Gene Expression\n(log₂ TPM)", "RdYlBu_r", 3, 12,
                  f"{expr.shape[1]} genes  (80 shown)")
    _draw_heatmap(fig.add_subplot(gs[1]),
                  meth, "DNA Methylation\n(β values)", "RdBu_r", 0.05, 0.95,
                  f"{meth.shape[1]} genes (HM450, gene-level)")
    _draw_heatmap(fig.add_subplot(gs[2]),
                  muts, "Somatic Mutations\n(binary)", "Greys", 0, 1,
                  f"{muts.shape[1]} driver genes  (all shown)")

    # Legend
    patches = [mpatches.Patch(color=SUB_COL[s], label=SUB_FULL[s])
               for s in ["PN","CL","MES"]]
    fig.legend(handles=patches, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 1.03), frameon=False, fontsize=9)

    fig.suptitle(f"TCGA-GBM Multi-Omics Data Overview  (N={expr.shape[0]})", y=1.07,
                 fontsize=13, fontweight="bold")
    fig.savefig(figpath("01_data_overview.png"))
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE 2 – VARIANCE EXPLAINED
# ═══════════════════════════════════════════════════════════════════════════════

def fig_variance_explained(r2: pd.DataFrame, r2t: pd.Series):
    """Stacked bar of R² per factor + total R² bar."""
    print("  → Figure 2: Variance explained…")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5),
                             gridspec_kw={"width_ratios": [3, 1]})

    # ── Left: stacked bar per factor ──────────────────────────────────────
    ax = axes[0]
    views  = ["Expression", "Methylation", "Mutations"]
    colors = [VIEW_COL[v] for v in views]
    bottom = np.zeros(len(r2))

    for v, c in zip(views, colors):
        vals = r2[v].values
        bars = ax.bar(r2.index, vals, bottom=bottom, color=c,
                      label=v, edgecolor="white", linewidth=0.5)
        # Add value labels for large bars
        for bar, val, bot in zip(bars, vals, bottom):
            if val > 2.5:
                ax.text(bar.get_x() + bar.get_width()/2,
                        bot + val/2, f"{val:.0f}%",
                        ha="center", va="center", fontsize=8.5,
                        fontweight="bold", color="white")
        bottom += vals

    ax.set_xlabel("MOFA+ Factor", fontsize=10)
    ax.set_ylabel("Marginal variance explained (%)", fontsize=10)
    ax.set_title("Variance Decomposition by Factor and View",
                 fontsize=11, fontweight="bold")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.set_ylim(0, bottom.max() * 1.18)
    ax.tick_params(axis="x", labelsize=9)

    # Annotation: factor labels
    factor_annot = {
        "Factor1": "Mutation\nlandscape",
        "Factor2": "Transcriptional\nsubtype",
        "Factor3": "Expression\nmodule A",
        "Factor4": "Expression\nmodule B",
        "Factor5": "Methylation\naxis",
    }
    for i, (fac, lab) in enumerate(factor_annot.items()):
        if fac in r2.index:
            ax.text(i, bottom[i] + 1.5, lab, ha="center",
                    va="bottom", fontsize=7.5, color="#444444",
                    style="italic")

    # ── Right: total R² per view ──────────────────────────────────────────
    ax2   = axes[1]
    vnames = list(r2t.index)
    vals   = r2t.values
    colors2 = [VIEW_COL[v] for v in vnames]

    bars = ax2.barh(vnames, vals, color=colors2, edgecolor="white")
    for bar, val in zip(bars, vals):
        ax2.text(val + 0.8, bar.get_y() + bar.get_height()/2,
                 f"{val:.1f}%", va="center", fontsize=9, fontweight="bold")

    ax2.set_xlim(0, 115)
    ax2.set_xlabel("Total R² explained (%)", fontsize=10)
    ax2.set_title("Total per View", fontsize=11, fontweight="bold")
    ax2.invert_yaxis()

    plt.tight_layout()
    fig.savefig(figpath("02_variance_explained.png"))
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE 3 – UMAP OF FACTOR SPACE
# ═══════════════════════════════════════════════════════════════════════════════

def fig_factor_umap(Z: pd.DataFrame, clin: pd.DataFrame):
    """UMAP of MOFA+ factor scores, coloured by subtype and Factor1."""
    print("  → Figure 3: UMAP of factor space…")

    reducer = umap.UMAP(n_neighbors=15, min_dist=0.2, random_state=42)
    emb     = reducer.fit_transform(Z.values)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # ── Left: by subtype ──────────────────────────────────────────────────
    ax = axes[0]
    for st in ["PN", "CL", "MES"]:
        mask = clin["subtype"].values == st
        ax.scatter(emb[mask, 0], emb[mask, 1],
                   c=SUB_COL[st], label=SUB_FULL[st],
                   s=45, alpha=0.80, edgecolors="white", linewidth=0.4)

    ax.set_xlabel("UMAP 1", fontsize=10)
    ax.set_ylabel("UMAP 2", fontsize=10)
    ax.set_title("MOFA+ Factor Space – Molecular Subtype",
                 fontsize=11, fontweight="bold")
    ax.legend(frameon=True, framealpha=0.9, fontsize=9,
              edgecolor="#cccccc")
    ax.set_xticks([]); ax.set_yticks([])

    # ── Right: by Factor1 score (continuous) ─────────────────────────────
    ax2 = axes[1]
    f1  = Z["Factor1"].values
    sc  = ax2.scatter(emb[:, 0], emb[:, 1], c=f1,
                      cmap="RdBu_r", s=45, alpha=0.85,
                      edgecolors="white", linewidth=0.4,
                      vmin=np.percentile(f1, 5), vmax=np.percentile(f1, 95))
    plt.colorbar(sc, ax=ax2, shrink=0.75, label="Factor 1 score")
    ax2.set_xlabel("UMAP 1", fontsize=10)
    ax2.set_ylabel("UMAP 2", fontsize=10)
    ax2.set_title("MOFA+ Factor Space – Factor 1 (Mutation landscape)",
                  fontsize=11, fontweight="bold")
    ax2.set_xticks([]); ax2.set_yticks([])

    plt.tight_layout()
    fig.savefig(figpath("03_factor_umap.png"))
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE 4 – FACTOR DISTRIBUTIONS BY SUBTYPE
# ═══════════════════════════════════════════════════════════════════════════════

def fig_factor_distributions(Z: pd.DataFrame, clin: pd.DataFrame):
    """Violin plots of factor scores stratified by molecular subtype."""
    print("  → Figure 4: Factor distributions by subtype…")

    K     = min(5, len(Z.columns))
    fig, axes = plt.subplots(1, K, figsize=(3.5*K, 5), sharey=False)

    factor_titles = [
        "Factor 1\nMutation landscape",
        "Factor 2\nTranscriptional\nsubtype",
        "Factor 3\nExpression\nmodule A",
        "Factor 4\nExpression\nmodule B",
        "Factor 5\nMethylation\naxis",
    ]

    for i, (fac, ax) in enumerate(zip(Z.columns[:K], axes)):
        subtypes_arr = clin["subtype"].values
        data_by_st = [Z.loc[Z.index[subtypes_arr == st], fac].values
                      for st in ["PN", "CL", "MES"]]
        labels     = [SUB_FULL[s] for s in ["PN", "CL", "MES"]]

        vp = ax.violinplot(data_by_st, positions=[0, 1, 2],
                           showmedians=True, showextrema=False)
        for j, (body, st) in enumerate(zip(vp["bodies"], ["PN","CL","MES"])):
            body.set_facecolor(SUB_COL[st])
            body.set_alpha(0.75)
            body.set_edgecolor(SUB_COL[st])

        vp["cmedians"].set_color("#333333")
        vp["cmedians"].set_linewidth(2)

        # Strip points
        for j, st in enumerate(["PN","CL","MES"]):
            mask = subtypes_arr == st
            y = Z.loc[Z.index[mask], fac].values
            x = np.random.normal(j, 0.04, len(y))
            ax.scatter(x, y, c=SUB_COL[st], s=8, alpha=0.55,
                       edgecolors="white", linewidth=0.2, zorder=3)

        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(labels, fontsize=8.5, rotation=30, ha="right")
        ax.set_title(factor_titles[i] if i < len(factor_titles) else fac,
                     fontsize=9.5, fontweight="bold")
        ax.axhline(0, ls="--", lw=0.8, c="#aaaaaa")
        if i == 0:
            ax.set_ylabel("Factor score", fontsize=10)

    fig.suptitle("MOFA+ Factor Scores by GBM Molecular Subtype",
                 fontsize=12, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(figpath("04_factor_distributions.png"))
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE 5 – TOP FEATURE WEIGHTS
# ═══════════════════════════════════════════════════════════════════════════════

def fig_top_weights(W_e, W_m, W_u):
    """Top loadings for Factor1 and Factor3 across all three views."""
    print("  → Figure 5: Top feature weights…")

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    configs = [
        # (factor, W, view_name, color, top_n, row, col)
        ("Factor1", W_e, "Expression",   VIEW_COL["Expression"],  16, 0, 0),
        ("Factor1", W_m, "Methylation",  VIEW_COL["Methylation"], 16, 0, 1),
        ("Factor1", W_u, "Mutations",    VIEW_COL["Mutations"],   15, 0, 2),
        ("Factor2", W_e, "Expression",   VIEW_COL["Expression"],  16, 1, 0),
        ("Factor2", W_m, "Methylation",  VIEW_COL["Methylation"], 16, 1, 1),
        ("Factor2", W_u, "Mutations",    VIEW_COL["Mutations"],   15, 1, 2),
    ]

    factor_labels = {
        "Factor1": "Factor 1  (Mutation landscape)",
        "Factor2": "Factor 2  (Transcriptional subtype axis)",
    }

    for (fac, W, view, col, top_n, row, c) in configs:
        ax   = axes[row][c]
        if fac not in W.columns:
            ax.set_visible(False)
            continue

        weights = W[fac].sort_values()
        n_top   = min(top_n // 2, len(weights))
        top     = pd.concat([weights.head(n_top), weights.tail(n_top)])
        top     = top.sort_values()

        colors = [col if w > 0 else "#aaaaaa" for w in top.values]
        bars   = ax.barh(range(len(top)), top.values, color=colors,
                         edgecolor="white", linewidth=0.4)

        # Clean gene/CpG labels
        labels = [g.replace("_", " ").replace("GENE ", "").replace("cg ", "")
                  for g in top.index]
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(labels, fontsize=7.5)
        ax.axvline(0, color="#666666", lw=0.8)
        ax.set_title(f"{view}", fontsize=10, fontweight="bold",
                     color=col, pad=4)
        ax.set_xlabel("Factor weight", fontsize=8.5)
        ax.tick_params(axis="x", labelsize=8)

    # Row super-titles
    for row_idx, fac in enumerate(["Factor1", "Factor2"]):
        fig.text(0.5, 0.97 - row_idx * 0.50,
                 factor_labels[fac], ha="center", fontsize=11,
                 fontweight="bold", va="top",
                 bbox=dict(facecolor="#f5f5f5", edgecolor="#dddddd",
                           boxstyle="round,pad=0.3"))

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(figpath("05_top_weights.png"))
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("  MOFA+ Visualisation  –  Figures 1–5")
    print("=" * 60)

    (expr, meth, muts, clin,
     Z, r2, r2t, W_e, W_m, W_u) = load_all()

    fig_data_overview(expr, meth, muts, clin)
    fig_variance_explained(r2, r2t)
    fig_factor_umap(Z, clin)
    fig_factor_distributions(Z, clin)
    fig_top_weights(W_e, W_m, W_u)

    print(f"\n  ✓  All figures saved to results/figures/")


if __name__ == "__main__":
    main()
