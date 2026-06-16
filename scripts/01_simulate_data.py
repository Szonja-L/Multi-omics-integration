#!/usr/bin/env python3
"""
01_simulate_data.py
====================
Simulate biologically realistic TCGA-GBM multi-omics data.

Generates three paired data matrices + clinical data:
  - Gene expression (RNA-seq log2-normalised, 400 genes × 150 patients)
  - DNA methylation (450K beta values, 300 CpG sites × 150 patients)
  - Somatic mutations (WES binary, 50 genes × 150 patients)

Biological grounding
---------------------
GBM is stratified into three molecular subtypes (Verhaak et al. 2010):
  • Proneural (PN)  – IDH-mutant, G-CIMP hypermethylation, MGMT methylated
  • Classical (CL)  – EGFR-amplified, CDKN2A deletion, TERT mutant
  • Mesenchymal (MES) – NF1/PTEN mutant, CHI3L1 high, worst prognosis

The simulation is driven by 5 biological latent factors so MOFA+ can
recover them meaningfully during integration.
"""

import numpy as np
import pandas as pd
from scipy.stats import weibull_min
import os, warnings

warnings.filterwarnings("ignore")

# ── Reproducibility ──────────────────────────────────────────────────────────
SEED = 42
rng  = np.random.default_rng(SEED)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data", "simulated")
os.makedirs(DATA, exist_ok=True)

# ── Dimensions ────────────────────────────────────────────────────────────────
N_PER   = 50          # patients per subtype
N       = N_PER * 3   # 150 total
N_GENES = 400         # expression features
N_CPGS  = 300         # methylation features
N_MUTS  = 50          # mutation features
K_TRUE  = 5           # biological latent factors


# ═══════════════════════════════════════════════════════════════════════════════
#  1.  GENE / CpG / MUTATION NAMES
# ═══════════════════════════════════════════════════════════════════════════════

def make_gene_names() -> list:
    """Curated GBM-relevant gene names (real + padded noise)."""
    proneural = [
        "IDH1","IDH2","PDGFRA","OLIG2","OLIG1","SOX2","SOX11","TCF7L2",
        "DLL3","NOTCH1","SEMA3C","NXPH1","CADPS","DSCAML1","PRICKLE1",
        "PTPRD","NRXN1","NR2E1","GPR17","PLCL1","CAMK2A","SHANK2","SYT1",
        "GABRA1","GRIN1","ASCL1","DNMT3A","LDHB","ATP1A2","RAB3A",
        "NR4A1","DPYSL3","SLC17A6","STMN2","GAD1","GAD2","MAP2","TUBB3",
        "SCN2A","SNAP25",
    ]
    classical = [
        "EGFR","CDK4","MDM2","MDM4","CCND1","CDK6","E2F3","GLI1","SHH",
        "PTCH1","VEGFA","FGF2","FGFR1","IGFBP2","NES","FBLN1","COL4A1",
        "LAMC1","SPARC","HIF1A","CA9","CXCL12","PDGFA","TGFB1","TGFB2",
        "MMP9","MMP2","TIMP1","TWIST1","SNAI1","CDK6","E2F2","TGFB3",
        "ITGB3","THBS2","TNFRSF12A","LTBP1","LTBP2","FSCN1","CXCL1",
    ]
    mesenchymal = [
        "CHI3L1","MET","CD44","FN1","SERPINE1","LGALS3","TGFB1","ACTN1",
        "VIM","TNFRSF1A","TRADD","RELB","CASP1","CASP3","BCL2A1","TNC",
        "COL5A2","CXCR4","PLAUR","TRAF3","CD163","BMP2","SLC2A3","LOXL2",
        "NRP1","CTGF","CYR61","ANXA2","SPP1","ITGA3","IQGAP1","MYH9",
        "ARHGEF2","ROCK1","SNAI2","ZEB1","ZEB2","COL1A1","COL3A1","MXRA8",
    ]
    proliferation = [
        "MKI67","TOP2A","PCNA","MCM2","MCM6","CCNB1","CCNB2","CDK1","AURKB",
        "BUB1","CCNE1","CDC20","CDC25C","H2AFX","INCENP","PTTG1","SMC4",
        "TPX2","BIRC5","PLK1","CENPF","KIF11","KIF2C","MELK","TACC3",
        "UHRF1","CKS2","FOXM1","E2F1","RRM2",
    ]
    suppressor_signal = [
        "PTEN","RB1","ATRX","TP53","NF1","CDKN2A","CDKN2B","TSC1","TSC2",
        "FBXW7","AKT1","PIK3CA","PIK3R1","MTOR","PARP1","ATR","CHK1",
        "CHK2","RAD51","BRCA1","BRCA2","STAG2","MGMT","MLH1","MSH2",
        "ALDH1A1","SOX9","GFAP","S100B","NESTIN",
    ]
    noise_genes = [f"GENE_{i:03d}" for i in range(160)]
    return proneural + classical + mesenchymal + proliferation + suppressor_signal + noise_genes


def make_cpg_names() -> list:
    """CpG site identifiers (GCMP + MGMT + background + noise)."""
    gcimp = [f"cg_GCIMP_{i:04d}" for i in range(40)]   # G-CIMP panel
    mgmt  = [f"cg_MGMT_{i:02d}"  for i in range(20)]   # MGMT promoter
    prom  = [f"cg_PROM_{i:03d}"  for i in range(80)]   # cancer-related promoters
    noise = [f"cg_BACK_{i:03d}"  for i in range(160)]  # background
    return gcimp + mgmt + prom + noise


def make_mut_names() -> list:
    """Driver gene panel (real GBM recurrently mutated genes)."""
    drivers = [
        "IDH1","IDH2","EGFR","NF1","PTEN","TP53","ATRX","RB1","PIK3CA",
        "PIK3R1","TERT","CDK4","CDKN2A","CDKN2B","PDGFRA","STAG2","FUBP1",
        "CIC","NOTCH1","PTCH1","GLI1","AKT1","AKT2","CCND2","RPS6KA3",
        "BRAF","KRAS","NRAS","HRAS","MET","ERBB2","ERBB3","ERBB4","FGFR1",
        "FGFR3","ALK","RET","MYC","MYCN","MDM2","MDM4","BCL2","MCL1",
        "VHL","SMARCA4","ARID1A","KMT2D","KMT2C","SETD2","EP300",
    ]
    return drivers


# ═══════════════════════════════════════════════════════════════════════════════
#  2.  LATENT FACTOR GENERATION
# ═══════════════════════════════════════════════════════════════════════════════

def make_latent_factors(N: int, rng: np.random.Generator) -> np.ndarray:
    """
    Build 5 biological latent factors (N × 5):
      F1 – Proneural / IDH-mutant axis
      F2 – Mesenchymal / immune-rich axis
      F3 – Proliferation (continuous, grade-related)
      F4 – MGMT/G-CIMP methylation
      F5 – EGFR / Classical axis
    """
    idx = {
        "PN":  np.arange(0, N_PER),
        "CL":  np.arange(N_PER, 2*N_PER),
        "MES": np.arange(2*N_PER, 3*N_PER),
    }

    def subtype_factor(high, mid, low, noise=0.45):
        f = np.zeros(N)
        for st, val in [("PN", high), ("CL", mid), ("MES", low)]:
            f[idx[st]] = val + rng.normal(0, noise, N_PER)
        return f

    F1 = subtype_factor( 2.0,  0.3, -1.5)   # Proneural
    F2 = subtype_factor(-1.5,  0.4,  2.0)   # Mesenchymal
    F3 = rng.normal(0, 1.0, N)              # Proliferation (continuous)
    # Boost F3 slightly in MES (GBM III→IV progression)
    F3[idx["MES"]] += 0.5
    F4 = subtype_factor( 2.0,  0.0, -1.8)   # MGMT/G-CIMP
    F5 = subtype_factor(-0.8,  2.0, -0.3)   # EGFR/Classical

    return np.column_stack([F1, F2, F3, F4, F5])   # (N, 5)


# ═══════════════════════════════════════════════════════════════════════════════
#  3.  DATA MODALITY SIMULATION
# ═══════════════════════════════════════════════════════════════════════════════

def simulate_expression(Z: np.ndarray, rng: np.random.Generator, genes: list) -> pd.DataFrame:
    """
    RNA-seq log2(TPM+1) expression matrix.
    Genes 0–39  : Proneural markers  (heavy loading on F1)
    Genes 40–79 : Classical markers  (heavy loading on F5)
    Genes 80–119: Mesenchymal markers(heavy loading on F2)
    Genes 120–149: Proliferation     (heavy loading on F3)
    Genes 150–179: Suppressor/signal (light mixed loadings)
    Genes 180+  : Noise
    """
    D = len(genes)
    W = np.zeros((D, 5))

    # Proneural marker genes
    W[:40, 0]  =  np.abs(rng.normal(1.8, 0.4, 40))   # F1 high
    W[:40, 1]  = -np.abs(rng.normal(0.6, 0.2, 40))   # F2 low
    # Classical marker genes
    W[40:80, 4] =  np.abs(rng.normal(1.8, 0.4, 40))  # F5 high
    W[40:80, 1] = -np.abs(rng.normal(0.4, 0.2, 40))  # F2 low
    # Mesenchymal marker genes
    W[80:120, 1] =  np.abs(rng.normal(1.8, 0.4, 40)) # F2 high
    W[80:120, 0] = -np.abs(rng.normal(0.5, 0.2, 40)) # F1 low
    # Proliferation genes
    W[120:150, 2] =  np.abs(rng.normal(1.6, 0.4, 30))# F3 high
    # Suppressor/signal genes
    W[150:180, 0] =  rng.normal(0, 0.6, 30)           # mixed
    W[150:180, 1] =  rng.normal(0, 0.5, 30)

    # Observed expression = signal + noise + baseline
    signal  = Z @ W.T                              # (N, D)
    noise   = rng.normal(0, 0.8, (Z.shape[0], D))  # biological noise
    expr    = signal + noise + 6.5                 # log2 scale baseline
    expr    = np.clip(expr, 0.0, 16.0)             # realistic range

    return pd.DataFrame(expr, columns=genes,
                        index=[f"TCGA-{i:03d}" for i in range(Z.shape[0])])


def simulate_methylation(Z: np.ndarray, rng: np.random.Generator, cpgs: list) -> pd.DataFrame:
    """
    DNA methylation 450K beta values in [0, 1].
    CpGs 0–39  : G-CIMP panel (F4-driven, high in PN)
    CpGs 40–59 : MGMT promoter (F4-driven)
    CpGs 60–139: Mixed cancer-related CpGs (F1, F2 loadings)
    CpGs 140+  : Background noise
    """
    N, D = Z.shape[0], len(cpgs)
    W = np.zeros((D, 5))

    # G-CIMP hypermethylation (IDH-mutant)
    W[:40,  3] =  np.abs(rng.normal(0.9, 0.15, 40))  # F4 high → high beta
    # MGMT promoter methylation
    W[40:60, 3] =  np.abs(rng.normal(0.8, 0.12, 20)) # F4 high → methylated
    W[40:60, 0] =  rng.normal(0.2, 0.1, 20)
    # Mixed promoter CpGs
    W[60:100, 0] = rng.normal(0.0, 0.4, 40)           # some F1 driven
    W[60:100, 1] = rng.normal(0.0, 0.3, 40)           # some F2 driven
    W[100:140, 4] = rng.normal(0.0, 0.3, 40)          # some F5 driven

    # Raw signal → logistic → beta in [0,1]
    signal = Z @ W.T
    noise  = rng.normal(0, 0.3, (N, D))
    raw    = signal + noise
    beta   = 1.0 / (1.0 + np.exp(-raw))   # sigmoid mapping
    beta   = np.clip(beta, 0.01, 0.99)

    return pd.DataFrame(beta, columns=cpgs,
                        index=[f"TCGA-{i:03d}" for i in range(N)])


def simulate_mutations(Z: np.ndarray, rng: np.random.Generator, mut_genes: list,
                       subtype_labels: np.ndarray) -> pd.DataFrame:
    """
    Binary somatic mutation matrix.
    Mutations are drawn from subtype-specific Bernoulli probabilities
    grounded in TCGA-GBM mutation frequencies.
    """
    N   = Z.shape[0]
    D   = len(mut_genes)

    # Base mutation probability per gene per subtype
    # Rows: genes, Columns: [PN, CL, MES]
    probs = {
        "IDH1":   [0.82, 0.03, 0.04],
        "IDH2":   [0.06, 0.01, 0.01],
        "EGFR":   [0.05, 0.88, 0.22],
        "NF1":    [0.04, 0.06, 0.42],
        "PTEN":   [0.12, 0.38, 0.54],
        "TP53":   [0.68, 0.12, 0.62],
        "ATRX":   [0.72, 0.02, 0.08],
        "RB1":    [0.08, 0.20, 0.36],
        "PIK3CA": [0.08, 0.12, 0.14],
        "PIK3R1": [0.06, 0.10, 0.12],
        "TERT":   [0.28, 0.72, 0.42],
        "CDK4":   [0.04, 0.38, 0.12],
        "CDKN2A": [0.06, 0.52, 0.28],
        "CDKN2B": [0.06, 0.44, 0.24],
        "PDGFRA": [0.32, 0.08, 0.06],
        "STAG2":  [0.06, 0.06, 0.18],
        "FUBP1":  [0.14, 0.02, 0.04],
        "CIC":    [0.10, 0.02, 0.04],
        "NOTCH1": [0.08, 0.12, 0.14],
        "PTCH1":  [0.04, 0.12, 0.06],
    }

    mut = np.zeros((N, D))
    st_idx = {"PN": 0, "CL": 1, "MES": 2}

    for i, gene in enumerate(mut_genes):
        if gene in probs:
            p = probs[gene]
        else:
            # Background noise mutations
            p = [0.04, 0.06, 0.06]

        for j, st in enumerate(["PN", "CL", "MES"]):
            mask = subtype_labels == st
            mut[mask, i] = rng.binomial(1, p[j], mask.sum())

    return pd.DataFrame(mut.astype(int), columns=mut_genes,
                        index=[f"TCGA-{i:03d}" for i in range(N)])


def simulate_clinical(Z: np.ndarray, subtype_labels: np.ndarray,
                      rng: np.random.Generator) -> pd.DataFrame:
    """
    Simulate clinical data with survival times driven by latent factors.
    log-hazard ~ −0.5·F1 + 0.6·F2 + 0.3·F3 + ε

    Median OS by subtype (reflecting literature):
      PN  → ~26 months   CL  → ~14 months   MES → ~9 months
    """
    N    = Z.shape[0]
    F1   = Z[:, 0]   # Proneural
    F2   = Z[:, 1]   # Mesenchymal
    F3   = Z[:, 2]   # Proliferation

    # Reduced loadings so OS ranges are realistic: PN~26m, CL~15m, MES~8m
    log_hazard = -0.15*F1 + 0.20*F2 + 0.10*F3 + rng.normal(0, 0.3, N)
    scale      = np.exp(-log_hazard) * 18.5   # months scale
    os_time    = weibull_min.rvs(2.0, scale=scale, size=N, random_state=SEED)
    os_time    = np.clip(os_time, 0.3, 72)

    # ~25% censoring
    censored  = rng.uniform(0, 1, N) < 0.25
    cens_time = os_time * rng.uniform(0.3, 0.85, N)
    os_time_f = np.where(censored, cens_time, os_time)
    os_event  = (~censored).astype(int)

    age = rng.normal(58, 12, N).clip(18, 85)

    return pd.DataFrame({
        "patient_id":   [f"TCGA-{i:03d}" for i in range(N)],
        "subtype":      subtype_labels,
        "age_at_dx":    age.round(1),
        "os_months":    os_time_f.round(2),
        "os_event":     os_event,
        "karnofsky":    rng.choice([70, 80, 90, 100], N,
                                   p=[0.15, 0.35, 0.35, 0.15]),
        "surgery":      rng.choice(["GTR","STR","Biopsy"], N,
                                   p=[0.45, 0.40, 0.15]),
    })


# ═══════════════════════════════════════════════════════════════════════════════
#  4.  MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("  TCGA-GBM Multi-Omics Data Simulation")
    print("=" * 60)

    subtypes = np.array(["PN"]*N_PER + ["CL"]*N_PER + ["MES"]*N_PER)

    print(f"  Patients : {N}  ({N_PER} per subtype)")
    print(f"  Subtypes : Proneural | Classical | Mesenchymal")
    print(f"  Omics    : Expression({N_GENES}) · Methylation({N_CPGS}) · Mutations({N_MUTS})")

    # 1. Latent factors
    Z = make_latent_factors(N, rng)

    # 2. Name lists
    genes    = make_gene_names()[:N_GENES]
    cpgs     = make_cpg_names()[:N_CPGS]
    mut_genes = make_mut_names()[:N_MUTS]

    # 3. Simulate each modality
    print("\n  Simulating modalities…")
    expr   = simulate_expression(Z, rng, genes)
    meth   = simulate_methylation(Z, rng, cpgs)
    muts   = simulate_mutations(Z, rng, mut_genes, subtypes)
    clin   = simulate_clinical(Z, subtypes, rng)

    # 4. Save
    print("\n  Saving data…")
    expr.to_csv(os.path.join(DATA, "expression_matrix.csv"))
    meth.to_csv(os.path.join(DATA, "methylation_matrix.csv"))
    muts.to_csv(os.path.join(DATA, "mutation_matrix.csv"))
    clin.to_csv(os.path.join(DATA, "clinical_data.csv"), index=False)
    np.save(os.path.join(DATA, "true_factors.npy"), Z)

    # 5. QC summary
    print("\n  ── Data Summary ──────────────────────────────────")
    print(f"  Expression  : {expr.shape[0]} × {expr.shape[1]}  "
          f"| mean={expr.values.mean():.2f} | std={expr.values.std():.2f}")
    print(f"  Methylation : {meth.shape[0]} × {meth.shape[1]}  "
          f"| beta mean={meth.values.mean():.2f}")
    print(f"  Mutations   : {muts.shape[0]} × {muts.shape[1]}  "
          f"| mut rate={muts.values.mean():.3f}")
    print(f"  Clinical    : {clin.shape[0]} patients | "
          f"median OS={clin.os_months.median():.1f} months")
    print(f"\n  Subtype OS (median months):")
    for st in ["PN","CL","MES"]:
        m = clin.loc[clin.subtype==st, "os_months"].median()
        print(f"    {st:4s}: {m:.1f}")

    print("\n  ✓  Simulation complete.  Data saved to data/simulated/")


if __name__ == "__main__":
    main()
