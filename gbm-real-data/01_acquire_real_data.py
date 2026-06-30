#!/usr/bin/env python3
"""
01_acquire_real_data.py
========================
Replaces 01_simulate_data.py. Pulls REAL TCGA-GBM data from the public
cBioPortal API for the complete-case cohort (patients with expression +
methylation + mutation data, a valid molecular subtype, and survival data)
and writes the same four CSVs the existing pipeline already expects:

    data/simulated/expression_matrix.csv   (patients x genes, RNA-seq RSEM)
    data/simulated/methylation_matrix.csv  (patients x genes, HM450 beta, gene-level)
    data/simulated/mutation_matrix.csv     (patients x genes, binary somatic)
    data/simulated/clinical_data.csv       (patient_id, subtype, age_at_dx,
                                             os_months, os_event, karnofsky)

NOTE ON FOLDER NAME: still writing into data/simulated/ on purpose, so that
02_mofa_integration.py, 03_visualize.py, 04_survival_analysis.py and
05_fusion_comparison.py need ZERO code changes. We'll rename the folder to
data/real/ (and update those 4 path constants) once results are verified —
that's a pure rename, not a logic change.

Data sources:
  gbm_tcga          (TCGA, Firehose Legacy) -> expression, methylation,
                     mutations, AGE, KARNOFSKY_PERFORMANCE_SCORE, OS_MONTHS,
                     OS_STATUS
  gbm_tcga_pub2013  (TCGA, Cell 2013)       -> EXPRESSION_SUBTYPE (real
                     Proneural/Classical/Mesenchymal/Neural/G-CIMP labels)

Subtype harmonization (documented, citable choices):
  - "Neural" dropped — Wang et al. 2017 (Cancer Cell) showed this subtype is
    likely driven by contaminating normal brain tissue; removed from the
    consensus classification.
  - "G-CIMP" merged into "Proneural" — G-CIMP tumors are IDH-mutant and
    transcriptionally Proneural (Noushmehr et al. 2010; Brennan et al. 2013).

Known limitation (documented for the README): "methylation" here is
cBioPortal's gene-level HM450 summary, not raw CpG-probe-level data, because
the standard molecular-data API returns one value per gene. TERT is kept in
the mutation panel but standard exome MAF calling often misses TERT promoter
mutations, so expect a near-zero rate for it here (a known TCGA limitation,
not a bug).

Usage:
    python 01_acquire_real_data.py

Copy the FULL output and send it back to Claude before re-running 02-05.
"""

import os
import requests
import pandas as pd
import numpy as np
from collections import Counter

BASE_API = "https://www.cbioportal.org/api"
STUDY_MAIN = "gbm_tcga"
STUDY_SUBTYPE = "gbm_tcga_pub2013"

# ── EDIT THIS if your local repo root differs ────────────────────────────────
REPO_ROOT = r"C:\Users\Szonja\Desktop\GitHub\8. Multi-omics-integration"
DATA_DIR = os.path.join(REPO_ROOT, "data", "simulated")

# ── Curated GBM gene panel (same real genes as 01_simulate_data.py, minus
#    the synthetic GENE_### noise placeholders, which aren't real symbols) ──
PRONEURAL = ["IDH1","IDH2","PDGFRA","OLIG2","OLIG1","SOX2","SOX11","TCF7L2",
    "DLL3","NOTCH1","SEMA3C","NXPH1","CADPS","DSCAML1","PRICKLE1","PTPRD",
    "NRXN1","NR2E1","GPR17","PLCL1","CAMK2A","SHANK2","SYT1","GABRA1",
    "GRIN1","ASCL1","DNMT3A","LDHB","ATP1A2","RAB3A","NR4A1","DPYSL3",
    "SLC17A6","STMN2","GAD1","GAD2","MAP2","TUBB3","SCN2A","SNAP25"]
CLASSICAL = ["EGFR","CDK4","MDM2","MDM4","CCND1","CDK6","E2F3","GLI1","SHH",
    "PTCH1","VEGFA","FGF2","FGFR1","IGFBP2","NES","FBLN1","COL4A1","LAMC1",
    "SPARC","HIF1A","CA9","CXCL12","PDGFA","TGFB1","TGFB2","MMP9","MMP2",
    "TIMP1","TWIST1","SNAI1","E2F2","TGFB3","ITGB3","THBS2","TNFRSF12A",
    "LTBP1","LTBP2","FSCN1","CXCL1"]
MESENCHYMAL = ["CHI3L1","MET","CD44","FN1","SERPINE1","LGALS3","ACTN1","VIM",
    "TNFRSF1A","TRADD","RELB","CASP1","CASP3","BCL2A1","TNC","COL5A2",
    "CXCR4","PLAUR","TRAF3","CD163","BMP2","SLC2A3","LOXL2","NRP1","CTGF",
    "CYR61","ANXA2","SPP1","ITGA3","IQGAP1","MYH9","ARHGEF2","ROCK1","SNAI2",
    "ZEB1","ZEB2","COL1A1","COL3A1","MXRA8"]
PROLIFERATION = ["MKI67","TOP2A","PCNA","MCM2","MCM6","CCNB1","CCNB2","CDK1",
    "AURKB","BUB1","CCNE1","CDC20","CDC25C","H2AFX","INCENP","PTTG1","SMC4",
    "TPX2","BIRC5","PLK1","CENPF","KIF11","KIF2C","MELK","TACC3","UHRF1",
    "CKS2","FOXM1","E2F1","RRM2"]
SUPPRESSOR = ["PTEN","RB1","ATRX","TP53","NF1","CDKN2A","CDKN2B","TSC1","TSC2",
    "FBXW7","AKT1","PIK3CA","PIK3R1","MTOR","PARP1","ATR","CHK1","CHK2",
    "RAD51","BRCA1","BRCA2","STAG2","MGMT","MLH1","MSH2","ALDH1A1","SOX9",
    "GFAP","S100B","NESTIN"]
EXPR_METH_PANEL = sorted(set(PRONEURAL + CLASSICAL + MESENCHYMAL + PROLIFERATION + SUPPRESSOR))

MUT_PANEL = ["IDH1","IDH2","EGFR","NF1","PTEN","TP53","ATRX","RB1","PIK3CA",
    "PIK3R1","TERT","CDK4","CDKN2A","CDKN2B","PDGFRA","STAG2","FUBP1","CIC",
    "NOTCH1","PTCH1","GLI1","AKT1","AKT2","CCND2","RPS6KA3","BRAF","KRAS",
    "NRAS","HRAS","MET","ERBB2","ERBB3","ERBB4","FGFR1","FGFR3","ALK","RET",
    "MYC","MYCN","MDM2","MDM4","BCL2","MCL1","VHL","SMARCA4","ARID1A",
    "KMT2D","KMT2C","SETD2","EP300"]

SUBTYPE_MAP = {"Proneural": "PN", "G-CIMP": "PN", "Classical": "CL", "Mesenchymal": "MES"}
# "Neural" deliberately absent -> dropped


def get(path, params=None):
    r = requests.get(f"{BASE_API}{path}", params=params, timeout=60)
    r.raise_for_status()
    return r.json()


def post(path, body, params=None):
    r = requests.post(f"{BASE_API}{path}", json=body, params=params, timeout=120)
    r.raise_for_status()
    return r.json()


def section(title):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def to_patient_id(sample_id: str) -> str:
    return "-".join(sample_id.split("-")[:3])


# ═══════════════════════════════════════════════════════════════════════════
# 1. BUILD COHORT
# ═══════════════════════════════════════════════════════════════════════════

def build_cohort():
    section("STEP 1: Building the complete-case cohort")

    sample_lists = get(f"/studies/{STUDY_MAIN}/sample-lists", params={"projection": "DETAILED"})
    list_ids = {"expression": "gbm_tcga_rna_seq_v2_mrna",
                "methylation": "gbm_tcga_methylation_hm450",
                "mutation": "gbm_tcga_sequenced"}
    patient_sets = {}
    for label, lid in list_ids.items():
        sl = next(s for s in sample_lists if s["sampleListId"] == lid)
        primary = [s for s in sl["sampleIds"] if s.endswith("-01")]
        patient_sets[label] = {to_patient_id(s) for s in primary}

    complete_case = patient_sets["expression"] & patient_sets["methylation"] & patient_sets["mutation"]
    print(f"  Complete-case (expr+meth+mut) patients: {len(complete_case)}")

    # Subtype labels (sample-level, gbm_tcga_pub2013)
    clin_sub = get(f"/studies/{STUDY_SUBTYPE}/clinical-data",
                   params={"clinicalDataType": "SAMPLE", "projection": "DETAILED",
                           "pageSize": 100000, "pageNumber": 0})
    raw_subtype = {}
    for r in clin_sub:
        if r["clinicalAttributeId"] == "EXPRESSION_SUBTYPE" and r["sampleId"].endswith("-01"):
            raw_subtype[to_patient_id(r["sampleId"])] = r["value"]

    subtype_map = {pid: SUBTYPE_MAP[val] for pid, val in raw_subtype.items() if val in SUBTYPE_MAP}
    dropped_neural = sum(1 for v in raw_subtype.values() if v == "Neural")
    print(f"  Patients with a usable subtype (Neural dropped, G-CIMP merged): {len(subtype_map)}")
    print(f"  ({dropped_neural} patients excluded for being 'Neural' subtype)")

    # Clinical / survival data (patient-level, gbm_tcga)
    clin_main = get(f"/studies/{STUDY_MAIN}/clinical-data",
                    params={"clinicalDataType": "PATIENT", "projection": "DETAILED",
                            "pageSize": 100000, "pageNumber": 0})
    clin_records = {}
    for r in clin_main:
        pid, attr, val = r["patientId"], r["clinicalAttributeId"], r["value"]
        clin_records.setdefault(pid, {})[attr] = val

    has_os = {pid for pid, d in clin_records.items()
              if "OS_MONTHS" in d and "OS_STATUS" in d}

    final_cohort = sorted(complete_case & subtype_map.keys() & has_os)
    print(f"  FINAL cohort (complete-case + subtype + OS data): {len(final_cohort)}")

    subtype_counts = Counter(subtype_map[p] for p in final_cohort)
    print(f"  Subtype breakdown: {dict(subtype_counts)}")

    return final_cohort, subtype_map, clin_records


# ═══════════════════════════════════════════════════════════════════════════
# 2. GENE RESOLUTION
# ═══════════════════════════════════════════════════════════════════════════

def resolve_genes(symbols):
    resp = post("/genes/fetch", body=symbols, params={"geneIdType": "HUGO_GENE_SYMBOL"})
    symbol_to_entrez = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in resp}
    entrez_to_symbol = {v: k for k, v in symbol_to_entrez.items()}
    missing = sorted(set(symbols) - symbol_to_entrez.keys())
    if missing:
        print(f"  !! {len(missing)} symbols did not resolve and will be skipped: {missing}")
    return symbol_to_entrez, entrez_to_symbol


# ═══════════════════════════════════════════════════════════════════════════
# 3. CONTINUOUS DATA (expression / methylation)
# ═══════════════════════════════════════════════════════════════════════════

def fetch_continuous_matrix(profile_id, entrez_to_symbol, sample_ids):
    body = {"entrezGeneIds": list(entrez_to_symbol.keys()), "sampleIds": sample_ids}
    records = post(f"/molecular-profiles/{profile_id}/molecular-data/fetch", body)
    df = pd.DataFrame.from_records(records)
    if df.empty:
        raise RuntimeError(f"No data returned for profile {profile_id}")
    df["patient_id"] = df["sampleId"].apply(to_patient_id)
    df["gene"] = df["entrezGeneId"].map(entrez_to_symbol)
    wide = df.pivot_table(index="patient_id", columns="gene", values="value", aggfunc="first")
    return wide


# ═══════════════════════════════════════════════════════════════════════════
# 4. MUTATIONS (binary)
# ═══════════════════════════════════════════════════════════════════════════

def fetch_mutation_matrix(profile_id, entrez_to_symbol, sample_ids, patient_ids):
    body = {"entrezGeneIds": list(entrez_to_symbol.keys()), "sampleIds": sample_ids}
    records = post(f"/molecular-profiles/{profile_id}/mutations/fetch", body)

    mut = pd.DataFrame(0, index=patient_ids, columns=list(entrez_to_symbol.values()))
    n_calls = 0
    for r in records:
        if r.get("mutationStatus") != "Somatic":
            continue
        pid = to_patient_id(r["sampleId"])
        gene = entrez_to_symbol.get(r["entrezGeneId"])
        if pid in mut.index and gene is not None:
            mut.loc[pid, gene] = 1
            n_calls += 1
    print(f"  {n_calls} somatic mutation calls applied across the panel")
    return mut


# ═══════════════════════════════════════════════════════════════════════════
# 5. CLINICAL TABLE
# ═══════════════════════════════════════════════════════════════════════════

def build_clinical_df(final_cohort, subtype_map, clin_records):
    rows = []
    for pid in final_cohort:
        d = clin_records[pid]
        os_status = d.get("OS_STATUS", "")
        rows.append({
            "patient_id": pid,
            "subtype": subtype_map[pid],
            "age_at_dx": float(d["AGE"]) if "AGE" in d else np.nan,
            "os_months": float(d["OS_MONTHS"]),
            "os_event": 1 if os_status.startswith("1") else 0,
            "karnofsky": float(d["KARNOFSKY_PERFORMANCE_SCORE"]) if "KARNOFSKY_PERFORMANCE_SCORE" in d else np.nan,
        })
    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("  TCGA-GBM REAL DATA ACQUISITION  (cBioPortal)")
    print("=" * 70)

    final_cohort, subtype_map, clin_records = build_cohort()
    sample_ids = [f"{pid}-01" for pid in final_cohort]

    section("STEP 2: Resolving gene panels to Entrez IDs")
    print(f"  Expression/methylation panel: {len(EXPR_METH_PANEL)} genes")
    em_sym2ez, em_ez2sym = resolve_genes(EXPR_METH_PANEL)
    print(f"  Mutation panel: {len(MUT_PANEL)} genes")
    mut_sym2ez, mut_ez2sym = resolve_genes(MUT_PANEL)

    section("STEP 3: Pulling expression data (RNA-seq RSEM)")
    expr = fetch_continuous_matrix("gbm_tcga_rna_seq_v2_mrna", em_ez2sym, sample_ids)
    expr = expr.reindex(final_cohort)
    # Raw RSEM is unbounded and right-skewed (e.g. EGFR ~41,000). MOFA+'s
    # Gaussian likelihood and the heatmap color scale in 03_visualize.py both
    # assume a log2 scale (the synthetic version generated data already on
    # that scale) -> standard RNA-seq log2(RSEM+1) transform, applied once here.
    print(f"  Raw RSEM range before log2: min={expr.min().min():.1f}, max={expr.max().max():.1f}")
    expr = np.log2(expr.clip(lower=0) + 1.0)
    print(f"  Log2(RSEM+1) range: min={expr.min().min():.2f}, max={expr.max().max():.2f}")
    print(f"  Expression matrix: {expr.shape}, NaNs: {expr.isna().sum().sum()}")

    section("STEP 4: Pulling methylation data (HM450, gene-level)")
    meth = fetch_continuous_matrix("gbm_tcga_methylation_hm450", em_ez2sym, sample_ids)
    meth = meth.reindex(final_cohort)
    print(f"  Methylation matrix: {meth.shape}, NaNs: {meth.isna().sum().sum()}")

    section("STEP 5: Pulling mutation data (binary, somatic)")
    muts = fetch_mutation_matrix("gbm_tcga_mutations", mut_ez2sym, sample_ids, final_cohort)
    zero_var_genes = muts.columns[muts.var(axis=0) == 0].tolist()
    if zero_var_genes:
        print(f"  Dropping {len(zero_var_genes)} genes with zero mutations observed "
              f"in this cohort: {zero_var_genes}")
        muts = muts.drop(columns=zero_var_genes)
    print(f"  Mutation matrix: {muts.shape}")

    section("STEP 6: Building clinical table")
    clin = build_clinical_df(final_cohort, subtype_map, clin_records)
    print(clin.describe(include="all").to_string())

    section("STEP 7: Saving CSVs")
    os.makedirs(DATA_DIR, exist_ok=True)
    expr.to_csv(os.path.join(DATA_DIR, "expression_matrix.csv"))
    meth.to_csv(os.path.join(DATA_DIR, "methylation_matrix.csv"))
    muts.to_csv(os.path.join(DATA_DIR, "mutation_matrix.csv"))
    clin.to_csv(os.path.join(DATA_DIR, "clinical_data.csv"), index=False)
    print(f"  Saved 4 CSVs to: {DATA_DIR}")

    section("STEP 8: Sanity checks (compare to known GBM biology)")
    for gene in ["IDH1", "EGFR", "NF1", "TP53", "TERT"]:
        if gene in muts.columns:
            rate_by_subtype = muts[gene].groupby(clin.set_index("patient_id")["subtype"]).mean()
            print(f"  {gene:<8} mutation rate by subtype: {rate_by_subtype.to_dict()}")

    print("\n" + "=" * 70)
    print("DONE. Copy everything above this line and send it back.")
    print("=" * 70)


if __name__ == "__main__":
    main()
