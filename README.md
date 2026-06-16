# Multi-Omics Integration in Glioblastoma (TCGA-GBM)
### MOFA+ · lifelines · scikit-learn

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![MOFA+](https://img.shields.io/badge/MOFA%2B-mofapy2%200.7.4-F28E2B)
![lifelines](https://img.shields.io/badge/lifelines-0.30-4E79A7)
![License](https://img.shields.io/badge/license-MIT-green)

> **Portfolio project** - Bioinformatics · Multi-omics integration · Survival analysis

---

## Overview

Glioblastoma Multiforme (GBM) is the most aggressive primary brain tumour in adults, with a median overall survival of ~15 months. It exhibits profound molecular heterogeneity across three well-characterised subtypes - **Proneural**, **Classical**, and **Mesenchymal** - each with distinct transcriptomic, epigenetic, and mutational landscapes (Verhaak et al., *Cancer Cell* 2010).

This project integrates three molecular modalities from TCGA-GBM (simulated with published biological parameters) using **Multi-Omics Factor Analysis Plus (MOFA+)** to:

- Discover latent factors that capture coordinated variation across modalities
- Link multi-omics signatures to patient survival endpoints
- Compare early fusion, late fusion, and MOFA+ integration for survival prediction

---

## Data & Biology

| Modality | Features | Biological content |
|----------|----------|-------------------|
| Gene expression | 340 genes | Oncogenes, tumour suppressors, proliferation markers, EMT factors |
| DNA methylation | 300 CpG sites | G-CIMP signature, MGMT promoter methylation |
| Somatic mutations | 50 genes | IDH1, EGFR, NF1, TP53, PTEN, ATRX, CDKN2A, TERT |

**N = 150 patients** (50 per subtype), simulated from published GBM biology:

| Subtype | Signature mutations | Median OS |
|---------|-------------------|-----------|
| Proneural (PN) | IDH1 (82%), ATRX (72%), TP53 (68%) | 19.6 months |
| Classical (CL) | EGFR (88%), TERT (72%), CDKN2A (52%) | 13.5 months |
| Mesenchymal (MES) | NF1 (42%), PTEN (54%), TP53 (62%) | 6.6 months |

---

## Pipeline

```
01_simulate_data.py       →  Simulate multi-omics matrices with GBM biology
02_mofa_integration.py    →  Fit MOFA+ model, extract 5 latent factors
03_visualize.py           →  Data overview, variance explained, UMAP, weights
04_survival_analysis.py   →  Kaplan-Meier curves, log-rank tests, Cox PH
05_fusion_comparison.py   →  Compare early / late / MOFA+ fusion by C-index
```

Run the full pipeline:
```bash
python run_pipeline.py
```

Or re-run from a specific step:
```bash
python run_pipeline.py --step 3
```

---

## Key Results

### MOFA+ Latent Factors

MOFA+ recovered **5 active factors** from 10 initialised (ARD sparsity prior):

| Factor | Top modality | Variance explained | Biological interpretation |
|--------|--------------|--------------------|--------------------------|
| Factor 1 | Methylation (70.5%) + Expression (65.2%) | Multi-modal | Proneural/Mesenchymal axis |
| Factor 2 | Expression (14.3%) + Mutations (9.3%) | Multi-modal | Mutation signature |
| Factor 3 | Methylation (18.5%) only | Uni-modal | G-CIMP / MGMT methylation |
| Factor 4 | Expression (5.5%) | Uni-modal | Expression module |
| Factor 5 | Minor signal | Uni-modal | Noise / minor variation |

**Total R²:** Expression 86.8% · Methylation 81.6% · Mutations 15.3%

### Survival Analysis

| Test | p-value | Significance |
|------|---------|-------------|
| Log-rank (subtype) | 1.7 × 10⁻¹⁵ | *** |
| Log-rank (Factor 1) | 1.7 × 10⁻⁸ | *** |
| Log-rank (Factor 3) | 3.4 × 10⁻² | * |
| Cox PH — Factor 1 (univariate) | HR = 0.63, p = 1.2 × 10⁻¹² | *** |

### Fusion Strategy Comparison (C-index, 5-fold CV)

| Method | C-index |
|--------|---------|
| Expression only | 0.690 ± 0.017 |
| Methylation only | 0.726 ± 0.026 |
| Mutations only | 0.602 ± 0.029 |
| Early fusion | 0.611 ± 0.035 |
| **MOFA+** | **0.688 ± 0.034** |
| Late fusion | **0.730** |

> Late fusion (per-view CPH ensemble) achieved the best C-index, while MOFA+ matched methylation-only performance with far fewer features (5 factors vs 300 CpG sites), demonstrating the value of interpretable dimensionality reduction.

---

## Figures

| Figure | Description |
|--------|-------------|
| `01_data_overview.png` | Feature distributions and subtype composition |
| `02_variance_explained.png` | R² per factor per modality (MOFA+ output) |
| `03_factor_umap.png` | UMAP of factor scores coloured by subtype |
| `04_factor_distributions.png` | Violin plots of top factors by subtype |
| `05_top_weights.png` | Top-weighted genes per factor |
| `06_km_curves.png` | Kaplan-Meier curves (subtype + factor terciles) |
| `07_cox_forest.png` | Univariate & multivariate Cox HR forest plot |
| `08_fusion_comparison.png` | C-index comparison across integration strategies |
| `09_modality_contribution.png` | Grouped bars: modality R² per factor |

---

## Setup

```bash
git clone https://github.com/<your-username>/tcga-gbm-multiomics.git
cd tcga-gbm-multiomics
pip install -r requirements.txt
python run_pipeline.py
```

**Python 3.10+ required.** All outputs regenerated in < 5 minutes on a standard laptop.

---

## References

- Verhaak RGW et al. (2010). Integrated genomic analysis identifies clinically relevant subtypes of glioblastoma. *Cancer Cell*, 17(1), 98–110.
- Argelaguet R et al. (2020). MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biology*, 21, 111.
- Davidson-Pilon C (2019). lifelines: survival analysis in Python. *Journal of Open Source Software*, 4(40), 1317.


