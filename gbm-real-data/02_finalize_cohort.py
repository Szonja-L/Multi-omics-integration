#!/usr/bin/env python3
"""
02_finalize_cohort.py
======================
Last recon step before the real data pull. Resolves the EXPRESSION_SUBTYPE
lookup (patient-level returned 0 rows, so we try sample-level), restricts
everything to primary tumor samples only (barcode suffix "-01") to avoid
double-counting patients with multiple tumor samples, and computes the
actual cohort sizes we'll be working with.

Usage:
    python 02_finalize_cohort.py

Copy the FULL output and send it back to Claude.
"""

import requests
from collections import Counter

BASE = "https://www.cbioportal.org/api"
STUDY_MAIN = "gbm_tcga"
STUDY_SUBTYPE = "gbm_tcga_pub2013"


def get(path, params=None):
    r = requests.get(f"{BASE}{path}", params=params, timeout=60)
    r.raise_for_status()
    return r.json()


def section(title):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def to_patient_id(sample_id: str) -> str:
    """TCGA-06-5415-01 -> TCGA-06-5415"""
    return "-".join(sample_id.split("-")[:3])


def main():
    # ── 1. Sample lists, filtered to primary tumor (-01) only ───────────
    section(f"STEP 1: Primary-tumor-only sample lists for '{STUDY_MAIN}'")
    sample_lists = get(f"/studies/{STUDY_MAIN}/sample-lists", params={"projection": "DETAILED"})
    wanted_lists = {
        "expression": "gbm_tcga_rna_seq_v2_mrna",
        "methylation": "gbm_tcga_methylation_hm450",
        "mutation": "gbm_tcga_sequenced",
    }
    patient_sets = {}
    for label, list_id in wanted_lists.items():
        sl = next((s for s in sample_lists if s["sampleListId"] == list_id), None)
        if sl is None:
            print(f"  !! Could not find sample list '{list_id}'")
            continue
        all_samples = sl.get("sampleIds", [])
        primary_only = [s for s in all_samples if s.endswith("-01")]
        patient_ids = {to_patient_id(s) for s in primary_only}
        patient_sets[label] = patient_ids
        print(f"  {label:<12}: {len(all_samples)} samples total -> "
              f"{len(primary_only)} primary-tumor -> {len(patient_ids)} unique patients")

    # ── 2. Overlaps (patient-level, primary tumor only) ─────────────────
    section("STEP 2: Patient-level overlaps (primary tumor only)")
    if len(patient_sets) == 3:
        e, m, u = patient_sets["expression"], patient_sets["methylation"], patient_sets["mutation"]
        print(f"  Expression only universe        : {len(e)}")
        print(f"  Methylation only universe       : {len(m)}")
        print(f"  Mutation only universe          : {len(u)}")
        print(f"  Expression ∩ Methylation         : {len(e & m)}")
        print(f"  Expression ∩ Mutation            : {len(e & u)}")
        print(f"  Methylation ∩ Mutation           : {len(m & u)}")
        print(f"  All three (complete case)        : {len(e & m & u)}")
        print(f"  Union (at least one view)        : {len(e | m | u)}")

    # ── 3. EXPRESSION_SUBTYPE: try SAMPLE level this time ───────────────
    section(f"STEP 3: EXPRESSION_SUBTYPE at SAMPLE level in '{STUDY_SUBTYPE}'")
    clin_sub_sample = get(
        f"/studies/{STUDY_SUBTYPE}/clinical-data",
        params={"clinicalDataType": "SAMPLE", "projection": "DETAILED",
                "pageSize": 100000, "pageNumber": 0},
    )
    print(f"  Total sample-level clinical records: {len(clin_sub_sample)}")
    subtype_rows = [r for r in clin_sub_sample if r["clinicalAttributeId"] == "EXPRESSION_SUBTYPE"]
    print(f"  Records with EXPRESSION_SUBTYPE: {len(subtype_rows)}")
    counts = Counter(r["value"] for r in subtype_rows)
    for val, cnt in counts.most_common():
        print(f"    {val:<20}: {cnt}")

    # Build patientId -> subtype map (primary tumor samples only)
    subtype_map = {}
    for r in subtype_rows:
        sid = r.get("sampleId", "")
        if sid.endswith("-01"):
            subtype_map[to_patient_id(sid)] = r["value"]
    print(f"\n  Usable patientId -> subtype mappings (primary tumor): {len(subtype_map)}")

    # If still empty, fall back and show what SAMPLE-level attributes DO exist
    if not subtype_rows:
        all_attr_ids = sorted({r["clinicalAttributeId"] for r in clin_sub_sample})
        print("\n  EXPRESSION_SUBTYPE not found at sample level either.")
        print("  Full list of SAMPLE-level attribute IDs actually present:")
        for a in all_attr_ids:
            print(f"    {a}")

    # ── 4. Cross-check: how many of our cohort patients have a subtype label? ─
    section("STEP 4: Final usable cohort (has data + has subtype + has OS)")
    if len(patient_sets) == 3 and subtype_map:
        union_patients = patient_sets["expression"] | patient_sets["methylation"] | patient_sets["mutation"]
        with_subtype = union_patients & subtype_map.keys()
        print(f"  Union cohort (>=1 omics view)        : {len(union_patients)}")
        print(f"  ...of those, with a subtype label    : {len(with_subtype)}")

        complete_case = patient_sets["expression"] & patient_sets["methylation"] & patient_sets["mutation"]
        complete_with_subtype = complete_case & subtype_map.keys()
        print(f"  Complete-case cohort (all 3 views)   : {len(complete_case)}")
        print(f"  ...of those, with a subtype label    : {len(complete_with_subtype)}")

        print(f"\n  Example subtype label values for our cohort:")
        sample_labels = [subtype_map[p] for p in list(with_subtype)[:10]]
        print(f"  {sample_labels}")

    print("\n" + "=" * 70)
    print("DONE. Copy everything above this line and send it back.")
    print("=" * 70)


if __name__ == "__main__":
    main()
