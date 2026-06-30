#!/usr/bin/env python3
"""
01_recon.py
===========
Second recon pass, now that we know which studies/profiles to use.

Checks, without pulling the full data yet:
  1. Sample-list overlap (which patients have expression / methylation /
     mutation data) in gbm_tcga
  2. Clinical data preview from gbm_tcga (OS_MONTHS, OS_STATUS, AGE,
     KARNOFSKY_PERFORMANCE_SCORE)
  3. EXPRESSION_SUBTYPE values from gbm_tcga_pub2013 (value counts)
  4. Patient-barcode overlap between gbm_tcga and gbm_tcga_pub2013
  5. ONE test gene (EGFR) pulled from expression + methylation + mutation
     profiles, to confirm the API call shape works before the big pull

Usage:
    python 01_recon.py

Copy the FULL output and send it back to Claude.
"""

import requests
import json

BASE = "https://www.cbioportal.org/api"

STUDY_MAIN = "gbm_tcga"            # expression, methylation, mutations, OS/AGE/KARNOFSKY
STUDY_SUBTYPE = "gbm_tcga_pub2013" # EXPRESSION_SUBTYPE labels


def get(path, params=None):
    r = requests.get(f"{BASE}{path}", params=params, timeout=60)
    r.raise_for_status()
    return r.json()


def post(path, body, params=None):
    r = requests.post(f"{BASE}{path}", json=body, params=params, timeout=60)
    return r  # NOT raising here — we want to inspect failures directly


def section(title):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def main():
    # ── 1. Sample lists per study ───────────────────────────────────────
    section(f"STEP 1: Sample lists for '{STUDY_MAIN}'")
    sample_lists = get(f"/studies/{STUDY_MAIN}/sample-lists", params={"projection": "DETAILED"})
    interesting = [sl for sl in sample_lists if any(
        k in sl["sampleListId"].lower() for k in
        ["rna_seq_v2_mrna", "methylation_hm450", "sequenced", "_all", "complete"]
    )]
    sample_id_sets = {}
    for sl in interesting:
        sids = sl.get("sampleIds", [])
        sample_id_sets[sl["sampleListId"]] = set(sids)
        print(f"  sampleListId : {sl['sampleListId']}")
        print(f"  name         : {sl.get('name')}")
        print(f"  sampleCount  : {len(sids)}")
        print("  " + "-" * 60)

    # Overlap across expression / methylation / mutation sample lists
    expr_list = next((k for k in sample_id_sets if "rna_seq_v2_mrna" in k), None)
    meth_list = next((k for k in sample_id_sets if "methylation_hm450" in k), None)
    mut_list = next((k for k in sample_id_sets if "sequenced" in k), None)

    if expr_list and meth_list and mut_list:
        overlap = sample_id_sets[expr_list] & sample_id_sets[meth_list] & sample_id_sets[mut_list]
        print(f"\n  >>> Samples with ALL THREE (expr + meth + mut): {len(overlap)}")
        print(f"  >>> Example sample IDs: {list(overlap)[:5]}")
    else:
        print("\n  Could not auto-identify all three relevant sample lists by name —")
        print("  inspect the printed list above manually.")
        overlap = set()

    # ── 2. Clinical data preview: gbm_tcga ──────────────────────────────
    section(f"STEP 2: Clinical data preview for '{STUDY_MAIN}'")
    clin_main = get(
        f"/studies/{STUDY_MAIN}/clinical-data",
        params={"clinicalDataType": "PATIENT", "projection": "DETAILED",
                "pageSize": 100000, "pageNumber": 0},
    )
    print(f"  Total clinical-data records (patient-level): {len(clin_main)}")
    wanted = {"OS_MONTHS", "OS_STATUS", "AGE", "KARNOFSKY_PERFORMANCE_SCORE"}
    sample_rows = [r for r in clin_main if r["clinicalAttributeId"] in wanted][:12]
    for r in sample_rows:
        print(f"  patientId={r['patientId']:<14} attr={r['clinicalAttributeId']:<28} value={r['value']}")

    os_status_values = sorted({r["value"] for r in clin_main if r["clinicalAttributeId"] == "OS_STATUS"})
    print(f"\n  Unique OS_STATUS values: {os_status_values}")

    # ── 3. EXPRESSION_SUBTYPE from gbm_tcga_pub2013 ─────────────────────
    section(f"STEP 3: EXPRESSION_SUBTYPE values from '{STUDY_SUBTYPE}'")
    clin_sub = get(
        f"/studies/{STUDY_SUBTYPE}/clinical-data",
        params={"clinicalDataType": "PATIENT", "projection": "DETAILED",
                "pageSize": 100000, "pageNumber": 0},
    )
    subtype_rows = [r for r in clin_sub if r["clinicalAttributeId"] == "EXPRESSION_SUBTYPE"]
    print(f"  Patients with EXPRESSION_SUBTYPE recorded: {len(subtype_rows)}")
    from collections import Counter
    counts = Counter(r["value"] for r in subtype_rows)
    for val, cnt in counts.most_common():
        print(f"    {val:<20}: {cnt}")

    # ── 4. Patient barcode overlap between the two studies ──────────────
    section("STEP 4: Patient barcode overlap between the two studies")
    patients_main = {r["patientId"] for r in clin_main}
    patients_sub = {r["patientId"] for r in clin_sub}
    common_patients = patients_main & patients_sub
    print(f"  Patients in {STUDY_MAIN}        : {len(patients_main)}")
    print(f"  Patients in {STUDY_SUBTYPE} : {len(patients_sub)}")
    print(f"  Patients in BOTH (barcode match)  : {len(common_patients)}")
    print(f"  Example shared patient IDs        : {list(common_patients)[:5]}")
    print(f"  Example {STUDY_MAIN} patientId format     : {list(patients_main)[:3]}")
    print(f"  Example {STUDY_SUBTYPE} patientId format : {list(patients_sub)[:3]}")

    # ── 5. Test pulling ONE gene's molecular data (EGFR) ────────────────
    section("STEP 5: Test molecular-data fetch for gene EGFR")

    # 5a. Resolve Hugo symbol -> Entrez Gene ID
    gene_resp = post("/genes/fetch", body=["EGFR"], params={"geneIdType": "HUGO_GENE_SYMBOL"})
    print(f"  /genes/fetch status: {gene_resp.status_code}")
    print(f"  /genes/fetch body  : {gene_resp.text[:500]}")

    if gene_resp.status_code == 200 and gene_resp.json():
        entrez_id = gene_resp.json()[0]["entrezGeneId"]
        print(f"\n  EGFR entrezGeneId = {entrez_id}")

        test_samples = list(overlap)[:5] if overlap else list(patients_main)[:5]

        for profile_suffix, label in [
            ("rna_seq_v2_mrna", "EXPRESSION"),
            ("methylation_hm450", "METHYLATION"),
        ]:
            profile_id = f"{STUDY_MAIN}_{profile_suffix}"
            body = {"entrezGeneIds": [entrez_id], "sampleIds": test_samples}
            resp = post(f"/molecular-profiles/{profile_id}/molecular-data/fetch", body=body)
            print(f"\n  [{label}] POST /molecular-profiles/{profile_id}/molecular-data/fetch")
            print(f"  status: {resp.status_code}")
            print(f"  body  : {resp.text[:600]}")
    else:
        print("  Could not resolve EGFR to an Entrez ID — inspect response above.")

    # 5b. Test mutation fetch for EGFR
    mut_profile_id = f"{STUDY_MAIN}_mutations"
    if gene_resp.status_code == 200 and gene_resp.json():
        entrez_id = gene_resp.json()[0]["entrezGeneId"]
        body = {"entrezGeneIds": [entrez_id], "sampleIds": list(overlap)[:5] if overlap else list(patients_main)[:5]}
        resp = post(f"/molecular-profiles/{mut_profile_id}/mutations/fetch", body=body)
        print(f"\n  [MUTATIONS] POST /molecular-profiles/{mut_profile_id}/mutations/fetch")
        print(f"  status: {resp.status_code}")
        print(f"  body  : {resp.text[:600]}")

    print("\n" + "=" * 70)
    print("DONE. Copy everything above this line and send it back.")
    print("=" * 70)


if __name__ == "__main__":
    main()
