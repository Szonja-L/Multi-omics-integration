#!/usr/bin/env python3
"""
00_discover_cbioportal.py
==========================
Run this FIRST, before any real-data pull. It only reads metadata from the
public cBioPortal REST API (no auth needed) and prints what's available, so
the next script (real data extraction) uses correct study/profile/attribute
IDs instead of guessed ones.

Usage:
    pip install requests
    python 00_discover_cbioportal.py

Copy the FULL printed output and send it back to Claude.
"""

import requests
import sys

BASE = "https://www.cbioportal.org/api"


def get(path, params=None):
    url = f"{BASE}{path}"
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    return r.json()


def main():
    # ── 1. Find GBM studies ─────────────────────────────────────────────
    print("=" * 70)
    print("STEP 1: Searching for GBM studies on cBioPortal")
    print("=" * 70)
    try:
        all_studies = get("/studies")
    except Exception as e:
        print(f"  ERROR fetching studies list: {e}")
        sys.exit(1)

    gbm_studies = [
        s for s in all_studies
        if "gbm" in s.get("studyId", "").lower()
        or "glioblastoma" in s.get("name", "").lower()
    ]

    if not gbm_studies:
        print("  No GBM studies found — cBioPortal API may have changed.")
        sys.exit(1)

    print(f"  Found {len(gbm_studies)} GBM-related studies:\n")
    for s in gbm_studies:
        print(f"  studyId        : {s.get('studyId')}")
        print(f"  name           : {s.get('name')}")
        print(f"  description    : {(s.get('description') or '')[:160]}")
        print(f"  allSampleCount : {s.get('allSampleCount')}")
        print("  " + "-" * 60)

    # Candidates we care most about — adjust if the printed list above differs
    candidates = ["gbm_tcga_pan_can_atlas_2018", "gbm_tcga", "gbm_tcga_pub2013"]
    target_studies = [c for c in candidates if c in [s["studyId"] for s in gbm_studies]]

    if not target_studies:
        print("\n  None of the expected candidate study IDs were found.")
        print("  Use one of the studyId values printed above instead.")
        target_studies = [gbm_studies[0]["studyId"]]

    # ── 2. For each candidate study, list molecular profiles ───────────
    for study_id in target_studies:
        print("\n" + "=" * 70)
        print(f"STEP 2: Molecular profiles for '{study_id}'")
        print("=" * 70)
        try:
            profiles = get(f"/studies/{study_id}/molecular-profiles")
        except Exception as e:
            print(f"  ERROR fetching molecular profiles: {e}")
            continue

        for p in profiles:
            print(f"  molecularProfileId      : {p.get('molecularProfileId')}")
            print(f"  molecularAlterationType : {p.get('molecularAlterationType')}")
            print(f"  datatype                : {p.get('datatype')}")
            print(f"  name                    : {p.get('name')}")
            print("  " + "-" * 60)

        # ── 3. Clinical attributes for the same study ───────────────────
        print(f"\nSTEP 3: Clinical attributes for '{study_id}'")
        print("-" * 70)
        try:
            attrs = get(f"/studies/{study_id}/clinical-attributes")
        except Exception as e:
            print(f"  ERROR fetching clinical attributes: {e}")
            continue

        # Highlight the ones we actually need
        keywords = ["OS_", "SURVIVAL", "SUBTYPE", "AGE", "KARNOFSKY",
                    "SURGERY", "STATUS", "IDH", "MGMT", "GRADE"]
        relevant = [a for a in attrs
                    if any(k in a.get("clinicalAttributeId", "").upper() for k in keywords)]

        print(f"  {len(attrs)} total attributes. Relevant-looking ones:\n")
        for a in relevant:
            print(f"  clinicalAttributeId : {a.get('clinicalAttributeId')}")
            print(f"  displayName         : {a.get('displayName')}")
            print(f"  datatype            : {a.get('datatype')}")
            print("  " + "-" * 60)

        # ── 4. Sample counts per molecular profile (rough overlap check) ─
        print(f"\nSTEP 4: Sample counts per profile for '{study_id}'")
        print("-" * 70)
        try:
            samples = get(f"/studies/{study_id}/samples")
            print(f"  Total samples in study: {len(samples)}")
        except Exception as e:
            print(f"  ERROR fetching sample list: {e}")

    print("\n" + "=" * 70)
    print("DONE. Copy everything above this line and send it back.")
    print("=" * 70)


if __name__ == "__main__":
    main()
