#!/usr/bin/env python3
"""
patch_unicode_fix.py
====================
Fixes the 4 cosmetic patches in 03_visualize.py that the first patch script
missed because the file contains literal Unicode characters (₂, ↔, –) while
the previous script searched for escape sequences.

Run from the repo root:
    python gbm-real-data\\patch_unicode_fix.py

Then re-run:
    python scripts\\03_visualize.py
"""

import os

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VIZ = os.path.join(REPO_ROOT, "scripts", "03_visualize.py")

with open(VIZ, "r", encoding="utf-8") as f:
    src = f.read()

original = src
changes = 0

PATCHES = [
    # 1. Expression heatmap xlabel: log₂ TPM → log₂ RSEM+1
    # The ₂ character (U+2082) appears literally in the file
    (
        "Gene Expression\n(log\u2082 TPM)",
        "Gene Expression\n(log\u2082 RSEM+1)",
    ),

    # 2. UMAP right-panel title: Factor 1 description
    # Uses – (U+2013 en-dash) and ↔ (U+2194)
    (
        "MOFA+ Factor Space \u2013 Factor 1 (Proneural \u2194 Mesenchymal)",
        "MOFA+ Factor Space \u2013 Factor 1 (Mutation landscape)",
    ),

    # 3. Violin-plot factor_titles list (uses ↔ U+2194 and / characters)
    (
        '"Factor 1\\nProneural \u2194 Mes",\n'
        '        "Factor 2\\nMutation signature",\n'
        '        "Factor 3\\nG-CIMP / MGMT",\n'
        '        "Factor 4\\nExpression module",\n'
        '        "Factor 5\\nMinor signal",',
        '"Factor 1\\nMutation landscape",\n'
        '        "Factor 2\\nTranscriptional\\nsubtype",\n'
        '        "Factor 3\\nExpression\\nmodule A",\n'
        '        "Factor 4\\nExpression\\nmodule B",\n'
        '        "Factor 5\\nMethylation\\naxis",',
    ),

    # 4. Top-weights factor_labels dict (uses ↔ U+2194)
    (
        '"Factor1": "Factor 1  (Proneural \u2194 Mesenchymal axis)",\n'
        '        "Factor3": "Factor 3  (G-CIMP / MGMT methylation)",',
        '"Factor1": "Factor 1  (Mutation landscape)",\n'
        '        "Factor2": "Factor 2  (Transcriptional subtype axis)",',
    ),
]

for old, new in PATCHES:
    if old in src:
        src = src.replace(old, new, 1)
        changes += 1
        print(f"  \u2713  Patched: {repr(old[:60])} \u2026")
    else:
        print(f"  !!  Could not find: {repr(old[:60])}")
        # Debug: show surrounding characters as hex to spot encoding issues
        # Find approximate location by searching for a stable nearby substring
        probe = old[:20]
        idx = src.find(probe)
        if idx >= 0:
            window = src[max(0, idx-5):idx+80]
            print(f"      Nearby text (hex): {window.encode('utf-8').hex()}")

if src != original:
    with open(VIZ, "w", encoding="utf-8") as f:
        f.write(src)
    print(f"\n  \u2192 Saved {VIZ} ({changes} patches applied)")
else:
    print("\n  (no changes written \u2014 all patches may already be applied or not found)")

print("\nDone. Re-run:")
print("  python scripts\\03_visualize.py")
