#!/usr/bin/env python3
"""
patch_fig5_fix.py
=================
Fixes the remaining KeyError: 'Factor3' in fig_top_weights.
The previous patches changed the configs rows to Factor2 but missed the
hardcoded ["Factor1", "Factor3"] list used for the row super-titles.

Run from repo root:
    python gbm-real-data\\patch_fig5_fix.py

Then re-run:
    python scripts\\03_visualize.py
"""

import os, re

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VIZ = os.path.join(REPO_ROOT, "scripts", "03_visualize.py")

with open(VIZ, "r", encoding="utf-8") as f:
    src = f.read()

original = src

# ── Fix 1: super-title loop list (pure ASCII, must always be fixable) ────────
OLD_LOOP = '["Factor1", "Factor3"]'
NEW_LOOP = '["Factor1", "Factor2"]'
if OLD_LOOP in src:
    src = src.replace(OLD_LOOP, NEW_LOOP, 1)
    print(f"  \u2713  Fixed row super-title loop: {OLD_LOOP} \u2192 {NEW_LOOP}")
else:
    print(f"  (loop list already patched or not found: {OLD_LOOP})")

# ── Fix 2: factor_labels dict — replace entire dict using regex so we avoid  ─
# ──         any Unicode encoding issues with the arrow characters             ─
#
# Strategy: find the dict by its unique surrounding context (pure ASCII) and
# rewrite it entirely, so we never need to match the Unicode arrow chars.
DICT_PATTERN = re.compile(
    r'factor_labels\s*=\s*\{[^}]+\}',
    re.DOTALL
)
NEW_DICT = (
    'factor_labels = {\n'
    '        "Factor1": "Factor 1  (Mutation landscape)",\n'
    '        "Factor2": "Factor 2  (Transcriptional subtype axis)",\n'
    '    }'
)
match = DICT_PATTERN.search(src)
if match:
    src = src[:match.start()] + NEW_DICT + src[match.end():]
    print(f"  \u2713  Replaced factor_labels dict with Factor1/Factor2 keys")
else:
    print("  !! Could not find factor_labels dict via regex — check file manually")

if src != original:
    with open(VIZ, "w", encoding="utf-8") as f:
        f.write(src)
    print(f"\n  \u2192 Saved {VIZ}")
else:
    print("\n  (no changes written)")

print("\nNow run:")
print("  python scripts\\03_visualize.py")
