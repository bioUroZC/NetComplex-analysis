#!/usr/bin/env python3
"""Protein-Min: minimum member-protein abundance per complex per sample.

Calls function/min.py's score_min() unchanged on the stable, zero-missing
protein matrix -- same real CORUM table used throughout the project.
"""

import sys
from pathlib import Path

sys.path.append("/proj/c.zihao/work3/function")
from min import score_min

CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
STABLE_DIR = Path("/proj/c.zihao/work3/01protein/02cor/01stable")
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
OUT_BASE = Path("/proj/c.zihao/work3/01protein/02cor/02scored/min")

for cancer in CANCERS:
    expr_path = STABLE_DIR / f"{cancer}_protein_stable.csv"
    out_dir = OUT_BASE / cancer
    print(f"[{cancer}] protein-min ...", flush=True)
    score_min(str(expr_path), CORUM, str(out_dir))

print("Done.")
