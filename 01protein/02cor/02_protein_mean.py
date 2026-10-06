#!/usr/bin/env python3
"""Protein-Mean: mean member-protein abundance per complex per sample.

One of the nine protein-side references matching 01stringTest/01benchmark.
Calls function/mean.py's score_mean() unchanged on the stable, zero-missing
protein matrix from 01_stable_subset.py -- same real CORUM table used
throughout the project, no reimplementation.
"""

import sys
from pathlib import Path

sys.path.append("/proj/c.zihao/work3/function")
from mean import score_mean

CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
STABLE_DIR = Path("/proj/c.zihao/work3/01protein/02cor/01stable")
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
OUT_BASE = Path("/proj/c.zihao/work3/01protein/02cor/02scored/mean")

for cancer in CANCERS:
    expr_path = STABLE_DIR / f"{cancer}_protein_stable.csv"
    out_dir = OUT_BASE / cancer
    print(f"[{cancer}] protein-mean ...", flush=True)
    score_mean(str(expr_path), CORUM, str(out_dir))

print("Done.")
