#!/usr/bin/env python3
"""Protein-Zscore: combined z-score of member-protein abundance per complex.

Calls function/zscore.py's score_zscore() unchanged on the stable,
zero-missing protein matrix -- same real CORUM table used throughout
the project.
"""

import sys
from pathlib import Path

sys.path.append("/proj/c.zihao/work3/function")
from zscore import score_zscore

CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
STABLE_DIR = Path("/proj/c.zihao/work3/01protein/02cor/01stable")
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
OUT_BASE = Path("/proj/c.zihao/work3/01protein/02cor/02scored/zscore")

for cancer in CANCERS:
    expr_path = STABLE_DIR / f"{cancer}_protein_stable.csv"
    out_dir = OUT_BASE / cancer
    print(f"[{cancer}] protein-zscore ...", flush=True)
    score_zscore(str(expr_path), CORUM, str(out_dir))

print("Done.")
