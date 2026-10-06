#!/usr/bin/env python3
"""Protein-NetMean: degree-weighted mean of member-protein z-scores.

Calls function/netmean.py's score_netmean() unchanged on the stable,
zero-missing protein matrix. Uses the same STRING links file as the
RNA-side benchmark; edges are unweighted inside score_netmean().
"""

import sys
from pathlib import Path

sys.path.append("/proj/c.zihao/work3/function")
from netmean import score_netmean

CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
STABLE_DIR = Path("/proj/c.zihao/work3/01protein/02cor/01stable")
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
LINKS = "/proj/c.zihao/work3/00data/string/links.csv"
OUT_BASE = Path("/proj/c.zihao/work3/01protein/02cor/02scored/netmean")

for cancer in CANCERS:
    expr_path = STABLE_DIR / f"{cancer}_protein_stable.csv"
    out_dir = OUT_BASE / cancer
    print(f"[{cancer}] protein-netmean ...", flush=True)
    score_netmean(str(expr_path), CORUM, LINKS, str(out_dir))

print("Done.")
