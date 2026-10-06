#!/usr/bin/env python3
"""Protein-PerPAS-adapted: topology-weighted local protein activity.

Calls function/perpas_adapted.py's score_perpas_adapted() unchanged on
the stable, zero-missing protein matrix. Uses the same STRING links
file as the RNA-side benchmark; edges are unweighted inside the scorer.
"""

import sys
from pathlib import Path

sys.path.append("/proj/c.zihao/work3/function")
from perpas_adapted import score_perpas_adapted

CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
STABLE_DIR = Path("/proj/c.zihao/work3/01protein/02cor/01stable")
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
LINKS = "/proj/c.zihao/work3/00data/string/links.csv"
OUT_BASE = Path("/proj/c.zihao/work3/01protein/02cor/02scored/perpas_adapted")

for cancer in CANCERS:
    expr_path = STABLE_DIR / f"{cancer}_protein_stable.csv"
    out_dir = OUT_BASE / cancer
    print(f"[{cancer}] protein-perpas_adapted ...", flush=True)
    score_perpas_adapted(str(expr_path), CORUM, LINKS, str(out_dir))

print("Done.")
