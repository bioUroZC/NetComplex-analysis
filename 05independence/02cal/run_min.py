"""Minimum member gene expression on KIRC cohort-size subsets."""

import os
import shutil
import sys

sys.path.insert(0, "/proj/c.zihao/work3/function")
from min import score_min

LABELS = ["frac0.20", "frac0.40", "frac0.60", "frac0.80", "frac1.00"]
EXPR_DIR = "/proj/c.zihao/work3/05independence/01data"
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
OUT_BASE = "/proj/c.zihao/work3/05independence/results/bench/min"


def reset_out_dir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


reset_out_dir(OUT_BASE)

for label in LABELS:
    expr_path = os.path.join(EXPR_DIR, f"KIRC_{label}_exprSet.csv")
    out_dir = os.path.join(OUT_BASE, label)
    os.makedirs(out_dir, exist_ok=True)
    print(f"[{label}] min ...", flush=True)
    score_min(expr_path, CORUM, out_dir)

print("Done.")
