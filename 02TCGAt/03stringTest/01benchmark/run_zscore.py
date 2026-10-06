"""
Gene set z-score scoring across samples.
Loops over all CPTAC cancers, saves zscore.csv per cancer.
"""

import os
import shutil
import sys

sys.path.insert(0, "/proj/c.zihao/work3/function")
from zscore import score_zscore

CANCERS = [
    "BLCA", "BRCA", "CRC", "GBM", "KIRC", "LGG", "LUAD", "LUSC",
    "MESO", "OV", "PAAD", "PRAD", "SKCM", "UVM",
]
EXPR_DIR = "/proj/c.zihao/work3/00data/TCGAt/exprset"
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
OUT_BASE = "/proj/c.zihao/work3/03TCGAt/03stringTest/bench/zscore"


def reset_out_dir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


reset_out_dir(OUT_BASE)

for cancer in CANCERS:
    expr_path = os.path.join(EXPR_DIR, f"{cancer}_exprSet_filtered.csv")
    out_dir = os.path.join(OUT_BASE, cancer)
    os.makedirs(out_dir, exist_ok=True)
    print(f"[{cancer}] zscore ...", flush=True)
    score_zscore(expr_path, CORUM, out_dir)

print("Done.")
