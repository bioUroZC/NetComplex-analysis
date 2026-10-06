"""
Static network-aware baseline using complex-internal PPI node strength
as weights for the mean of member gene z-scores.
Loops over all CPTAC cancers and saves netmean.csv per cancer.
"""

import os
import shutil
import sys

sys.path.insert(0, "/proj/c.zihao/work3/function")
from netmean import score_netmean


CANCERS = [
    "BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC",
]
EXPR_DIR = "/proj/c.zihao/work3/00data/cptacT/exprset"
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
LINKS = "/proj/c.zihao/work3/00data/HuRI/links.csv"
OUT_BASE = "/proj/c.zihao/work3/01protein/01HuriTest/bench/netmean"


def reset_out_dir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


reset_out_dir(OUT_BASE)

for cancer in CANCERS:
    expr_path = os.path.join(EXPR_DIR, f"{cancer}_exprSet_filtered.csv")
    out_dir = os.path.join(OUT_BASE, cancer)
    os.makedirs(out_dir, exist_ok=True)
    print(f"[{cancer}] netmean ...", flush=True)
    score_netmean(expr_path, CORUM, LINKS, out_dir)

print("Done.")
