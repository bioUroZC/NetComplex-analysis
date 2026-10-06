"""
PerPAS-inspired topology baseline adapted to undirected PPI-defined
protein complexes.
Loops over GSE96583 samples and saves perpas_adapted.csv per sample.
"""

import os
import shutil
import sys

import pandas as pd

sys.path.insert(0, "/proj/c.zihao/work3/function")
from perpas_adapted import run_perpas_adapted

SAMPLES = ["GSM2560248", "GSM2560249"]
EXPR_DIR = "/proj/c.zihao/work3/00data/singleCell/GSE96583/out"
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
LINKS = "/proj/c.zihao/work3/00data/string/links.csv"
OUT_BASE = "/proj/c.zihao/work3/06singleCell/GSE96583/bench/perpas_adapted"


def reset_out_dir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


reset_out_dir(OUT_BASE)
corum = pd.read_csv(CORUM)
links = pd.read_csv(LINKS)
links = links.loc[:, ~links.columns.str.startswith("Unnamed")]

for sample in SAMPLES:
    expr_path = os.path.join(EXPR_DIR, sample, "logNorm.txt")
    out_dir = os.path.join(OUT_BASE, sample)
    os.makedirs(out_dir, exist_ok=True)
    print(f"[{sample}] perpas_adapted ...", flush=True)
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    print(f"  expr dim: {expr.shape}", flush=True)
    scores = run_perpas_adapted(expr, corum, links).round(5)
    scores.to_csv(os.path.join(out_dir, "perpas_adapted.csv"))
    print(f"  Saved perpas_adapted: {scores.shape}", flush=True)

print("Done.")
