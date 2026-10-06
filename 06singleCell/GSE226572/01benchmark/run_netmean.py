"""
Static network-aware baseline using complex-internal PPI node strength
as weights for the mean of member gene z-scores.
Loops over GSE226572 samples and saves netmean.csv per sample.
"""

import os
import shutil
import sys

import pandas as pd

sys.path.insert(0, "/proj/c.zihao/work3/function")
from netmean import run_netmean

SAMPLES = ["GSM7079694", "GSM7079695", "GSM7079696", "GSM7079702", "GSM7079703", "GSM7079705", "GSM7079710", "GSM7079711", "GSM7079713"]
EXPR_DIR = "/proj/c.zihao/work3/00data/singleCell/GSE226572/out"
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
LINKS = "/proj/c.zihao/work3/00data/string/links.csv"
OUT_BASE = "/proj/c.zihao/work3/06singleCell/GSE226572/bench/netmean"


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
    print(f"[{sample}] netmean ...", flush=True)
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    print(f"  expr dim: {expr.shape}", flush=True)
    scores = run_netmean(expr, corum, links).round(5)
    scores.to_csv(os.path.join(out_dir, "netmean.csv"))
    print(f"  Saved netmean: {scores.shape}", flush=True)

print("Done.")
