"""
Classic NetComplex on single-cell data (cell = sample column).

Three modules: within-cell rank -> restart RWR on full STRING -> complex mean.
No SC-specific expression smoothing. Saves netComplex.csv per sample.
"""

import os
import shutil
import sys

import pandas as pd

sys.path.insert(0, "/proj/c.zihao/work3/06singleCell/function")
sys.path.append("/proj/c.zihao/work3/function")
import netComplex

SAMPLES = ["GSM7079694", "GSM7079695", "GSM7079696", "GSM7079702", "GSM7079703", "GSM7079705", "GSM7079710", "GSM7079711", "GSM7079713"]
EXPR_DIR = "/proj/c.zihao/work3/00data/singleCell/GSE226572/out"
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
LINKS = "/proj/c.zihao/work3/00data/string/links.csv"
OUT_BASE = "/proj/c.zihao/work3/06singleCell/GSE226572/netcomplex/bulk"
RWR_ALPHA = 0.3


def reset_out_dir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


reset_out_dir(OUT_BASE)
links = netComplex.read_links(LINKS)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))

for sample in SAMPLES:
    expr_path = os.path.join(EXPR_DIR, sample, "logNorm.txt")
    out_dir = os.path.join(OUT_BASE, sample)
    os.makedirs(out_dir, exist_ok=True)
    print(f"[{sample}] netComplex (classic) ...", flush=True)
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    print(f"  expr dim: {expr.shape}", flush=True)
    nodes, smoothing, degree = netComplex.prepare_network(links, expr)
    ranks = netComplex.within_sample_rank(expr, nodes)
    node_scores, n_iter, final_delta = netComplex.propagate_with_restart(
        ranks, smoothing, degree, alpha=RWR_ALPHA
    )
    scores, _ = netComplex.score_complexes(node_scores, complexes)
    scores.round(5).to_csv(os.path.join(out_dir, "netComplex.csv"))
    print(
        f"  Saved netComplex: {scores.shape} "
        f"nodes={len(nodes)} iter={n_iter} "
        f"delta={final_delta:.2e}",
        flush=True,
    )

print("Done.")
