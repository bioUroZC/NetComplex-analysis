"""
NetComplex SC variant: KNN-PCA expression smoothing, then classic NetComplex.

Default: k=20, lam=0.5, n_pcs=50. Saves smooth.csv per sample.
"""

import os
import shutil
import sys

import pandas as pd

sys.path.insert(0, "/proj/c.zihao/work3/06singleCell/function")
sys.path.append("/proj/c.zihao/work3/function")
import netComplex
from netComplex_smooth import score_cells_smooth

SAMPLES = ["GSM7079694", "GSM7079695", "GSM7079696", "GSM7079702", "GSM7079703", "GSM7079705", "GSM7079710", "GSM7079711", "GSM7079713"]
EXPR_DIR = "/proj/c.zihao/work3/00data/singleCell/GSE226572/out"
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
LINKS = "/proj/c.zihao/work3/00data/string/links.csv"
OUT_BASE = "/proj/c.zihao/work3/06singleCell/GSE226572/netcomplex/smooth"
RWR_ALPHA = 0.3
N_NEIGHBORS = 20
LAM = 0.5
N_PCS = 50


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
    print(f"[{sample}] smooth ...", flush=True)
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    print(f"  expr dim: {expr.shape}", flush=True)
    scores, _, _, meta = score_cells_smooth(
        links,
        expr,
        complexes,
        alpha=RWR_ALPHA,
        n_neighbors=N_NEIGHBORS,
        lam=LAM,
        n_pcs=N_PCS,
    )
    scores.round(5).to_csv(os.path.join(out_dir, "smooth.csv"))
    print(
        f"  Saved smooth: {scores.shape} "
        f"nodes={meta['n_nodes']} iter={meta['n_iter']} "
        f"delta={meta['final_delta']:.2e} "
        f"k={N_NEIGHBORS} lam={LAM} pcs={N_PCS}",
        flush=True,
    )

print("Done.")
