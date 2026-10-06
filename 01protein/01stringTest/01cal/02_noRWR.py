#!/usr/bin/env python3
"""No-RWR baseline: within-sample rank followed by direct complex averaging."""

import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.append("/proj/c.zihao/work3/function")
import netComplex


CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
EXPR_DIR = Path("/proj/c.zihao/work3/00data/cptacT/exprset")
CORUM = Path("/proj/c.zihao/work3/00data/corum/complexData.csv")
LINKS = Path("/proj/c.zihao/work3/00data/string/links.csv")
OUT_ROOT = Path("/proj/c.zihao/work3/01protein/01stringTest/ablation/noRWR")


if OUT_ROOT.exists():
    shutil.rmtree(OUT_ROOT)

links = netComplex.read_links(LINKS)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))
for cancer in CANCERS:
    expression = pd.read_csv(EXPR_DIR / f"{cancer}_exprSet_filtered.csv", index_col=0)
    nodes, _, _ = netComplex.prepare_network(links, expression)
    ranks = netComplex.within_sample_rank(expression, nodes)
    complex_scores, _ = netComplex.score_complexes(ranks, complexes)

    output_dir = OUT_ROOT / cancer
    output_dir.mkdir(parents=True, exist_ok=True)
    complex_scores.round(6).to_csv(output_dir / f"{cancer}_noRWR_complex_score.csv")
    print(f"[{cancer}] nodes={len(nodes)} complexes={complex_scores.shape[0]}", flush=True)
