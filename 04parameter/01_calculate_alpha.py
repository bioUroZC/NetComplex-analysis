#!/usr/bin/env python3
"""Calculate STRING NetComplex scores across RWR restart alpha values.

Uses the same full STRING network and CPTAC cohorts as the primary
01protein/01stringTest proteome-concordance benchmark.
The only model parameter swept here is alpha; PPI, rank, and complex
aggregation stay fixed.
"""

import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.append("/proj/c.zihao/work3/function")
import netComplex


CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
ALPHAS = (0.1, 0.3, 0.5, 0.7, 0.9)
EXPR_DIR = Path("/proj/c.zihao/work3/00data/cptacT/exprset")
CORUM = Path("/proj/c.zihao/work3/00data/corum/complexData.csv")
LINKS = Path("/proj/c.zihao/work3/00data/string/links.csv")
OUT_ROOT = Path("/proj/c.zihao/work3/04parameter/results/complex_score")


if OUT_ROOT.exists():
    shutil.rmtree(OUT_ROOT)

links = netComplex.read_links(LINKS)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))
for cancer in CANCERS:
    expression = pd.read_csv(EXPR_DIR / f"{cancer}_exprSet_filtered.csv", index_col=0)
    nodes, smoothing, degree = netComplex.prepare_network(links, expression)
    ranks = netComplex.within_sample_rank(expression, nodes)

    for alpha in ALPHAS:
        node_scores, n_iter, final_delta = netComplex.propagate_with_restart(
            ranks, smoothing, degree, alpha=alpha
        )
        complex_scores, _ = netComplex.score_complexes(node_scores, complexes)
        output_dir = OUT_ROOT / f"alpha_{alpha:.1f}" / cancer
        output_dir.mkdir(parents=True, exist_ok=True)
        output = output_dir / f"{cancer}_complex_score.csv"
        complex_scores.round(6).to_csv(output)
        print(f"[alpha={alpha:.1f}][{cancer}] nodes={len(nodes)} "
              f"complexes={complex_scores.shape[0]} "
              f"iterations={n_iter} delta={final_delta:.2e}", flush=True)
