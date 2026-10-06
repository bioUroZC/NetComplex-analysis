#!/usr/bin/env python3
"""netComplex for GSE13507."""

import sys
from pathlib import Path

import pandas as pd

sys.path.append("/proj/c.zihao/work3/function")
import netComplex

EXPR = Path("/proj/c.zihao/work3/03Survival/BLCA/GSE13507/data/exprSet.csv")
CORUM = Path("/proj/c.zihao/work3/00data/corum/complexData.csv")
LINKS = Path("/proj/c.zihao/work3/00data/string/links.csv")
RWR_ALPHA = 0.3

links = netComplex.read_links(LINKS)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))
expression = pd.read_csv(EXPR, index_col=0)
nodes, smoothing, degree = netComplex.prepare_network(links, expression)
ranks = netComplex.within_sample_rank(expression, nodes)
node_scores, n_iter, final_delta = netComplex.propagate_with_restart(
    ranks, smoothing, degree, alpha=RWR_ALPHA
)
complex_scores, coverage = netComplex.score_complexes(node_scores, complexes)

out_dir = Path("/proj/c.zihao/work3/03Survival/BLCA/GSE13507/bench/netcomplex")
out_dir.mkdir(parents=True, exist_ok=True)
complex_scores.round(6).to_csv(out_dir / "GSE13507_netcomplex_complex_score.csv")

print(f"[GSE13507] nodes={len(nodes)} complexes={complex_scores.shape[0]} "
      f"RWR iterations={n_iter} delta={final_delta:.2e} alpha={RWR_ALPHA}", flush=True)
