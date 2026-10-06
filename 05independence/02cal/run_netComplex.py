#!/usr/bin/env python3
"""Three-module NetComplex on KIRC cohort-size subsets."""

import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.append("/proj/c.zihao/work3/function")
import netComplex


LABELS = ["frac0.20", "frac0.40", "frac0.60", "frac0.80", "frac1.00"]
EXPR_DIR = Path("/proj/c.zihao/work3/05independence/01data")
CORUM = Path("/proj/c.zihao/work3/00data/corum/complexData.csv")
LINKS = Path("/proj/c.zihao/work3/00data/string/links.csv")
OUT_ROOT = Path("/proj/c.zihao/work3/05independence/results/netcomplex")
RWR_ALPHA = 0.3


if OUT_ROOT.exists():
    shutil.rmtree(OUT_ROOT)

links = netComplex.read_links(LINKS)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))
for label in LABELS:
    expression = pd.read_csv(EXPR_DIR / f"KIRC_{label}_exprSet.csv", index_col=0)
    nodes, smoothing, degree = netComplex.prepare_network(links, expression)
    ranks = netComplex.within_sample_rank(expression, nodes)
    node_scores, n_iter, final_delta = netComplex.propagate_with_restart(
        ranks, smoothing, degree, alpha=RWR_ALPHA
    )
    complex_scores, coverage = netComplex.score_complexes(node_scores, complexes)

    for kind, matrix, filename in (
        ("node_score", node_scores, f"{label}_netcomplex_node_score.csv"),
        ("complex_score", complex_scores, f"{label}_netcomplex_complex_score.csv"),
        ("coverage", coverage, f"{label}_netcomplex_coverage.csv"),
    ):
        out_dir = OUT_ROOT / kind / label
        out_dir.mkdir(parents=True, exist_ok=True)
        matrix.round(6).to_csv(out_dir / filename)

    print(f"[{label}] nodes={len(nodes)} complexes={complex_scores.shape[0]} "
          f"samples={expression.shape[1]} RWR iterations={n_iter} "
          f"delta={final_delta:.2e} alpha={RWR_ALPHA}", flush=True)
