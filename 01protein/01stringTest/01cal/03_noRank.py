#!/usr/bin/env python3
"""No-rank ablation: raw expression -> RWR -> complex mean.

This removes only NetComplex module 1 (within-sample rank transformation).
The PPI, restart propagation, and complex aggregation are identical to the
main NetComplex model.
"""

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
OUT_ROOT = Path("/proj/c.zihao/work3/01protein/01stringTest/ablation/noRank")
RWR_ALPHA = 0.3


if OUT_ROOT.exists():
    shutil.rmtree(OUT_ROOT)

links = netComplex.read_links(LINKS)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))
for cancer in CANCERS:
    expression = pd.read_csv(EXPR_DIR / f"{cancer}_exprSet_filtered.csv", index_col=0)
    nodes, smoothing, degree = netComplex.prepare_network(links, expression)
    node_values = expression.loc[nodes].astype(float)
    if node_values.isna().any().any():
        raise ValueError(f"[{cancer}] expression has missing values at PPI nodes")

    node_scores, n_iter, final_delta = netComplex.propagate_with_restart(
        node_values, smoothing, degree, alpha=RWR_ALPHA
    )
    complex_scores, coverage = netComplex.score_complexes(node_scores, complexes)

    for kind, matrix, filename in (
        ("node_score", node_scores, f"{cancer}_noRank_node_score.csv"),
        ("complex_score", complex_scores, f"{cancer}_noRank_complex_score.csv"),
        ("coverage", coverage, f"{cancer}_noRank_coverage.csv"),
    ):
        out_dir = OUT_ROOT / kind / cancer
        out_dir.mkdir(parents=True, exist_ok=True)
        matrix.round(6).to_csv(out_dir / filename)

    print(f"[{cancer}] nodes={len(nodes)} complexes={complex_scores.shape[0]} "
          f"RWR iterations={n_iter} delta={final_delta:.2e} raw_expression=True "
          f"alpha={RWR_ALPHA}", flush=True)
