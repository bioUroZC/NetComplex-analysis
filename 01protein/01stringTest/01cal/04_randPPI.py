#!/usr/bin/env python3
"""Random-PPI ablation: rank -> RWR on a node/edge-matched random PPI -> mean."""

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.append("/proj/c.zihao/work3/function")
import netComplex


CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
EXPR_DIR = Path("/proj/c.zihao/work3/00data/cptacT/exprset")
CORUM = Path("/proj/c.zihao/work3/00data/corum/complexData.csv")
LINKS = Path("/proj/c.zihao/work3/00data/string/links.csv")
OUT_ROOT = Path("/proj/c.zihao/work3/01protein/01stringTest/ablation/randPPI")
SEED = 1
RWR_ALPHA = 0.3


if OUT_ROOT.exists():
    shutil.rmtree(OUT_ROOT)

links = netComplex.read_links(LINKS)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))
for cancer_idx, cancer in enumerate(CANCERS):
    expression = pd.read_csv(EXPR_DIR / f"{cancer}_exprSet_filtered.csv", index_col=0)
    nodes, _, real_degree = netComplex.prepare_network(links, expression)
    connected = np.flatnonzero(real_degree > 0)
    real_edges = links[links.protein1.isin(nodes) & links.protein2.isin(nodes)]
    real_edges = real_edges[real_edges.protein1 != real_edges.protein2].drop_duplicates()
    smoothing, degree = netComplex.random_network_smoothing(
        len(nodes), len(real_edges), seed=SEED + cancer_idx, active=connected
    )
    ranks = netComplex.within_sample_rank(expression, nodes)
    node_scores, n_iter, final_delta = netComplex.propagate_with_restart(
        ranks, smoothing, degree, alpha=RWR_ALPHA
    )
    complex_scores, coverage = netComplex.score_complexes(node_scores, complexes)

    for kind, matrix, filename in (
        ("node_score", node_scores, f"{cancer}_randPPI_node_score.csv"),
        ("complex_score", complex_scores, f"{cancer}_randPPI_complex_score.csv"),
        ("coverage", coverage, f"{cancer}_randPPI_coverage.csv"),
    ):
        out_dir = OUT_ROOT / kind / cancer
        out_dir.mkdir(parents=True, exist_ok=True)
        matrix.round(6).to_csv(out_dir / filename)

    print(f"[{cancer}] nodes={len(nodes)} connected={len(connected)} "
          f"random_edges={len(real_edges)} "
          f"complexes={complex_scores.shape[0]} RWR iterations={n_iter} "
          f"delta={final_delta:.2e} seed={SEED + cancer_idx}", flush=True)
