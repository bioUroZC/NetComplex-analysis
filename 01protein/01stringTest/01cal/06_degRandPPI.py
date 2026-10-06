#!/usr/bin/env python3
"""Degree-preserving-PPI ablation: rank -> rewired-PPI RWR -> complex mean."""

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
OUT_ROOT = Path("/proj/c.zihao/work3/01protein/01stringTest/ablation/degRandPPI")
SEED = 1
SWAPS_PER_EDGE = 5
RWR_ALPHA = 0.3


if OUT_ROOT.exists():
    shutil.rmtree(OUT_ROOT)

links = netComplex.read_links(LINKS)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))
for cancer_idx, cancer in enumerate(CANCERS):
    expression = pd.read_csv(EXPR_DIR / f"{cancer}_exprSet_filtered.csv", index_col=0)
    nodes, _, _ = netComplex.prepare_network(links, expression)
    node_pos = {gene: i for i, gene in enumerate(nodes)}
    real_edges = links[links.protein1.isin(nodes) & links.protein2.isin(nodes)]
    real_edges = real_edges[real_edges.protein1 != real_edges.protein2].drop_duplicates()
    edge_pairs = [(node_pos[p1], node_pos[p2])
                  for p1, p2 in real_edges[["protein1", "protein2"]].itertuples(index=False)]
    smoothing, degree = netComplex.degree_preserving_smoothing(
        len(nodes), edge_pairs, seed=SEED + cancer_idx, swaps_per_edge=SWAPS_PER_EDGE
    )
    ranks = netComplex.within_sample_rank(expression, nodes)
    node_scores, n_iter, final_delta = netComplex.propagate_with_restart(
        ranks, smoothing, degree, alpha=RWR_ALPHA
    )
    complex_scores, coverage = netComplex.score_complexes(node_scores, complexes)

    for kind, matrix, filename in (
        ("node_score", node_scores, f"{cancer}_degRandPPI_node_score.csv"),
        ("complex_score", complex_scores, f"{cancer}_degRandPPI_complex_score.csv"),
        ("coverage", coverage, f"{cancer}_degRandPPI_coverage.csv"),
    ):
        out_dir = OUT_ROOT / kind / cancer
        out_dir.mkdir(parents=True, exist_ok=True)
        matrix.round(6).to_csv(out_dir / filename)

    print(f"[{cancer}] nodes={len(nodes)} edges={len(edge_pairs)} "
          f"swaps_per_edge={SWAPS_PER_EDGE} complexes={complex_scores.shape[0]} "
          f"RWR iterations={n_iter} delta={final_delta:.2e} seed={SEED + cancer_idx}", flush=True)
