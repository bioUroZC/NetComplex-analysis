#!/usr/bin/env python3
"""noRWR / degRandPPI ablations for both single-cell NetComplex branches.

Four variants:

noRWR_bulk       within-cell rank -> complex mean                    (no RWR, no KNN)
degRandPPI_bulk     rank -> RWR on a degree-preserving rewired PPI -> mean
noRWR_smooth     KNN-PCA smoothing -> rank -> complex mean           (no RWR)
degRandPPI_smooth   KNN smoothing -> rank -> RWR on a rewired PPI -> mean

Together with netcomplex/{bulk,smooth} these separate the three contributions
the single-cell model mixes:

    RWR              netComplex        - noRWR
    real topology    netComplex        - degRandPPI
    KNN smoothing    netComplex_smooth - netComplex

degRandPPI keeps the KNN smoothing untouched and randomises only the PPI, so the
contrast isolates the topology rather than the amount of smoothing.  The control
is a degree-preserving rewiring (every node keeps its exact degree, only which
pairs are joined changes), the same construction the rest of the project uses in
*/01cal/06_degRandPPI.py.  An Erdos-Renyi control would flatten the degree
distribution as well, so its contrast would mix topology with degree
heterogeneity instead of isolating topology.

The KNN step is the expensive one (a (n_cells, k, n_genes) neighbour array), so
it is computed once per sample and shared by both smooth variants.

Network: full STRING (combined_score > 700), matching netcomplex/{bulk,smooth}
and the two PPI benchmarks.  HuRI connected only ~47% of measured CORUM members
at single-cell gene coverage, which left propagation nothing to work with; full
STRING connects 99%.  All PPI-using methods share one network so the benchmark
comparison stays fair.
"""

import gc
import os
import shutil
import sys

import pandas as pd

sys.path.insert(0, "/proj/c.zihao/work3/06singleCell/function")
sys.path.append("/proj/c.zihao/work3/function")
import netComplex
from netComplex_smooth import smooth_expression_knn

SAMPLES = ["GSM7079694", "GSM7079695", "GSM7079696", "GSM7079702", "GSM7079703", "GSM7079705", "GSM7079710", "GSM7079711", "GSM7079713"]
EXPR_DIR = "/proj/c.zihao/work3/00data/singleCell/GSE226572/out"
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
LINKS = "/proj/c.zihao/work3/00data/string/links.csv"
OUT_BASE = "/proj/c.zihao/work3/06singleCell/GSE226572/ablation"
RWR_ALPHA = 0.3
N_NEIGHBORS = 20
LAM = 0.5
N_PCS = 50
SEED = 1
SWAPS_PER_EDGE = 5

VARIANTS = ["noRWR_bulk", "degRandPPI_bulk", "noRWR_smooth", "degRandPPI_smooth"]


def reset_out_dirs():
    for variant in VARIANTS:
        path = os.path.join(OUT_BASE, variant)
        if os.path.isdir(path):
            shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)


def save(scores, variant, sample):
    out_dir = os.path.join(OUT_BASE, variant, sample)
    os.makedirs(out_dir, exist_ok=True)
    scores.round(5).to_csv(os.path.join(out_dir, f"{variant}.csv"))
    print(f"  Saved {variant}: {scores.shape}", flush=True)


def rewired_smoothing_like(links, nodes, seed):
    """Degree-preserving rewiring of the real PPI, as in */01cal/06_degRandPPI.py.

    Genes the PPI does not cover keep degree 0 and stay isolated, exactly as in
    the real network, so the only thing that changes is which pairs are joined.
    """
    node_pos = {gene: i for i, gene in enumerate(nodes)}
    real_edges = links[links.protein1.isin(nodes) & links.protein2.isin(nodes)]
    real_edges = real_edges[real_edges.protein1 != real_edges.protein2].drop_duplicates()
    edge_pairs = [(node_pos[p1], node_pos[p2])
                  for p1, p2 in real_edges[["protein1", "protein2"]].itertuples(index=False)]
    smoothing, rewired_degree = netComplex.degree_preserving_smoothing(
        len(nodes), edge_pairs, seed=seed, swaps_per_edge=SWAPS_PER_EDGE
    )
    return smoothing, rewired_degree, len(edge_pairs)


def score_branch(expression, links, complexes, sample, suffix, seed):
    """Score the noRWR and degRandPPI variants of one branch (bulk or smooth)."""
    nodes, _, _ = netComplex.prepare_network(links, expression)
    ranks = netComplex.within_sample_rank(expression, nodes)

    scores, _ = netComplex.score_complexes(ranks, complexes)
    save(scores, f"noRWR_{suffix}", sample)
    del scores

    rand_smoothing, rand_degree, n_edges = rewired_smoothing_like(links, nodes, seed)
    node_scores, n_iter, final_delta = netComplex.propagate_with_restart(
        ranks, rand_smoothing, rand_degree, alpha=RWR_ALPHA
    )
    scores, _ = netComplex.score_complexes(node_scores, complexes)
    save(scores, f"degRandPPI_{suffix}", sample)
    print(f"    nodes={len(nodes)} rewired_edges={n_edges} "
          f"swaps_per_edge={SWAPS_PER_EDGE} iter={n_iter} "
          f"delta={final_delta:.2e} seed={seed}", flush=True)
    del scores, node_scores, rand_smoothing, ranks
    gc.collect()


reset_out_dirs()
links = netComplex.read_links(LINKS)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))

for sample_idx, sample in enumerate(SAMPLES):
    expr_path = os.path.join(EXPR_DIR, sample, "logNorm.txt")
    print(f"[{sample}] ablations ...", flush=True)
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    print(f"  expr dim: {expr.shape}", flush=True)

    score_branch(expr, links, complexes, sample, "bulk", SEED + sample_idx)

    print(f"  KNN smoothing (k={N_NEIGHBORS} lam={LAM} pcs={N_PCS}) ...", flush=True)
    smoothed = smooth_expression_knn(
        expr, n_neighbors=N_NEIGHBORS, lam=LAM, n_pcs=N_PCS
    )
    del expr
    gc.collect()
    score_branch(smoothed, links, complexes, sample, "smooth", SEED + sample_idx)
    del smoothed
    gc.collect()

print("Done.")
