#!/usr/bin/env python3
"""Degree-preserving-PPI control, swept over K and over rewiring seeds.

The matched control for 01_netComplex.py: same K grid, same alpha.  Repeating
the rewiring under many seeds turns each K into a *null distribution* rather
than a single point, so the real network's advantage can be read as its
position inside that distribution instead of a mean +/- cross-cancer SD.

Two structural points:

* All CPTAC cohorts share one filtered gene list, so a seed defines the same
  rewired graph in every cohort.  The graph is therefore built once per seed,
  outside the cancer loop, and the shared gene list is asserted.  This is also
  what makes pooling the null across cohorts within a seed legitimate.
* Only ``complex_score`` is written.  ``node_score`` would cost ~100 GB across
  20 seeds and nothing downstream reads it; ``coverage`` is identical for
  every seed and every K.  Both are still written by 01_netComplex.py for the
  real network.

Run one seed per Slurm array task (``sbatch run_null.sh``), or every seed
in one process when SLURM_ARRAY_TASK_ID is unset.
"""

import os
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/proj/c.zihao/work3/04propagation")
import netComplex


CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
K_VALUES = [0, 1, 2, 3, 5, 8, 12, 20, 50]
N_SEEDS = 20
EXPR_DIR = Path("/proj/c.zihao/work3/00data/cptacT/exprset")
CORUM = Path("/proj/c.zihao/work3/00data/corum/complexData.csv")
LINKS = "/proj/c.zihao/work3/00data/string/links.csv"
OUT_ROOT = Path("/proj/c.zihao/work3/04propagation/ablation/degRandPPI")
SWAPS_PER_EDGE = 5
RWR_ALPHA = 0.3

_task = os.environ.get("SLURM_ARRAY_TASK_ID")
SEEDS = [int(_task)] if _task else list(range(1, N_SEEDS + 1))


links = netComplex.read_links(LINKS)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))

reference = pd.read_csv(EXPR_DIR / f"{CANCERS[0]}_exprSet_filtered.csv", index_col=0)
nodes, _, _ = netComplex.prepare_network(links, reference)
node_pos = {gene: i for i, gene in enumerate(nodes)}
real_edges = links[links.protein1.isin(nodes) & links.protein2.isin(nodes)]
real_edges = real_edges[real_edges.protein1 != real_edges.protein2].drop_duplicates()
edge_pairs = [(node_pos[p1], node_pos[p2])
              for p1, p2 in real_edges[["protein1", "protein2"]].itertuples(index=False)]

for seed in SEEDS:
    seed_root = OUT_ROOT / f"seed_{seed:03d}"
    if seed_root.exists():
        shutil.rmtree(seed_root)

    smoothing, degree = netComplex.degree_preserving_smoothing(
        len(nodes), edge_pairs, seed=seed, swaps_per_edge=SWAPS_PER_EDGE
    )

    for cancer in CANCERS:
        expression = pd.read_csv(EXPR_DIR / f"{cancer}_exprSet_filtered.csv", index_col=0)
        if not nodes.equals(pd.Index(sorted(expression.index))):
            raise ValueError(f"[{cancer}] gene list differs from {CANCERS[0]}; "
                             "the rewired network cannot be shared across cohorts.")
        ranks = netComplex.within_sample_rank(expression, nodes)

        for k in K_VALUES:
            node_scores, n_iter, final_delta = netComplex.propagate_k_steps(
                ranks, smoothing, degree, alpha=RWR_ALPHA, n_steps=k
            )
            complex_scores, _ = netComplex.score_complexes(node_scores, complexes)

            out_dir = seed_root / f"K_{k}" / "complex_score" / cancer
            out_dir.mkdir(parents=True, exist_ok=True)
            complex_scores.round(6).to_csv(
                out_dir / f"{cancer}_degRandPPI_complex_score.csv")

            print(f"[seed={seed}][K={k}][{cancer}] nodes={len(nodes)} "
                  f"edges={len(edge_pairs)} swaps_per_edge={SWAPS_PER_EDGE} "
                  f"complexes={complex_scores.shape[0]} steps={n_iter} "
                  f"delta={final_delta:.2e}", flush=True)
