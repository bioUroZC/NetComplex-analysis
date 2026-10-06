#!/usr/bin/env python3
"""Three-module NetComplex, swept over a fixed number of propagation steps.

Same model as the production script; the only change is that propagation
stops after exactly K steps instead of running to convergence.  K = 0 is the
untouched rank matrix (the noRWR baseline) and K = 50 is already the
converged solution at alpha = 0.3, so the whole curve is bracketed.
"""

import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/proj/c.zihao/work3/04propagation")
import netComplex


CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
K_VALUES = [0, 1, 2, 3, 5, 8, 12, 20, 50]
EXPR_DIR = Path("/proj/c.zihao/work3/00data/cptacT/exprset")
CORUM = Path("/proj/c.zihao/work3/00data/corum/complexData.csv")
LINKS = "/proj/c.zihao/work3/00data/string/links.csv"
OUT_ROOT = Path("/proj/c.zihao/work3/04propagation/netcomplex")
RWR_ALPHA = 0.3


if OUT_ROOT.exists():
    shutil.rmtree(OUT_ROOT)

links = netComplex.read_links(LINKS)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))
for cancer in CANCERS:
    expression = pd.read_csv(EXPR_DIR / f"{cancer}_exprSet_filtered.csv", index_col=0)
    nodes, smoothing, degree = netComplex.prepare_network(links, expression)
    ranks = netComplex.within_sample_rank(expression, nodes)

    for k in K_VALUES:
        node_scores, n_iter, final_delta = netComplex.propagate_k_steps(
            ranks, smoothing, degree, alpha=RWR_ALPHA, n_steps=k
        )
        complex_scores, coverage = netComplex.score_complexes(node_scores, complexes)

        for kind, matrix, filename in (
            ("node_score", node_scores, f"{cancer}_netcomplex_node_score.csv"),
            ("complex_score", complex_scores, f"{cancer}_netcomplex_complex_score.csv"),
            ("coverage", coverage, f"{cancer}_netcomplex_coverage.csv"),
        ):
            out_dir = OUT_ROOT / f"K_{k}" / kind / cancer
            out_dir.mkdir(parents=True, exist_ok=True)
            matrix.round(6).to_csv(out_dir / filename)

        print(f"[K={k}][{cancer}] nodes={len(nodes)} complexes={complex_scores.shape[0]} "
              f"steps={n_iter} delta={final_delta:.2e} alpha={RWR_ALPHA}", flush=True)
