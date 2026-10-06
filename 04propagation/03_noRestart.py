#!/usr/bin/env python3
"""No-restart control: pure diffusion on the real network, swept over K.

Identical to 01_netComplex.py except that the restart coefficient is zero, so
the update is ``X <- S @ X`` with no term pulling the original ranks back in.
This is the control for the *necessity of the restart*, which the K sweep at
alpha = 0.3 cannot address on its own.

Without a restart the iteration has no fixed point that retains the original
signal: within each connected component it converges to a degree-proportional
constant times a per-sample scalar, so every gene ends up with the same
across-sample profile and the complex-score matrix collapses towards rank one.
The expected reading is therefore a curve that *falls* with K, in contrast to
the alpha = 0.3 curve which rises and then plateaus.

Genes the PPI does not cover stay isolated and keep their rank here too, so
the collapse is confined to the connected part of the network.
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
OUT_ROOT = Path("/proj/c.zihao/work3/04propagation/ablation/noRestart")
RWR_ALPHA = 0.0


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
        complex_scores, _ = netComplex.score_complexes(node_scores, complexes)

        out_dir = OUT_ROOT / f"K_{k}" / "complex_score" / cancer
        out_dir.mkdir(parents=True, exist_ok=True)
        complex_scores.round(6).to_csv(
            out_dir / f"{cancer}_noRestart_complex_score.csv")

        print(f"[K={k}][{cancer}] nodes={len(nodes)} connected={int((degree > 0).sum())} "
              f"complexes={complex_scores.shape[0]} steps={n_iter} "
              f"delta={final_delta:.2e} alpha={RWR_ALPHA}", flush=True)
