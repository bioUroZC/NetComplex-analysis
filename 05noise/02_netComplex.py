#!/usr/bin/env python3
"""Three-module NetComplex: within-sample rank, restart RWR, complex mean.

Same layout as 01protein/01stringTest/01cal/01_netComplex.py, but expression
comes from 05noise/01data/<cancer>/expr_level*_seed*.csv.
"""

import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.append("/proj/c.zihao/work3/function")
import netComplex


CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
EXPR_DIR = Path("/proj/c.zihao/work3/05noise/01data")
CORUM = Path("/proj/c.zihao/work3/00data/corum/complexData.csv")
LINKS = Path("/proj/c.zihao/work3/00data/string/links.csv")
OUT_ROOT = Path("/proj/c.zihao/work3/05noise/results/netcomplex")
RWR_ALPHA = 0.3


if OUT_ROOT.exists():
    shutil.rmtree(OUT_ROOT)

links = netComplex.read_links(LINKS)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))
for cancer in CANCERS:
    expr_files = sorted((EXPR_DIR / cancer).glob("expr_level*_seed*.csv"))
    if not expr_files:
        raise FileNotFoundError(f"No expression files under {EXPR_DIR / cancer}")

    clean = pd.read_csv(EXPR_DIR / cancer / "expr_level0_seed0.csv", index_col=0)
    nodes, smoothing, degree = netComplex.prepare_network(links, clean)

    for expr_path in expr_files:
        tag = expr_path.stem.replace("expr_", "")                     
        expression = pd.read_csv(expr_path, index_col=0)
        ranks = netComplex.within_sample_rank(expression, nodes)
        node_scores, n_iter, final_delta = netComplex.propagate_with_restart(
            ranks, smoothing, degree, alpha=RWR_ALPHA
        )
        complex_scores, coverage = netComplex.score_complexes(node_scores, complexes)

        for kind, matrix, filename in (
            ("node_score", node_scores,
             f"{cancer}_{tag}_netcomplex_node_score.csv"),
            ("complex_score", complex_scores,
             f"{cancer}_{tag}_netcomplex_complex_score.csv"),
            ("coverage", coverage,
             f"{cancer}_{tag}_netcomplex_coverage.csv"),
        ):
            out_dir = OUT_ROOT / kind / cancer
            out_dir.mkdir(parents=True, exist_ok=True)
            matrix.round(6).to_csv(out_dir / filename)

        print(f"[{cancer}/{tag}] nodes={len(nodes)} complexes={complex_scores.shape[0]} "
              f"RWR iterations={n_iter} delta={final_delta:.2e} alpha={RWR_ALPHA}",
              flush=True)
