#!/usr/bin/env python3
"""No-RWR baseline: within-sample rank followed by direct complex averaging.

Same layout as 01protein/01stringTest/01cal/02_noRWR.py, but expression
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
OUT_ROOT = Path("/proj/c.zihao/work3/05noise/results/noRWR")


if OUT_ROOT.exists():
    shutil.rmtree(OUT_ROOT)

links = netComplex.read_links(LINKS)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))
for cancer in CANCERS:
    expr_files = sorted((EXPR_DIR / cancer).glob("expr_level*_seed*.csv"))
    if not expr_files:
        raise FileNotFoundError(f"No expression files under {EXPR_DIR / cancer}")

    clean = pd.read_csv(EXPR_DIR / cancer / "expr_level0_seed0.csv", index_col=0)
    nodes, _, _ = netComplex.prepare_network(links, clean)

    for expr_path in expr_files:
        tag = expr_path.stem.replace("expr_", "")                     
        expression = pd.read_csv(expr_path, index_col=0)
        ranks = netComplex.within_sample_rank(expression, nodes)
        complex_scores, _ = netComplex.score_complexes(ranks, complexes)

        output_dir = OUT_ROOT / cancer
        output_dir.mkdir(parents=True, exist_ok=True)
        complex_scores.round(6).to_csv(
            output_dir / f"{cancer}_{tag}_noRWR_complex_score.csv"
        )
        print(f"[{cancer}/{tag}] nodes={len(nodes)} "
              f"complexes={complex_scores.shape[0]}", flush=True)
