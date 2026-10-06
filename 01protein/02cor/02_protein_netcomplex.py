#!/usr/bin/env python3
"""Protein-NetComplex: rank + RWR + mean, run on protein abundance instead
of RNA expression.

One of the nine protein-side references matching 01stringTest/01benchmark.
Calls the exact same function/netComplex.py pipeline used for the RNA-side
production benchmark -- prepare_network/within_sample_rank/
propagate_with_restart/score_complexes are generic over "some gene-by-sample
matrix", so passing the stable, zero-missing protein matrix from
01_stable_subset.py in place of RNA needs no code change. Because
prepare_network takes its node set directly from whatever expression matrix
it is given, the PPI network here is automatically restricted to the stable
protein subset -- there is no separate network-filtering step, and no
imputation.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.append("/proj/c.zihao/work3/function")
import netComplex


CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
STABLE_DIR = Path("/proj/c.zihao/work3/01protein/02cor/01stable")
CORUM = Path("/proj/c.zihao/work3/00data/corum/complexData.csv")
LINKS = Path("/proj/c.zihao/work3/00data/string/links.csv")
OUT_ROOT = Path("/proj/c.zihao/work3/01protein/02cor/02scored/netcomplex")
RWR_ALPHA = 0.3

links = netComplex.read_links(LINKS)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))

for cancer in CANCERS:
    protein = pd.read_csv(STABLE_DIR / f"{cancer}_protein_stable.csv", index_col=0)
    nodes, smoothing, degree = netComplex.prepare_network(links, protein)
    ranks = netComplex.within_sample_rank(protein, nodes)
    node_scores, n_iter, final_delta = netComplex.propagate_with_restart(
        ranks, smoothing, degree, alpha=RWR_ALPHA
    )
    complex_scores, coverage = netComplex.score_complexes(node_scores, complexes)

    for kind, matrix, filename in (
        ("node_score", node_scores, f"{cancer}_netcomplex_node_score.csv"),
        ("complex_score", complex_scores, f"{cancer}_netcomplex_complex_score.csv"),
        ("coverage", coverage, f"{cancer}_netcomplex_coverage.csv"),
    ):
        out_dir = OUT_ROOT / kind / cancer
        out_dir.mkdir(parents=True, exist_ok=True)
        matrix.round(6).to_csv(out_dir / filename)

    print(f"[{cancer}] stable_proteins={len(nodes)} complexes={complex_scores.shape[0]} "
          f"RWR iterations={n_iter} delta={final_delta:.2e} alpha={RWR_ALPHA}", flush=True)

print("Done.")
