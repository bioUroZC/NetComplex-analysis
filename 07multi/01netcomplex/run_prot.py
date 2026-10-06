#!/usr/bin/env python3
"""Protein-only NetComplex on matched multi-omics CPTAC samples."""

import os
import sys

import pandas as pd

sys.path.append("/proj/c.zihao/work3/07multi/function")
sys.path.append("/proj/c.zihao/work3/function")
import netComplex


CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]

CORUM_PATH = "/proj/c.zihao/work3/00data/corum/complexData.csv"
LINKS_PATH = "/proj/c.zihao/work3/00data/string/links.csv"
DATA_DIR = "/proj/c.zihao/work3/07multi/00data"
OUT_BASE = "/proj/c.zihao/work3/07multi/netcomplex"
RWR_ALPHA = 0.3

links = netComplex.read_links(LINKS_PATH)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM_PATH))

for cancer in CANCERS:
    prot = pd.read_csv(os.path.join(DATA_DIR, f"{cancer}_protein.csv"), index_col=0)
    common_samples = list(prot.columns)

    output_dir = os.path.join(OUT_BASE, cancer)
    os.makedirs(output_dir, exist_ok=True)
    pd.DataFrame({"sample_id": common_samples}).to_csv(
        os.path.join(output_dir, "matched_samples.csv"),
        index=False,
    )

    print(f"[{cancer}] matched={len(common_samples)} protein-only scoring ...", flush=True)
    nodes, smoothing, degree = netComplex.prepare_network(links, prot)
    ranks = netComplex.within_sample_rank(prot, nodes)
    node_scores, n_iter, final_delta = netComplex.propagate_with_restart(
        ranks, smoothing, degree, alpha=RWR_ALPHA
    )
    complex_scores, _ = netComplex.score_complexes(node_scores, complexes)
    complex_scores.round(6).to_csv(os.path.join(output_dir, "NetComplex_prot.csv"))
    print(
        f"[{cancer}] done complexes={complex_scores.shape[0]} "
        f"iter={n_iter} delta={final_delta:.2e}",
        flush=True,
    )
