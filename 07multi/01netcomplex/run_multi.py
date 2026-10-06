#!/usr/bin/env python3
"""NetComplex_multi: restart-fusion multi-omics on a fixed PPI."""

import os
import sys

import pandas as pd

sys.path.append("/proj/c.zihao/work3/07multi/function")
sys.path.append("/proj/c.zihao/work3/function")
import netComplex
import netComplex_multi


CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]

CORUM_PATH = "/proj/c.zihao/work3/00data/corum/complexData.csv"
LINKS_PATH = "/proj/c.zihao/work3/00data/string/links.csv"
DATA_DIR = "/proj/c.zihao/work3/07multi/00data"
OUT_BASE = "/proj/c.zihao/work3/07multi/netcomplex"
RWR_ALPHA = 0.3
WEIGHT_BACKBONE = 0.5

links = netComplex.read_links(LINKS_PATH)
complexes = netComplex.parse_complexes(pd.read_csv(CORUM_PATH))

for cancer in CANCERS:
    rna = pd.read_csv(os.path.join(DATA_DIR, f"{cancer}_RNA.csv"), index_col=0)
    prot = pd.read_csv(os.path.join(DATA_DIR, f"{cancer}_protein.csv"), index_col=0)
    common_samples = list(rna.columns)

    output_dir = os.path.join(OUT_BASE, cancer)
    os.makedirs(output_dir, exist_ok=True)
    pd.DataFrame({"sample_id": common_samples}).to_csv(
        os.path.join(output_dir, "matched_samples.csv"),
        index=False,
    )

    print(
        f"[{cancer}] matched={len(common_samples)} NetComplex_multi ...",
        flush=True,
    )
    complex_scores, _, _, n_iter, final_delta = netComplex_multi.score_multiomics(
        backbone_expression=rna,
        auxiliary_expression=prot,
        links=links,
        complexes=complexes,
        alpha=RWR_ALPHA,
        weight_backbone=WEIGHT_BACKBONE,
    )
    complex_scores.round(6).to_csv(os.path.join(output_dir, "NetComplex_multi.csv"))
    print(
        f"[{cancer}] done complexes={complex_scores.shape[0]} "
        f"iter={n_iter} delta={final_delta:.2e}",
        flush=True,
    )
