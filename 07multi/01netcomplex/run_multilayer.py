#!/usr/bin/env python3
"""NetComplex_multilayer: fixed dual-layer PPI + inter-layer edges."""

import os
import sys

import pandas as pd

sys.path.append("/proj/c.zihao/work3/07multi/function")
sys.path.append("/proj/c.zihao/work3/function")
import netComplex
import netComplex_multilayer


CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]

CORUM_PATH = "/proj/c.zihao/work3/00data/corum/complexData.csv"
LINKS_PATH = "/proj/c.zihao/work3/00data/string/links.csv"
DATA_DIR = "/proj/c.zihao/work3/07multi/00data"
OUT_BASE = "/proj/c.zihao/work3/07multi/netcomplex"
RWR_ALPHA = 0.3
WEIGHT_RNA = 0.5
INTERLAYER_WEIGHT = 1.0

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
        f"[{cancer}] matched={len(common_samples)} NetComplex_multilayer ...",
        flush=True,
    )
    complex_scores, _, _, n_iter, final_delta = netComplex_multilayer.score_multilayer(
        rna_expression=rna,
        protein_expression=prot,
        links=links,
        complexes=complexes,
        alpha=RWR_ALPHA,
        weight_rna=WEIGHT_RNA,
        interlayer_weight=INTERLAYER_WEIGHT,
    )
    complex_scores.round(6).to_csv(
        os.path.join(output_dir, "NetComplex_multilayer.csv")
    )
    print(
        f"[{cancer}] done complexes={complex_scores.shape[0]} "
        f"iter={n_iter} delta={final_delta:.2e}",
        flush=True,
    )
