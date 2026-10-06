#!/usr/bin/env python3
"""How much of each PPI actually touches the genes NetComplex scores.

Since f347e6c every expressed gene is a propagation node; genes the PPI does
not cover enter as isolated nodes and ``propagate_with_restart`` returns their
rank unchanged.  A complex whose members are all isolated therefore receives
exactly the mean of its members' ranks -- the RWR module contributes nothing.

This script quantifies that per network by calling the same
``netComplex.prepare_network`` the production scripts use, so the degrees
reported here are the degrees the scoring run actually saw.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.append("/proj/c.zihao/work3/function")
import netComplex


CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
EXPR_DIR = Path("/proj/c.zihao/work3/00data/cptacT/exprset")
CORUM = Path("/proj/c.zihao/work3/00data/corum/complexData.csv")
NETWORKS = {
    "HuRI": Path("/proj/c.zihao/work3/00data/HuRI/links.csv"),
    "string": Path("/proj/c.zihao/work3/00data/string/links.csv"),
    "biogrid": Path("/proj/c.zihao/work3/00data/biogrid/links.csv"),
}
OUT_DIR = Path("/proj/c.zihao/work3/01protein/count")


complexes = netComplex.parse_complexes(pd.read_csv(CORUM))
corum_genes = pd.Index(sorted({g for members in complexes.values() for g in members}))
links_by_network = {name: netComplex.read_links(path) for name, path in NETWORKS.items()}

rows = []
for cancer in CANCERS:
    expression = pd.read_csv(EXPR_DIR / f"{cancer}_exprSet_filtered.csv", index_col=0)
    for network, links in links_by_network.items():
        nodes, _, degree = netComplex.prepare_network(links, expression)
        degree_by_gene = pd.Series(degree, index=nodes)
        linked = degree_by_gene > 0

        corum_in_expr = nodes.intersection(corum_genes)
        corum_linked = degree_by_gene.loc[corum_in_expr] > 0

        isolated_fracs, n_all_isolated, n_scored = [], 0, 0
        for members in complexes.values():
            present = [g for g in members if g in degree_by_gene.index]
            if not present:
                continue
            n_scored += 1
            frac_isolated = float(np.mean([degree_by_gene[g] == 0 for g in present]))
            isolated_fracs.append(frac_isolated)
            if frac_isolated == 1.0:
                n_all_isolated += 1

        rows.append({
            "network": network,
            "cancer": cancer,
            "n_expressed_genes": len(nodes),
            "n_edges": int(degree.sum() // 2),
            "n_nonisolated": int(linked.sum()),
            "frac_nonisolated": round(float(linked.mean()), 4),
            "mean_degree_all_genes": round(float(degree.mean()), 3),
            "mean_degree_nonisolated": round(float(degree[degree > 0].mean()), 3),
            "n_corum_genes_in_expr": len(corum_in_expr),
            "n_corum_nonisolated": int(corum_linked.sum()),
            "frac_corum_nonisolated": round(float(corum_linked.mean()), 4),
            "n_complexes_scored": n_scored,
            "n_complexes_all_isolated": n_all_isolated,
            "frac_complexes_all_isolated": round(n_all_isolated / n_scored, 4),
            "mean_frac_isolated_members": round(float(np.mean(isolated_fracs)), 4),
        })
        print(f"[{cancer}/{network}] nonisolated={linked.mean():.3f} "
              f"corum_nonisolated={corum_linked.mean():.3f} "
              f"mean_degree={degree.mean():.2f}", flush=True)

table = pd.DataFrame(rows)
OUT_DIR.mkdir(parents=True, exist_ok=True)
table.to_csv(OUT_DIR / "network_degradation.csv", index=False)

summary = table.groupby("network")[[
    "n_expressed_genes", "n_edges", "frac_nonisolated", "mean_degree_all_genes",
    "mean_degree_nonisolated", "frac_corum_nonisolated",
    "frac_complexes_all_isolated", "mean_frac_isolated_members",
]].mean().round(4)
summary.to_csv(OUT_DIR / "network_degradation_summary.csv")
print("\n=== mean over the 10 CPTAC cohorts ===")
print(summary.to_string())
