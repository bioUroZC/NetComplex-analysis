#!/usr/bin/env python3
"""Decompose the selected KIRC NetComplex shift into RWR input sources."""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path("/proj/c.zihao/work3")
sys.path.append(str(ROOT / "function"))
import netComplex              

STUDY = ROOT / "09TCGAstudy"
RESULT_DIR = STUDY / "results"
EXPR_FILE = ROOT / "00data" / "TCGAnt" / "exprset" / "KIRC_exprSet_filtered.csv"
LINKS_FILE = ROOT / "00data" / "string" / "links.csv"
CORUM_FILE = ROOT / "00data" / "corum" / "complexData.csv"
NODE_SCORE_FILE = (ROOT / "02TCGAnt" / "02stringTest" / "netcomplex" / "node_score" /
                   "KIRC" / "KIRC_netcomplex_node_score.csv")
ALPHA = 0.3


def paired_samples(sample_ids: pd.Index):
    tumour = sample_ids.str.endswith("_01A")
    normal = sample_ids.str.endswith("_11A")
    patient = sample_ids.str.replace(r"_(01A|11A)$", "", regex=True)
    tumour_map = pd.Series(sample_ids[tumour].to_numpy(), index=patient[tumour])
    normal_map = pd.Series(sample_ids[normal].to_numpy(), index=patient[normal])
    common = tumour_map.index.intersection(normal_map.index).sort_values()
    if len(common) == 0:
        raise ValueError("No KIRC tumour-normal pairs found")
    return tumour_map.loc[common].tolist(), normal_map.loc[common].tolist()


def paired_delta(values: pd.DataFrame, tumour_ids, normal_ids) -> float:
    return float((values[tumour_ids].to_numpy() - values[normal_ids].to_numpy()).mean())


selected = pd.read_csv(RESULT_DIR / "kirc_selected_case.csv")
complex_name = selected.loc[0, "complex"]

expression = pd.read_csv(EXPR_FILE, index_col=0)
links = netComplex.read_links(str(LINKS_FILE))
nodes, smoothing, _ = netComplex.prepare_network(links, expression)
ranks = netComplex.within_sample_rank(expression, nodes)
node_scores = pd.read_csv(NODE_SCORE_FILE, index_col=0).loc[nodes]
node_scores = node_scores.loc[:, ranks.columns]
tumour_ids, normal_ids = paired_samples(ranks.columns)

complexes = netComplex.parse_complexes(pd.read_csv(CORUM_FILE))
members = [gene for gene in complexes[complex_name] if gene in nodes]

node_position = {gene: i for i, gene in enumerate(nodes)}
member_set = set(members)
n_members = len(members)
contribution = {}
own_term = 0.0

for member in members:
    own_term += paired_delta(ALPHA * ranks.loc[[member]], tumour_ids, normal_ids) / n_members
    i = node_position[member]
    start, end = smoothing.indptr[i], smoothing.indptr[i + 1]
    for edge_pos in range(start, end):
        neighbour = nodes[smoothing.indices[edge_pos]]
        weight = (1.0 - ALPHA) * float(smoothing.data[edge_pos]) / n_members
        values = weight * node_scores.loc[[neighbour]]
        contribution[neighbour] = contribution.get(neighbour, 0.0) + paired_delta(values, tumour_ids, normal_ids)

rows = []
for gene, value in contribution.items():
    rows.append({
        "gene": gene,
        "is_complex_member": gene in member_set,
        "rwr_contribution_to_score_delta": value,
        "raw_expression_delta": paired_delta(expression.loc[[gene]], tumour_ids, normal_ids),
        "rank_delta": paired_delta(ranks.loc[[gene]], tumour_ids, normal_ids),
        "propagated_score_delta": paired_delta(node_scores.loc[[gene]], tumour_ids, normal_ids),
    })
contributions = pd.DataFrame(rows).sort_values("rwr_contribution_to_score_delta", ascending=False)

member_network = float(contributions.loc[contributions.is_complex_member, "rwr_contribution_to_score_delta"].sum())
external_network = float(contributions.loc[~contributions.is_complex_member, "rwr_contribution_to_score_delta"].sum())
netcomplex_delta = paired_delta(node_scores.loc[members], tumour_ids, normal_ids)
components = pd.DataFrame({
    "component": ["Own member ranks", "Member-to-member PPI input", "External PPI-neighbour input", "Total NetComplex shift"],
    "score_delta": [own_term, member_network, external_network, netcomplex_delta],
})

for frame in (contributions, components):
    for column in frame.select_dtypes(include="number").columns:
        frame[column] = [
            f"{value:.5e}" if 0 < abs(value) < 1e-5 else f"{value:.5f}"
            for value in frame[column]
        ]

contributions.to_csv(RESULT_DIR / "kirc_signal_contributions.csv", index=False)
components.to_csv(RESULT_DIR / "kirc_signal_decomposition.csv", index=False)
external = contributions.loc[~contributions.is_complex_member].head(5)
print(f"Selected complex: {complex_name}")
print(components.to_string(index=False, float_format=lambda x: f"{x:.6f}"))
print("Top external PPI contributors:", "; ".join(external.gene.tolist()))
print(f"Decomposition check (components minus total): {own_term + member_network + external_network - netcomplex_delta:.2e}")
