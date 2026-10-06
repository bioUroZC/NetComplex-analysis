"""Multilayer NetComplex: RNA and protein layers on a fixed PPI.

Aligned with the single-omics three-module core in ``netComplex.py``:

1. Within-sample ranks on each omics layer.
2. Build a two-layer graph with **fixed** unweighted intra-layer PPI
   edges and constant inter-layer edges; run the same restart RWR.
3. Collapse layer node values to gene scores, then complex means.

Dual-omics information enters through per-layer restarts and the
inter-layer coupling. Intra-layer edges are not expression-modulated.
There is no coherence term.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

import sys
sys.path.append("/proj/c.zihao/work3/function")
import netComplex


def _validate_aligned(rna: pd.DataFrame, protein: pd.DataFrame) -> None:
    if not rna.index.equals(protein.index):
        raise ValueError(
            "RNA and protein matrices must have identical gene order."
        )
    if not rna.columns.equals(protein.columns):
        raise ValueError(
            "RNA and protein matrices must have identical sample order."
        )


def prepare_multilayer_network(
    links: pd.DataFrame,
    expression: pd.DataFrame,
    interlayer_weight: float = 1.0,
):
    """Build a row-normalised 2N multilayer transition on a fixed PPI.

    Node layout: ``[RNA::gene_0 .. RNA::gene_{N-1},
                    PROT::gene_0 .. PROT::gene_{N-1}]``.
    Intra-layer PPI edges have weight 1; same-gene inter-layer edges
    have weight ``interlayer_weight``.
    """
    if interlayer_weight <= 0:
        raise ValueError("interlayer_weight must be positive.")

    nodes, _, single_degree = netComplex.prepare_network(links, expression)
    n = len(nodes)
    node_pos = {gene: i for i, gene in enumerate(nodes)}

    link_genes = set(nodes)
    edges = links[
        links.protein1.isin(link_genes) & links.protein2.isin(link_genes)
    ].copy()
    edges = edges[edges.protein1 != edges.protein2].drop_duplicates()

    row, col, data = [], [], []

    def _add_undirected(i: int, j: int, weight: float) -> None:
        row.extend([i, j])
        col.extend([j, i])
        data.extend([weight, weight])

    for u, v in zip(edges.protein1.to_numpy(), edges.protein2.to_numpy()):
        i = node_pos[u]
        j = node_pos[v]
        _add_undirected(i, j, 1.0)                     
        _add_undirected(i + n, j + n, 1.0)                 

    for i in range(n):
        _add_undirected(i, i + n, interlayer_weight)

    adjacency = csr_matrix((data, (row, col)), shape=(2 * n, 2 * n))
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    inv_degree = np.divide(1.0, degree, out=np.zeros_like(degree), where=degree > 0)
    smoothing = adjacency.multiply(inv_degree[:, None]).tocsr()
    return nodes, smoothing, degree, single_degree


def collapse_layer_scores(
    multilayer_scores: pd.DataFrame,
    nodes: pd.Index,
    weight_rna: float = 0.5,
) -> pd.DataFrame:
    """Fold RNA/protein layer node values into gene-level scores."""
    if not 0.0 <= weight_rna <= 1.0:
        raise ValueError("weight_rna must lie in [0, 1].")
    n = len(nodes)
    if multilayer_scores.shape[0] != 2 * n:
        raise ValueError(
            f"Expected 2N={2 * n} multilayer rows, got {multilayer_scores.shape[0]}."
        )
    rna_part = multilayer_scores.iloc[:n].to_numpy(dtype=float)
    prot_part = multilayer_scores.iloc[n:].to_numpy(dtype=float)
    combined = weight_rna * rna_part + (1.0 - weight_rna) * prot_part
    return pd.DataFrame(combined, index=nodes, columns=multilayer_scores.columns)


def score_multilayer(
    rna_expression: pd.DataFrame,
    protein_expression: pd.DataFrame,
    links: pd.DataFrame,
    complexes: dict[str, list[str]],
    alpha: float = netComplex.ALPHA,
    weight_rna: float = 0.5,
    interlayer_weight: float = 1.0,
):
    """Score complexes with multilayer NetComplex on a fixed PPI.

    Parameters
    ----------
    rna_expression, protein_expression
        Pre-aligned gene x sample matrices.
    links
        PPI edge table as returned by ``netComplex.read_links``.
    complexes
        Parsed CORUM dict from ``netComplex.parse_complexes``.
    alpha
        RWR restart probability (same meaning as single-omics).
    weight_rna
        Weight on the RNA layer when collapsing to gene scores.
    interlayer_weight
        Constant weight of same-gene RNA–protein edges.

    Returns
    -------
    complex_scores : pd.DataFrame
        Complex x sample score matrix.
    coverage : pd.DataFrame
        Complex x sample coverage diagnostic.
    gene_scores : pd.DataFrame
        Collapsed gene x sample node values.
    n_iter : int
        RWR iterations to convergence.
    final_delta : float
        Final max absolute update.
    """
    _validate_aligned(rna_expression, protein_expression)

    nodes, smoothing, degree, _ = prepare_multilayer_network(
        links, rna_expression, interlayer_weight=interlayer_weight
    )
    rank_rna = netComplex.within_sample_rank(rna_expression, nodes)
    rank_prot = netComplex.within_sample_rank(protein_expression, nodes)

    restart = pd.concat([rank_rna, rank_prot], axis=0)
    restart.index = (
        [f"RNA::{g}" for g in nodes] + [f"PROT::{g}" for g in nodes]
    )

    multilayer_scores, n_iter, final_delta = netComplex.propagate_with_restart(
        restart, smoothing, degree, alpha=alpha
    )
    gene_scores = collapse_layer_scores(
        multilayer_scores, nodes, weight_rna=weight_rna
    )
    complex_scores, coverage = netComplex.score_complexes(gene_scores, complexes)
    return complex_scores, coverage, gene_scores, n_iter, final_delta
