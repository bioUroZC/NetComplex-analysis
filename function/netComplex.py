"""Rank/RWR NetComplex core.

The model has three modules:
1. Convert each sample's expression to ranks in ``(0, 1)``.
2. Propagate ranks on the fixed experimental PPI with restart.  Expressed
   genes absent from the PPI enter as isolated nodes and keep their rank.
3. Average propagated member values to score each CORUM complex.

Every expressed gene is a node, so which complexes can be scored depends on
the expression matrix alone and not on the PPI in use.
"""

from __future__ import annotations

import numpy as np
import networkx as nx
import pandas as pd
from scipy.sparse import csr_matrix


ALPHA = 0.3
TOLERANCE = 1e-4
MAX_ITER = 1000


def read_links(path: str) -> pd.DataFrame:
    """Read a PPI CSV, accepting the legacy unnamed index column."""
    links = pd.read_csv(path)
    links = links.loc[:, ~links.columns.str.startswith("Unnamed")]
    required = {"protein1", "protein2"}
    missing = required.difference(links.columns)
    if missing:
        raise ValueError(f"PPI links missing required column(s): {sorted(missing)}")
    return links.loc[:, ["protein1", "protein2"]].astype(str)


def prepare_network(links: pd.DataFrame, expression: pd.DataFrame):
    """Build the propagation graph over every expressed gene.

    Nodes are all genes in ``expression``.  Genes that the PPI does not cover
    simply carry no edges, so they are isolated and ``propagate_with_restart``
    leaves them at their original value.  Edges are kept only between genes
    present in both the PPI and the expression matrix.
    """
    link_genes = pd.Index(pd.unique(links[["protein1", "protein2"]].to_numpy().ravel()))
    linked = expression.index.intersection(link_genes)
    if len(linked) == 0:
        raise ValueError("No expression genes overlap the PPI network.")

    nodes = pd.Index(sorted(expression.index))
    node_pos = {gene: i for i, gene in enumerate(nodes)}
    edges = links[links.protein1.isin(linked) & links.protein2.isin(linked)].copy()
    edges = edges[edges.protein1 != edges.protein2].drop_duplicates()
    if edges.empty:
        raise ValueError("No non-self-loop PPI edges remain after expression filtering.")

    row = np.concatenate((edges.protein1.map(node_pos).to_numpy(),
                          edges.protein2.map(node_pos).to_numpy()))
    col = np.concatenate((edges.protein2.map(node_pos).to_numpy(),
                          edges.protein1.map(node_pos).to_numpy()))
    adjacency = csr_matrix((np.ones(len(row)), (row, col)), shape=(len(nodes), len(nodes)))
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    inv_degree = np.divide(1.0, degree, out=np.zeros_like(degree), where=degree > 0)
    smoothing = adjacency.multiply(inv_degree[:, None]).tocsr()
    return nodes, smoothing, degree


def within_sample_rank(expression: pd.DataFrame, nodes: pd.Index) -> pd.DataFrame:
    """Map each sample's expression ranks to ``(0, 1)``.

    ``nodes`` is the full node set from :func:`prepare_network`, i.e. every
    expressed gene, so a gene's percentile is taken against the whole matrix
    and does not shift when the PPI is swapped.
    """
    expr = expression.loc[nodes]
    if expr.isna().any().any():
        raise ValueError("Expression has missing values for one or more nodes.")
    rank = expr.rank(axis=0, method="average")
    return (rank - 0.5) / len(nodes)


def propagate_with_restart(rank_values: pd.DataFrame, smoothing: csr_matrix,
                           degree: np.ndarray, alpha: float = ALPHA):
    """Iterate ``X_new = alpha*R + (1-alpha)*S@X`` until convergence."""
    if not 0.0 < alpha <= 1.0:
        raise ValueError("alpha must lie in (0, 1].")
    original = rank_values.to_numpy(dtype=float)
    propagated = original.copy()
    isolated = degree == 0

    for iteration in range(1, MAX_ITER + 1):
        updated = alpha * original + (1.0 - alpha) * (smoothing @ propagated)
        updated[isolated, :] = original[isolated, :]
        delta = float(np.max(np.abs(updated - propagated)))
        propagated = updated
        if delta < TOLERANCE:
            return pd.DataFrame(propagated, index=rank_values.index,
                                columns=rank_values.columns), iteration, delta
    raise RuntimeError(f"Propagation did not converge within {MAX_ITER} iterations.")


def random_network_smoothing(n_nodes: int, n_edges: int, seed: int, active=None):
    """Create a simple undirected random graph with matched node/edge counts.

    ``active`` optionally restricts the random edges to those node positions.
    Pass the positions of the genes the real PPI connects so that the control
    randomises *which pairs* are joined while keeping *which genes* carry
    network information; the remaining nodes stay isolated, as in the real
    network.  With ``active=None`` every node may receive edges.
    """
    positions = np.arange(n_nodes) if active is None else np.asarray(active, dtype=int)
    n_active = len(positions)
    if n_active < 2:
        raise ValueError("At least two nodes are required for a random PPI.")
    max_edges = n_active * (n_active - 1) // 2
    if not 0 < n_edges <= max_edges:
        raise ValueError(f"n_edges must lie in [1, {max_edges}], got {n_edges}.")

    rng = np.random.default_rng(seed)
    pairs = set()
    while len(pairs) < n_edges:
        remaining = n_edges - len(pairs)
        draw_n = max(1024, 2 * remaining)
        left = rng.integers(0, n_active, size=draw_n)
        right = rng.integers(0, n_active, size=draw_n)
        for i, j in zip(left, right):
            if i != j:
                pairs.add((min(i, j), max(i, j)))
                if len(pairs) == n_edges:
                    break

    edges = np.asarray(sorted(pairs), dtype=int)
    row = np.concatenate((positions[edges[:, 0]], positions[edges[:, 1]]))
    col = np.concatenate((positions[edges[:, 1]], positions[edges[:, 0]]))
    adjacency = csr_matrix((np.ones(len(row)), (row, col)), shape=(n_nodes, n_nodes))
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    inv_degree = np.divide(1.0, degree, out=np.zeros_like(degree), where=degree > 0)
    return adjacency.multiply(inv_degree[:, None]).tocsr(), degree


def degree_preserving_smoothing(n_nodes: int, edges, seed: int,
                                swaps_per_edge: int = 5):
    """Rewire a simple undirected PPI while preserving every node degree."""
    edge_list = [(int(i), int(j)) for i, j in edges]
    if len(edge_list) < 2:
        raise ValueError("At least two edges are required for degree-preserving rewiring.")
    if swaps_per_edge < 1:
        raise ValueError("swaps_per_edge must be at least one.")

    graph = nx.Graph()
    graph.add_nodes_from(range(n_nodes))
    graph.add_edges_from(edge_list)
    original_degree = np.fromiter((graph.degree(i) for i in range(n_nodes)), dtype=float)
    n_swap = len(edge_list) * swaps_per_edge
    nx.double_edge_swap(graph, nswap=n_swap, max_tries=n_swap * 20, seed=seed)

    rewired_degree = np.fromiter((graph.degree(i) for i in range(n_nodes)), dtype=float)
    if not np.array_equal(original_degree, rewired_degree):
        raise RuntimeError("Degree-preserving PPI rewiring changed node degrees.")

    rewired_edges = np.asarray(list(graph.edges()), dtype=int)
    row = np.concatenate((rewired_edges[:, 0], rewired_edges[:, 1]))
    col = np.concatenate((rewired_edges[:, 1], rewired_edges[:, 0]))
    adjacency = csr_matrix((np.ones(len(row)), (row, col)), shape=(n_nodes, n_nodes))
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    inv_degree = np.divide(1.0, degree, out=np.zeros_like(degree), where=degree > 0)
    return adjacency.multiply(inv_degree[:, None]).tocsr(), degree


def parse_complexes(complex_data: pd.DataFrame) -> dict[str, list[str]]:
    """Parse CORUM ``Complex``/``Genes`` records into unique member lists."""
    required = {"Complex", "Genes"}
    missing = required.difference(complex_data.columns)
    if missing:
        raise ValueError(f"CORUM table missing required column(s): {sorted(missing)}")
    return {
        row.Complex: list(dict.fromkeys(g.strip() for g in str(row.Genes).split(";") if g.strip()))
        for row in complex_data.itertuples(index=False)
    }


def score_complexes(node_values: pd.DataFrame, complexes: dict[str, list[str]]):
    """Score each complex by the mean of its available network members.

    Coverage is returned as a diagnostic only.  It does not determine whether
    a partially represented complex receives a score.
    """
    node_set = set(node_values.index)
    score_rows, coverage_rows = {}, {}
    for name, members_original in complexes.items():
        members = [gene for gene in members_original if gene in node_set]
        coverage = len(members) / len(members_original) if members_original else np.nan
        coverage_rows[name] = coverage
        if len(members) == 0:
            score_rows[name] = pd.Series(np.nan, index=node_values.columns)
        else:
            score_rows[name] = node_values.loc[members].mean(axis=0)

    scores = pd.DataFrame(score_rows).T
    coverage = pd.DataFrame(
        np.repeat(np.array(list(coverage_rows.values()))[:, None], node_values.shape[1], axis=1),
        index=list(coverage_rows), columns=node_values.columns,
    )
    scores.index.name = coverage.index.name = "Complex"
    return scores, coverage
