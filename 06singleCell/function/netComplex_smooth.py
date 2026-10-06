"""NetComplex single-cell smooth variant.

KNN-PCA expression smoothing, then the classic three-module NetComplex:
within-sample rank -> restart RWR -> complex mean.
"""

from __future__ import annotations

import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

import sys
sys.path.append("/proj/c.zihao/work3/function")
import netComplex


def smooth_expression_knn(
    expression_data: pd.DataFrame,
    n_neighbors: int = 20,
    lam: float = 0.5,
    n_pcs: int = 50,
    random_state: int = 0,
) -> pd.DataFrame:
    """Smooth expression by mixing each cell with k nearest neighbors in PCA space.

    ``x_smooth(c) = lam * x(c) + (1 - lam) * mean_{c' in N(c)} x(c')``
    """
    if not 0.0 <= lam <= 1.0:
        raise ValueError(f"lam must lie in [0, 1], got {lam}")

    X = expression_data.to_numpy(dtype=float).T                 
    n_cells, n_genes = X.shape
    if n_cells < 2:
        return expression_data.copy()

    actual_pcs = min(n_pcs, n_cells - 1, n_genes - 1)
    if actual_pcs < 1:
        return expression_data.copy()

    X_pca = PCA(n_components=actual_pcs, random_state=random_state).fit_transform(X)
    k = min(n_neighbors, n_cells - 1)
    _, indices = NearestNeighbors(n_neighbors=k, metric="euclidean").fit(X_pca).kneighbors(X_pca)

    neighbor_mean = X[indices].mean(axis=1)
    X_smoothed = lam * X + (1.0 - lam) * neighbor_mean
    return pd.DataFrame(
        X_smoothed.T,
        index=expression_data.index,
        columns=expression_data.columns,
    )


def score_cells_smooth(
    links: pd.DataFrame,
    expression: pd.DataFrame,
    complexes: dict,
    alpha: float = netComplex.ALPHA,
    n_neighbors: int = 20,
    lam: float = 0.5,
    n_pcs: int = 50,
    random_state: int = 0,
):
    """KNN-PCA smooth expression, then score each cell with classic NetComplex."""
    smoothed = smooth_expression_knn(
        expression,
        n_neighbors=n_neighbors,
        lam=lam,
        n_pcs=n_pcs,
        random_state=random_state,
    )
    nodes, smoothing, degree = netComplex.prepare_network(links, smoothed)
    ranks = netComplex.within_sample_rank(smoothed, nodes)
    node_scores, n_iter, final_delta = netComplex.propagate_with_restart(
        ranks, smoothing, degree, alpha=alpha
    )
    complex_scores, coverage = netComplex.score_complexes(node_scores, complexes)
    meta = {
        "n_nodes": len(nodes),
        "n_columns": smoothed.shape[1],
        "n_complexes": complex_scores.shape[0],
        "n_iter": n_iter,
        "final_delta": final_delta,
        "alpha": alpha,
        "n_neighbors": n_neighbors,
        "lam": lam,
        "n_pcs": n_pcs,
    }
    return complex_scores, coverage, node_scores, meta
