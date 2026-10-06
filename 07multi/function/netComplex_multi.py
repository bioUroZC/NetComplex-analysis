"""Multi-omics NetComplex via restart fusion on a fixed PPI.

Aligned with the single-omics three-module core in ``netComplex.py``:

1. Within-sample ranks for backbone and auxiliary omics.
2. Fuse ranks into one restart vector; propagate on the **fixed**
   unweighted PPI (same ``prepare_network`` / ``propagate_with_restart``).
3. Score each complex by the mean of available member node values.

Dual-omics information enters **only** through the restart vector.
Edge weights are not expression-modulated. There is no coherence term.

Used by ``01netcomplex/run_multi.py`` to write ``NetComplex_multi`` scores
(``weight_backbone`` defaults to 0.5).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import sys
sys.path.append("/proj/c.zihao/work3/function")
import netComplex


def _validate_aligned(backbone: pd.DataFrame, auxiliary: pd.DataFrame) -> None:
    if not backbone.index.equals(auxiliary.index):
        raise ValueError(
            "Backbone and auxiliary matrices must have identical gene order."
        )
    if not backbone.columns.equals(auxiliary.columns):
        raise ValueError(
            "Backbone and auxiliary matrices must have identical sample order."
        )


def fuse_ranks(
    rank_backbone: pd.DataFrame,
    rank_auxiliary: pd.DataFrame,
    weight_backbone: float = 0.5,
) -> pd.DataFrame:
    """Convex combination of two within-sample rank matrices."""
    if not 0.0 <= weight_backbone <= 1.0:
        raise ValueError("weight_backbone must lie in [0, 1].")
    if not rank_backbone.index.equals(rank_auxiliary.index):
        raise ValueError("Rank matrices must share the same gene index.")
    if not rank_backbone.columns.equals(rank_auxiliary.columns):
        raise ValueError("Rank matrices must share the same sample columns.")
    return (
        weight_backbone * rank_backbone
        + (1.0 - weight_backbone) * rank_auxiliary
    )


def score_multiomics(
    backbone_expression: pd.DataFrame,
    auxiliary_expression: pd.DataFrame,
    links: pd.DataFrame,
    complexes: dict[str, list[str]],
    alpha: float = netComplex.ALPHA,
    weight_backbone: float = 0.5,
):
    """Score complexes with restart-fused multi-omics NetComplex.

    Parameters
    ----------
    backbone_expression, auxiliary_expression
        Pre-aligned gene x sample matrices (identical index/columns).
    links
        PPI edge table as returned by ``netComplex.read_links``.
    complexes
        Parsed CORUM dict from ``netComplex.parse_complexes``.
    alpha
        RWR restart probability (same meaning as single-omics).
    weight_backbone
        Weight on backbone ranks in the fused restart vector.

    Returns
    -------
    complex_scores : pd.DataFrame
        Complex x sample score matrix.
    coverage : pd.DataFrame
        Complex x sample coverage diagnostic.
    node_scores : pd.DataFrame
        Gene x sample propagated node values.
    n_iter : int
        RWR iterations to convergence.
    final_delta : float
        Final max absolute update.
    """
    _validate_aligned(backbone_expression, auxiliary_expression)

    nodes, smoothing, degree = netComplex.prepare_network(
        links, backbone_expression
    )
    rank_backbone = netComplex.within_sample_rank(backbone_expression, nodes)
    rank_auxiliary = netComplex.within_sample_rank(auxiliary_expression, nodes)
    restart = fuse_ranks(rank_backbone, rank_auxiliary, weight_backbone)

    node_scores, n_iter, final_delta = netComplex.propagate_with_restart(
        restart, smoothing, degree, alpha=alpha
    )
    complex_scores, coverage = netComplex.score_complexes(node_scores, complexes)
    return complex_scores, coverage, node_scores, n_iter, final_delta
