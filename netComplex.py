import numpy as np
import pandas as pd
import networkx as nx
from scipy.sparse import csr_matrix

# NetComplex current workflow:
# 1. retain genes shared by the PPI network and expression matrix;
# 2. rank-normalize one sample at a time;
# 3. update sample-specific PPI edge scores using endpoint expression ranks;
# 4. run RWR on the sample-specific weighted network to obtain node scores;
# 5. aggregate node scores into complex-level base activity, coherence, and final score.
#
# Current complex scoring definition:
# - BaseActivity = expression-weighted average of propagated node scores;
# - Coherence = average pairwise node-score similarity across retained complex members
#   multiplied by a size penalty k / (k + 1);
# - ComplexScore = BaseActivity * Coherence.


def rank_normalize_series(sample_expression):
    """Rank-normalize one sample so higher expression receives higher values in (0, 1]."""
    n_genes = len(sample_expression)
    if n_genes == 0:
        return sample_expression.copy()

    ranks = sample_expression.rank(method='average', ascending=True)
    return ranks / n_genes


def prepare_data_netCom(links, expression_data):
    """Retain genes shared by the PPI and expression matrix, then build the base PPI graph."""
    genes_in_links = pd.unique(links[['protein1', 'protein2']].values.ravel())
    common_genes = expression_data.index.intersection(genes_in_links)
    link_filtered = links[
        links['protein1'].isin(common_genes) & links['protein2'].isin(common_genes)
    ].drop_duplicates(subset=['protein1', 'protein2'])
    links = link_filtered
    used_genes = pd.unique(link_filtered[['protein1', 'protein2']].values.ravel())
    expression_data = expression_data.loc[expression_data.index.isin(used_genes)]

    G = nx.Graph()
    for _, row in links.iterrows():
        G.add_edge(row['protein1'], row['protein2'], weight=row['score'])

    real_original_weights = {
        (row['protein1'], row['protein2']): row['score']
        for _, row in links.iterrows()
    }

    return G, real_original_weights, expression_data


def run_rwr(sample_id, G, real_original_weights, expression_data, rwr_alpha=0.3, beta=0.5):
    """Build a sample-specific weighted PPI graph and run RWR to obtain node scores."""
    for (u, v), weight in real_original_weights.items():
        if G.has_edge(u, v):
            G[u][v]['weight'] = weight
            G[u][v]['adjusted_weight'] = weight
            G[u][v]['correction_score'] = 1.0

    sample_expression = expression_data[sample_id]
    expr_norm = rank_normalize_series(sample_expression).copy()

    for u, v, data in G.edges(data=True):
        zu = expr_norm.get(u, 0.0)
        zv = expr_norm.get(v, 0.0)
        BaseWeight = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))
        data['BaseWeight'] = BaseWeight
        # Sample-specific edge score = original PPI score scaled by the two endpoint expression ranks.
        adjusted_weight = BaseWeight * ((zu * zv) ** beta)
        data['adjusted_weight'] = adjusted_weight
        data['correction_score'] = (zu * zv) ** beta

    nodes = list(G.nodes)
    n = len(nodes)
    node_to_index = {node: i for i, node in enumerate(nodes)}

    row, col, edge_weights = [], [], []
    for u, v in G.edges():
        weight = G[u][v]['adjusted_weight']
        row.append(node_to_index[u])
        col.append(node_to_index[v])
        edge_weights.append(weight)
        row.append(node_to_index[v])
        col.append(node_to_index[u])
        edge_weights.append(weight)

    T_sparse = csr_matrix((edge_weights, (row, col)), shape=(n, n))
    col_sums = np.array(T_sparse.sum(axis=0)).flatten()
    col_sums[col_sums == 0] = 1
    T_sparse = T_sparse.multiply(1 / col_sums)

    P0 = np.array([expr_norm.get(node, 0.0) for node in nodes], dtype=float)
    if P0.sum() == 0:
        P0[:] = 1.0
    P0 = P0 / P0.sum()
    P = P0.copy()

    # Propagate the sample-specific restart signal over the weighted network until convergence.
    for _ in range(50):
        P_new = rwr_alpha * (T_sparse @ P) + (1 - rwr_alpha) * P0
        if np.linalg.norm(P - P_new) < 1e-4:
            P = P_new
            break
        P = P_new

    P_normalized = (
        pd.Series(P)
        .rank(method='average', ascending=True)
        .to_numpy(dtype=float) / len(P)
    )

    node_values_df = pd.DataFrame({
        "Sample": sample_id,
        "Node": nodes,
        "Expression": [sample_expression.get(node, np.nan) for node in nodes],
        "normExpr": [expr_norm.get(node, 0.0) for node in nodes],
        "NodeValue": P_normalized,
    })
    node_values_df["normExpr"] = node_values_df["normExpr"].round(5)
    node_values_df["NodeValue"] = node_values_df["NodeValue"].round(5)

    p_star_series = pd.Series(P_normalized, index=nodes)

    return node_values_df, p_star_series, expr_norm


def compute_complex_score(sample_id, p_star, expression_ranked, complex_data,
                          eps=1e-12, gamma=1.0):
    """Aggregate node scores into complex-level base activity, coherence, and final score.

    Parameters
    ----------
    gamma : float, default 1.0
        Exponent applied to the mean pairwise node-score similarity.
        Coherence = (mean_pair_score)^gamma * size_penalty.
        gamma=0 removes pairwise consistency from coherence (size penalty only).
        gamma=1 (default) uses the raw mean pairwise similarity.
        gamma>1 amplifies the contribution of high pairwise consistency.
    """
    complex_records = []

    for _, row in complex_data.iterrows():
        complex_name = row["Complex"]
        complex_genes = [
            gene.strip() for gene in str(row["Genes"]).split(";")
            if gene.strip()
        ]
        retained_genes = [gene for gene in complex_genes if gene in p_star.index]
        complex_size = len(retained_genes)
        size_penalty = complex_size / (complex_size + 1) if complex_size > 0 else 0.0
        base_activity = 0.0
        coherence = 0.0
        score = 0.0

        if retained_genes:
            node_scores = p_star.loc[retained_genes]
            expression_weights = expression_ranked.loc[retained_genes]

            if node_scores.sum() > eps and expression_weights.sum() > eps:
                # Base activity is the expression-weighted average of propagated node scores.
                base_activity = (node_scores * expression_weights).sum() / (expression_weights.sum() + eps)

                if complex_size == 1:
                    coherence = size_penalty
                else:
                    pair_scores = []
                    for i in range(complex_size):
                        for j in range(i + 1, complex_size):
                            gene_a = retained_genes[i]
                            gene_b = retained_genes[j]
                            node_score_a = p_star.get(gene_a, 0.0)
                            node_score_b = p_star.get(gene_b, 0.0)
                            pair_score = (2 * min(node_score_a, node_score_b)) / (node_score_a + node_score_b + eps)
                            pair_scores.append(pair_score)
                    # Coherence captures network-supported overall coordination across retained members.
                    mean_pair = float(np.mean(pair_scores)) if pair_scores else 0.0
                    coherence = (mean_pair ** gamma) * size_penalty

                score = base_activity * coherence

        complex_records.append({
            "Sample": sample_id,
            "Complex": complex_name,
            "Genes": ";".join(complex_genes),
            "RetainedGenes": ";".join(retained_genes),
            "ComplexSize": complex_size,
            "SizePenalty": size_penalty,
            "BaseActivity": base_activity,
            "Coherence": coherence,
            "ComplexScore": score,
        })

    complex_scores_df = pd.DataFrame(complex_records)
    for col in ["SizePenalty", "BaseActivity", "Coherence", "ComplexScore"]:
        complex_scores_df[col] = complex_scores_df[col].round(5)

    return complex_scores_df.reset_index(drop=True)


def compute_all_sample_complex_scores(links, expression_data, complex_data,
                                      rwr_alpha=0.3, beta=0.5, gamma=1.0, eps=1e-12):
    """Compute one long-format complex-score table containing all samples."""
    G, real_original_weights, expression_data = prepare_data_netCom(
        links=links,
        expression_data=expression_data
    )

    all_complex_scores = []
    for sample_id in expression_data.columns:
        _, p_star_series, expr_norm = run_rwr(
            sample_id=sample_id,
            G=G,
            real_original_weights=real_original_weights,
            expression_data=expression_data,
            rwr_alpha=rwr_alpha,
            beta=beta
        )
        complex_scores_df = compute_complex_score(
            sample_id=sample_id,
            p_star=p_star_series,
            expression_ranked=expr_norm,
            complex_data=complex_data,
            eps=eps,
            gamma=gamma
        )
        all_complex_scores.append(complex_scores_df)

    if not all_complex_scores:
        return pd.DataFrame(
            columns=[
                "Sample", "Complex", "Genes", "RetainedGenes", "ComplexSize",
                "SizePenalty", "BaseActivity", "Coherence", "ComplexScore"
            ]
        )

    return pd.concat(all_complex_scores, ignore_index=True)


def build_complex_score_matrix(all_complex_scores):
    """Convert the long-format complex-score table to a complex-by-sample matrix."""
    complex_score_matrix = all_complex_scores.pivot(
        index="Complex",
        columns="Sample",
        values="ComplexScore"
    )
    return complex_score_matrix.reset_index()
