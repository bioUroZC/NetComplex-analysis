
import os
import heapq
import numpy as np
import pandas as pd


def _parse_corum(corum_df, gene_set):
    geneset = {}
    for _, row in corum_df.iterrows():
        name = str(row.iloc[0])
        genes = [g.strip() for g in str(row.iloc[1]).split(";") if g.strip()]
        genes = [g for g in genes if g in gene_set]
        geneset[name] = genes
    return geneset


def _prepare_links(links_df, gene_set):
    """Keep edges between expressed genes and treat every edge as unweighted.

    Any ``score`` column in the input is deliberately ignored: ``netComplex``
    drops it in ``read_links``, so weighting this baseline by STRING
    confidence while NetComplex runs on an unweighted graph would confound
    every comparison between them.  The returned frame still carries a
    ``score`` column, now constant 1.0, so node strength reduces to degree
    and the weighted betweenness below reduces to hop-count betweenness.
    """
    links = links_df.loc[:, ["protein1", "protein2"]].copy()
    links = links[
        links["protein1"].isin(gene_set) & links["protein2"].isin(gene_set)
    ]
    links = links.drop_duplicates(subset=["protein1", "protein2"])
    links["score"] = 1.0
    return links


def _normalized(series):
    if series.empty:
        return series.astype(np.float64)
    max_val = float(series.max())
    if max_val <= 0:
        return pd.Series(0.0, index=series.index, dtype=np.float64)
    return series.astype(np.float64) / max_val


def _build_complex_graph(genes, links):
    graph = {gene: {} for gene in genes}
    if not genes or links.empty:
        return graph

    sub = links[
        links["protein1"].isin(genes) & links["protein2"].isin(genes)
    ][["protein1", "protein2", "score"]]
    for _, row in sub.iterrows():
        u = row["protein1"]
        v = row["protein2"]
        weight = float(row["score"])
        if weight <= 0:
            continue
        graph[u][v] = weight
        graph[v][u] = weight
    return graph


def _weighted_betweenness(graph):
    nodes = list(graph.keys())
    betweenness = dict.fromkeys(nodes, 0.0)

    for source in nodes:
        stack = []
        preds = {node: [] for node in nodes}
        sigma = dict.fromkeys(nodes, 0.0)
        sigma[source] = 1.0
        dist = dict.fromkeys(nodes, np.inf)
        dist[source] = 0.0

        queue = [(0.0, source)]
        while queue:
            dist_v, v = heapq.heappop(queue)
            if dist_v > dist[v]:
                continue
            stack.append(v)

            for nbr, weight in graph[v].items():
                edge_dist = 1.0 / (weight + 1e-12)
                alt = dist[v] + edge_dist
                if alt + 1e-12 < dist[nbr]:
                    dist[nbr] = alt
                    heapq.heappush(queue, (alt, nbr))
                    sigma[nbr] = sigma[v]
                    preds[nbr] = [v]
                elif abs(alt - dist[nbr]) <= 1e-12:
                    sigma[nbr] += sigma[v]
                    preds[nbr].append(v)

        delta = dict.fromkeys(nodes, 0.0)
        while stack:
            w = stack.pop()
            if sigma[w] > 0:
                coeff = (1.0 + delta[w]) / sigma[w]
                for v in preds[w]:
                    delta[v] += sigma[v] * coeff
            if w != source:
                betweenness[w] += delta[w]

    for node in betweenness:
        betweenness[node] /= 2.0

    return pd.Series(betweenness, dtype=np.float64)


def _topology_weights(graph):
    genes = list(graph.keys())
    if not genes:
        return pd.Series(dtype=np.float64)

    strength = pd.Series(0.0, index=genes, dtype=np.float64)
    edge_count = 0
    for gene in genes:
        if graph[gene]:
            strength.loc[gene] = sum(float(w) for w in graph[gene].values())
            edge_count += len(graph[gene])

    strength_norm = _normalized(strength)

    if edge_count == 0:
        betweenness_norm = pd.Series(0.0, index=genes, dtype=np.float64)
    else:
        bet = _weighted_betweenness(graph)
        betweenness_norm = _normalized(bet)

    topology = 0.5 * (strength_norm + betweenness_norm)
    if float(topology.sum()) <= 0:
        topology[:] = 1.0
    return topology


def _local_activity(graph, z_values, gene):
    numer = float(z_values.get(gene, 0.0))
    denom = 1.0
    for nbr, weight in graph.get(gene, {}).items():
        weight = float(weight)
        numer += weight * float(z_values.get(nbr, 0.0))
        denom += weight
    return numer / denom if denom > 0 else 0.0


def run_perpas_adapted(expr, corum, links):
    """
    Parameters
    ----------
    expr  : pd.DataFrame  — gene x sample (log2 TPM+1)
    corum : pd.DataFrame  — col 0 = complex name, col 1 = genes (;-separated)
    links : pd.DataFrame  — columns = protein1, protein2 (edges unweighted;
                            a `score` column, if present, is ignored)

    Returns
    -------
    pd.DataFrame — complex x sample topology-weighted local activity scores
    """
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    gene_set = set(expr.index)
    geneset = _parse_corum(corum, gene_set)
    links = _prepare_links(links, gene_set)

    E = expr.values.astype(np.float64)
    mu = E.mean(axis=1, keepdims=True)
    sd = E.std(axis=1, keepdims=True) + 1e-10
    Z = pd.DataFrame((E - mu) / sd, index=expr.index, columns=expr.columns)

    out = pd.DataFrame(0.0, index=list(geneset.keys()), columns=expr.columns)

    for complex_name, genes in geneset.items():
        if not genes:
            continue

        graph = _build_complex_graph(genes, links)
        topology = _topology_weights(graph).reindex(genes).fillna(0.0)
        topo_sum = float(topology.sum())
        if topo_sum <= 0:
            topology[:] = 1.0
            topo_sum = float(topology.sum())

        for sample in expr.columns:
            z_values = Z[sample]
            local_scores = np.array(
                [_local_activity(graph, z_values, gene) for gene in genes],
                dtype=np.float64
            )
            weights = topology.to_numpy(dtype=np.float64)
            out.loc[complex_name, sample] = float(
                np.dot(weights, local_scores) / topo_sum
            )

    return out


def score_perpas_adapted(expr_path, corum_path, links_path, out_dir):
    """
    Parameters
    ----------
    expr_path  : str — gene x sample CSV (first column = gene symbols)
    corum_path : str — two-column CSV/TSV
    links_path : str — edge list CSV with protein1, protein2 (any `score`
                       column is ignored; edges are unweighted)
    out_dir    : str — output directory; saves perpas_adapted.csv
    """
    os.makedirs(out_dir, exist_ok=True)

    expr = pd.read_csv(expr_path, index_col=0)
    sep = "\t" if corum_path.endswith((".tsv", ".txt")) else ","
    corum = pd.read_csv(corum_path, sep=sep)
    links = pd.read_csv(links_path)

    df = run_perpas_adapted(expr, corum, links)
    df.round(5).to_csv(os.path.join(out_dir, "perpas_adapted.csv"))
    print(f"Saved perpas_adapted: {df.shape}")
