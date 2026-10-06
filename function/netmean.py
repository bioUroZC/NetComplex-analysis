
import os
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
    ``score`` column, now constant 1.0, so node strength reduces to degree.
    """
    links = links_df.loc[:, ["protein1", "protein2"]].copy()
    links = links[
        links["protein1"].isin(gene_set) & links["protein2"].isin(gene_set)
    ]
    links = links.drop_duplicates(subset=["protein1", "protein2"])
    links["score"] = 1.0
    return links


def _node_strength_weights(genes, links, eps=1e-6):
    if not genes:
        return pd.Series(dtype=np.float64)

    weights = pd.Series(eps, index=genes, dtype=np.float64)
    if links.empty:
        return weights

    sub = links[
        links["protein1"].isin(genes) & links["protein2"].isin(genes)
    ][["protein1", "protein2", "score"]]
    if sub.empty:
        return weights

    strength = {}
    for _, row in sub.iterrows():
        u = row["protein1"]
        v = row["protein2"]
        w = float(row["score"])
        strength[u] = strength.get(u, 0.0) + w
        strength[v] = strength.get(v, 0.0) + w

    for gene, value in strength.items():
        weights.loc[gene] = value + eps

    return weights


def run_netmean(expr, corum, links):
    """
    Parameters
    ----------
    expr  : pd.DataFrame  — gene x sample (log2 TPM+1)
    corum : pd.DataFrame  — col 0 = complex name, col 1 = genes (;-separated)
    links : pd.DataFrame  — columns = protein1, protein2 (edges unweighted;
                            a `score` column, if present, is ignored)

    Returns
    -------
    pd.DataFrame — complex x sample static-centrality weighted mean z-scores
    """
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    gene_set = set(expr.index)
    geneset = _parse_corum(corum, gene_set)
    links = _prepare_links(links, gene_set)

    E = expr.values.astype(np.float64)
    mu = E.mean(axis=1, keepdims=True)
    sd = E.std(axis=1, keepdims=True) + 1e-10
    Z = (E - mu) / sd
    gene_idx = {g: i for i, g in enumerate(expr.index)}

    names = []
    rows = []
    n_samples = expr.shape[1]

    for name, genes in geneset.items():
        names.append(name)
        if not genes:
            rows.append(np.zeros(n_samples, dtype=np.float64))
            continue

        weights = _node_strength_weights(genes, links)
        gene_ids = [gene_idx[g] for g in genes]
        W = weights.loc[genes].to_numpy(dtype=np.float64).reshape(-1)
        Z_sub = np.atleast_2d(Z[gene_ids])
        denom = W.sum()
        if denom <= 0 or Z_sub.shape[0] != W.shape[0]:
            rows.append(Z_sub.mean(axis=0))
        else:
            rows.append((Z_sub * W[:, None]).sum(axis=0) / denom)

    scores = np.vstack(rows) if rows else np.empty((0, n_samples))
    return pd.DataFrame(scores, index=names, columns=expr.columns)


def score_netmean(expr_path, corum_path, links_path, out_dir):
    """
    Parameters
    ----------
    expr_path  : str — gene x sample CSV (first column = gene symbols)
    corum_path : str — two-column CSV/TSV
    links_path : str — edge list CSV with protein1, protein2 (any `score`
                       column is ignored; edges are unweighted)
    out_dir    : str — output directory; saves netmean.csv
    """
    os.makedirs(out_dir, exist_ok=True)

    expr = pd.read_csv(expr_path, index_col=0)
    sep = "\t" if corum_path.endswith((".tsv", ".txt")) else ","
    corum = pd.read_csv(corum_path, sep=sep)
    links = pd.read_csv(links_path)

    df = run_netmean(expr, corum, links)
    df.round(5).to_csv(os.path.join(out_dir, "netmean.csv"))
    print(f"Saved netmean: {df.shape}")
