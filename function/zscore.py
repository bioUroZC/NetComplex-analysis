
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


def run_zscore(expr, corum):
    """
    Parameters
    ----------
    expr  : pd.DataFrame  — gene x sample (log2 TPM+1)
    corum : pd.DataFrame  — col 0 = complex name, col 1 = genes (;-separated)

    Returns
    -------
    pd.DataFrame — complex x sample z-score activity
    """
    gene_set = set(expr.index)
    geneset = _parse_corum(corum, gene_set)
    E = expr.values.astype(np.float64)
    gene_idx = {g: i for i, g in enumerate(expr.index)}

    mu = E.mean(axis=1, keepdims=True)
    sd = E.std(axis=1, keepdims=True) + 1e-10
    Z = (E - mu) / sd

    names = list(geneset.keys())
    n_samples = expr.shape[1]
    scores = np.array([
        Z[[gene_idx[g] for g in genes]].sum(axis=0) / np.sqrt(len(genes))
        if genes else np.zeros(n_samples)
        for genes in geneset.values()
    ])
    return pd.DataFrame(scores, index=names, columns=expr.columns)


def score_zscore(expr_path, corum_path, out_dir):
    """
    Parameters
    ----------
    expr_path  : str — gene x sample CSV (first column = gene symbols)
    corum_path : str — two-column CSV/TSV, no header required
    out_dir    : str — output directory; saves zscore.csv
    """
    os.makedirs(out_dir, exist_ok=True)

    expr = pd.read_csv(expr_path, index_col=0)
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    sep = "\t" if corum_path.endswith((".tsv", ".txt")) else ","
    corum = pd.read_csv(corum_path, sep=sep)

    df = run_zscore(expr, corum)
    df.round(5).to_csv(os.path.join(out_dir, "zscore.csv"))
    print(f"Saved zscore: {df.shape}")
