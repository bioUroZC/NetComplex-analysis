#!/usr/bin/env python3
"""Build nested CPTAC-KIRC expression subsets for the independence test.

Source: full KIRC matrix in 00data/cptacT/exprset.
Writes nested sample subsets at 20 / 40 / 60 / 80% of the cohort (same random
order, so smaller subsets are contained in larger ones), plus a 100% copy as
the reference cohort.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


SRC = Path("/proj/c.zihao/work3/00data/cptacT/exprset/KIRC_exprSet_filtered.csv")
OUT_DIR = Path("/proj/c.zihao/work3/05independence/01data")
FRACTIONS = [0.20, 0.40, 0.60, 0.80, 1.00]
SEED = 2026


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    expr = pd.read_csv(SRC, index_col=0)
    n = expr.shape[1]
    rng = np.random.default_rng(SEED)
    order = rng.permutation(expr.columns.to_numpy())

    rows = []
    for frac in FRACTIONS:
        keep = max(1, int(round(frac * n)))
        keep = min(keep, n)
        samples = list(order[:keep])
        label = f"frac{frac:.2f}"
        out_path = OUT_DIR / f"KIRC_{label}_exprSet.csv"
        expr.loc[:, samples].to_csv(out_path)
        for s in samples:
            rows.append({"fraction": frac, "label": label, "sample": s, "n_samples": keep})
        print(f"[write] {out_path.name}  genes={expr.shape[0]} samples={keep}", flush=True)

    membership = pd.DataFrame(rows)
    membership.to_csv(OUT_DIR / "sample_membership.csv", index=False)

    summary = (
        membership.groupby(["fraction", "label", "n_samples"], as_index=False)
        .size()
        .rename(columns={"size": "n_rows"})
    )
    summary.to_csv(OUT_DIR / "subset_summary.csv", index=False)
    print(f"[done] nested subsets under {OUT_DIR}  (seed={SEED}, N={n})", flush=True)


if __name__ == "__main__":
    main()
