#!/usr/bin/env python3
"""Build each cohort's zero-missing "stable" protein subset.

Purpose (see README.md): the supplementary cross-reference analysis needs
protein-side complex scores from every method in the primary CPTAC
benchmark (Mean, Min, z-score, GSVA, ssGSEA, PLAGE, NetMean, PerPAS-adapted,
NetComplex). To keep the comparison fair, all nine must see the exact same
protein input -- otherwise a method-specific missing-value rule would be a
second, uncontrolled source of difference between them, on top of the
aggregation logic the analysis is actually trying to isolate.

Rule (deliberately the simplest one, not an imputation): production's own
convention already treats a proteomics value of 0 as "not measured"
(01stringTest/02analysis/03_calc.R: `prot[prot == 0] <- NA`). A gene is kept
here only if it has ZERO such missing entries across every matched sample in
a cohort -- no imputation, no partial-missingness tolerance. This keeps
54-69% of each cohort's measured proteins (verified directly, not assumed)
and leaves 3,300-4,150 of the 4,923 CORUM complexes scoreable per cohort.

Downstream, this stable matrix is also what NetComplex's PPI network gets
built from (netComplex.prepare_network takes its node set directly from the
expression matrix passed to it), so passing this file as the "expression"
input automatically restricts the protein-side network to the same stable
genes -- no separate network-filtering step is needed.
"""

from pathlib import Path

import pandas as pd

CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
PROT_DIR = Path("/proj/c.zihao/work3/00data/cptacT/protein")
OUT_DIR = Path("/proj/c.zihao/work3/01protein/02cor/01stable")

OUT_DIR.mkdir(parents=True, exist_ok=True)

summary_rows = []
for cancer in CANCERS:
    prot = pd.read_csv(PROT_DIR / f"{cancer}_proteomics.csv", index_col=0)
    missing = (prot == 0)
    stable = prot.loc[missing.sum(axis=1) == 0]

    out_path = OUT_DIR / f"{cancer}_protein_stable.csv"
    stable.to_csv(out_path)

    summary_rows.append({
        "cancer": cancer,
        "n_proteins_total": prot.shape[0],
        "n_samples": prot.shape[1],
        "n_proteins_stable": stable.shape[0],
        "pct_stable": round(100 * stable.shape[0] / prot.shape[0], 1),
    })
    print(f"[{cancer}] {prot.shape[0]} -> {stable.shape[0]} proteins "
          f"({summary_rows[-1]['pct_stable']}%), {prot.shape[1]} samples, "
          f"zero missing entries by construction", flush=True)

pd.DataFrame(summary_rows).to_csv(OUT_DIR / "stable_subset_summary.csv", index=False)
print(f"\n[done] stable protein matrices -> {OUT_DIR}")
