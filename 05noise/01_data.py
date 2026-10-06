#!/usr/bin/env python3

import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


EXPR_DIR = Path("/proj/c.zihao/work3/00data/cptacT/exprset")
OUT_ROOT = Path("/proj/c.zihao/work3/05noise/01data")

CANCERS = [
    "BRCA", "COAD", "GBM", "KIRC", "LUAD",
    "LUSC", "OV", "PDAC", "UCEC", "HNSCC",
]
NOISE_LEVELS = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0]
N_SEEDS = 3
SEED_BASE = 1000


parser = argparse.ArgumentParser(description="Build noisy expression matrices")
parser.add_argument("--cancer", required=True, choices=CANCERS)
args = parser.parse_args()
cancer = args.cancer


expression = pd.read_csv(
    EXPR_DIR / f"{cancer}_exprSet_filtered.csv",
    index_col=0,
)
print("[start]", cancer, "genes:", expression.shape[0],
      "samples:", expression.shape[1])


out_dir = OUT_ROOT / cancer
if out_dir.exists():
    shutil.rmtree(out_dir)
out_dir.mkdir(parents=True, exist_ok=True)


for level in NOISE_LEVELS:
    seeds = [0] if level == 0.0 else range(N_SEEDS)

    for seed in seeds:
        if level == 0.0:
            noisy = expression
        else:
            rng = np.random.default_rng(SEED_BASE + seed)
            sigma = expression.std(axis=1).to_numpy()[:, None] * level
            values = expression.to_numpy() + rng.normal(0.0, 1.0, expression.shape) * sigma
            noisy = pd.DataFrame(
                values,
                index=expression.index,
                columns=expression.columns,
            )

        if float(level) == int(level):
            level_str = str(int(level))
        else:
            level_str = str(level)
        out_path = out_dir / f"expr_level{level_str}_seed{seed}.csv"
        noisy.round(6).to_csv(out_path)
        print("  wrote:", out_path.name)

print("[done]", out_dir)
