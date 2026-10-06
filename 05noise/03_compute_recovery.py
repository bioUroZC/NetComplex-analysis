#!/usr/bin/env python3

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr


RESULT_DIR = Path("/proj/c.zihao/work3/05noise/results")
SCAN_DIR = RESULT_DIR / "scan"

CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
METHOD_DIRS = {
    "netcomplex": RESULT_DIR / "netcomplex" / "complex_score",
    "noRWR": RESULT_DIR / "noRWR",
    "RewiredPPI": RESULT_DIR / "RewiredPPI" / "complex_score",
}
NETWORK = "string"
TAG_RE = re.compile(r"level(?P<level>[0-9.]+)_seed(?P<seed>[0-9]+)")
LEVEL0_RECOVERY_LO = 0.999
LEVEL0_RECOVERY_HI = 1.001


def complex_recovery(noisy, clean):
    a, b = noisy.to_numpy(), clean.to_numpy()
    keep = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
    keep &= (a.std(axis=1) > 0) & (b.std(axis=1) > 0)
    if keep.sum() == 0:
        return np.nan
    rho = [spearmanr(x, y).statistic for x, y in zip(a[keep], b[keep])]
    return float(np.nanmedian(rho))


def geometry_recovery(noisy, clean):
    a, b = noisy.to_numpy(), clean.to_numpy()
    keep = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
    keep &= (a.std(axis=1) > 0) & (b.std(axis=1) > 0)
    if keep.sum() < 2:
        return np.nan

    def dists(m):
        z = (m - m.mean(axis=1, keepdims=True)) / m.std(axis=1, keepdims=True)
        return pdist(z.T, metric="euclidean")

    da, db = dists(a[keep]), dists(b[keep])
    if len(da) < 3:
        return np.nan
    return float(spearmanr(da, db).statistic)


def align_to_clean(mat, clean, path):
    """Reindex noisy scores onto the clean complex/sample axes; fail on mismatch."""
    if not mat.index.equals(clean.index) or not mat.columns.equals(clean.columns):
        missing_rows = clean.index.difference(mat.index)
        extra_rows = mat.index.difference(clean.index)
        missing_cols = clean.columns.difference(mat.columns)
        extra_cols = mat.columns.difference(clean.columns)
        if len(missing_rows) or len(extra_rows) or len(missing_cols) or len(extra_cols):
            raise ValueError(
                f"Score matrix axes differ from clean reference: {path}\n"
                f"  missing complexes: {len(missing_rows)}, "
                f"extra complexes: {len(extra_rows)}, "
                f"missing samples: {len(missing_cols)}, "
                f"extra samples: {len(extra_cols)}"
            )
        mat = mat.reindex(index=clean.index, columns=clean.columns)
    return mat


def score_file(method, cancer, tag):
    root = METHOD_DIRS[method]
    if method == "noRWR":
        return root / cancer / f"{cancer}_{tag}_noRWR_complex_score.csv"
    return root / cancer / f"{cancer}_{tag}_{method}_complex_score.csv"


def list_tags(cancer):
    score_dir = METHOD_DIRS["netcomplex"] / cancer
    tags = []
    for path in sorted(score_dir.glob(f"{cancer}_level*_seed*_netcomplex_complex_score.csv")):
        match = TAG_RE.search(path.name)
        if match:
            tags.append(f"level{match.group('level')}_seed{match.group('seed')}")
    return tags


def parse_tag(tag):
    match = TAG_RE.fullmatch(tag)
    if not match:
        raise ValueError(f"Bad tag: {tag}")
    return float(match.group("level")), int(match.group("seed"))


def check_level0_recovery(cancer, method, complex_rec, geometry_rec):
    for name, value in (("complex_recovery", complex_rec),
                        ("geometry_recovery", geometry_rec)):
        if not np.isfinite(value) or not (LEVEL0_RECOVERY_LO <= value <= LEVEL0_RECOVERY_HI):
            raise ValueError(
                f"level0_seed0 {name} out of [{LEVEL0_RECOVERY_LO}, {LEVEL0_RECOVERY_HI}] "
                f"for {cancer}/{method}: {value}"
            )


SCAN_DIR.mkdir(parents=True, exist_ok=True)
all_rows = []

for cancer in CANCERS:
    tags = list_tags(cancer)
    if not tags:
        print("[skip]", cancer, ": no score files")
        continue

    rows = []
    for method in METHOD_DIRS:
        clean_path = score_file(method, cancer, "level0_seed0")
        if not clean_path.exists():
            raise FileNotFoundError(f"Missing clean scores: {clean_path}")
        clean = pd.read_csv(clean_path, index_col=0)

        for tag in tags:
            level, seed = parse_tag(tag)
            path = score_file(method, cancer, tag)
            mat = align_to_clean(pd.read_csv(path, index_col=0), clean, path)
            c_rec = complex_recovery(mat, clean)
            g_rec = geometry_recovery(mat, clean)
            if tag == "level0_seed0":
                check_level0_recovery(cancer, method, c_rec, g_rec)
            rows.append({
                "network": NETWORK,
                "cancer": cancer,
                "noise_level": level,
                "seed": seed,
                "method": method,
                "n_complex": int(mat.notna().any(axis=1).sum()),
                "complex_recovery": c_rec,
                "geometry_recovery": g_rec,
            })
        print("[ok]", cancer, method)

    out = SCAN_DIR / f"{cancer}_noise_scan.csv"
    pd.DataFrame(rows).round(6).to_csv(out, index=False)
    all_rows.extend(rows)
    print("  saved", out)

if not all_rows:
    raise SystemExit("No recovery rows written — run 02_*.py scoring scripts first.")

print("[done]", len(all_rows), "rows under", SCAN_DIR)
