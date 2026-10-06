"""Bulk netmean runtime — same compute path as 01stringTest/01benchmark/run_netmean.py."""

import os
import sys
import time

import pandas as pd

sys.path.insert(0, "/proj/c.zihao/work3/function")
from netmean import run_netmean

DATA_DIR       = "/proj/c.zihao/work3/08runtime/bulk/data"
RESULTS_DIR    = "/proj/c.zihao/work3/08runtime/bulk/results"
LINKS_PATH     = os.path.join(DATA_DIR, "links.csv")
SAMPLE_SIZES   = [10, 20, 50, 100, 200, 500]
COMPLEX_COUNTS = [100, 200, 500, 1000, 2000, 5000]

os.makedirs(RESULTS_DIR, exist_ok=True)
links = pd.read_csv(LINKS_PATH)
links = links.loc[:, ~links.columns.str.startswith("Unnamed")]
records = []

for n in SAMPLE_SIZES:
    print(f"[sample_size={n}] netmean ...", flush=True)
    expr = pd.read_csv(os.path.join(DATA_DIR, f"expr_n{n}.csv"), index_col=0)
    complexes = pd.read_csv(os.path.join(DATA_DIR, "complexes_n1000.csv"))
    t0 = time.time()
    run_netmean(expr, complexes, links)
    elapsed = time.time() - t0
    print(f"  {elapsed:.3f}s")
    records.append(dict(
        sweep="sample_size", value=n,
        n_samples=n, n_complexes=1000,
        runtime_seconds=round(elapsed, 4),
    ))

for n in COMPLEX_COUNTS:
    print(f"[complex_count={n}] netmean ...", flush=True)
    expr = pd.read_csv(os.path.join(DATA_DIR, "expr_n100.csv"), index_col=0)
    complexes = pd.read_csv(os.path.join(DATA_DIR, f"complexes_n{n}.csv"))
    t0 = time.time()
    run_netmean(expr, complexes, links)
    elapsed = time.time() - t0
    print(f"  {elapsed:.3f}s")
    records.append(dict(
        sweep="complex_count", value=n,
        n_samples=100, n_complexes=n,
        runtime_seconds=round(elapsed, 4),
    ))

pd.DataFrame(records).to_csv(
    os.path.join(RESULTS_DIR, f"runtime_netmean_rep{os.environ.get('RUNTIME_REPLICATE', '1')}.csv"), index=False
)
print("Done.")
