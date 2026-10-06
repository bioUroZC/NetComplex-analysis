"""Cell NetComplex runtime — mirrors the classic Figure 6 scoring workflow.

Timing covers gene-set parsing as well as scoring: every baseline parses the
complex table inside its ``run_*`` call, so leaving ``parse_complexes`` outside
``t0`` would credit NetComplex with work the baselines are charged for.
"""

import os
import sys
import time

import pandas as pd

sys.path.insert(0, "/proj/c.zihao/work3/function")
import netComplex as nc

DATA_DIR       = "/proj/c.zihao/work3/08runtime/cell/data"
RESULTS_DIR    = "/proj/c.zihao/work3/08runtime/cell/results"
LINKS_PATH     = os.path.join(DATA_DIR, "links.csv")
CELL_COUNTS    = [500, 1000, 2000, 5000, 10000]
COMPLEX_COUNTS = [100, 200, 500, 1000, 2000, 5000]

RWR_ALPHA   = 0.3

os.makedirs(RESULTS_DIR, exist_ok=True)

links = nc.read_links(LINKS_PATH)
records = []

for n in CELL_COUNTS:
    print(f"[cell_count={n}] netcomplex ...", flush=True)
    expr = pd.read_csv(os.path.join(DATA_DIR, f"expr_n{n}.csv"), index_col=0)
    complex_table = pd.read_csv(os.path.join(DATA_DIR, "complexes_n1000.csv"))
    t0 = time.time()
    complexes = nc.parse_complexes(complex_table)
    nodes, smoothing, degree = nc.prepare_network(links, expr)
    ranks = nc.within_sample_rank(expr, nodes)
    node_scores, _, _ = nc.propagate_with_restart(
        ranks, smoothing, degree, alpha=RWR_ALPHA
    )
    nc.score_complexes(node_scores, complexes)
    elapsed = time.time() - t0
    print(f"  {elapsed:.3f}s")
    records.append(dict(
        sweep="cell_count", value=n,
        n_cells=n, n_complexes=1000,
        runtime_seconds=round(elapsed, 4),
    ))

for n in COMPLEX_COUNTS:
    print(f"[complex_count={n}] netcomplex ...", flush=True)
    expr = pd.read_csv(os.path.join(DATA_DIR, "expr_n2000.csv"), index_col=0)
    complex_table = pd.read_csv(os.path.join(DATA_DIR, f"complexes_n{n}.csv"))
    t0 = time.time()
    complexes = nc.parse_complexes(complex_table)
    nodes, smoothing, degree = nc.prepare_network(links, expr)
    ranks = nc.within_sample_rank(expr, nodes)
    node_scores, _, _ = nc.propagate_with_restart(
        ranks, smoothing, degree, alpha=RWR_ALPHA
    )
    nc.score_complexes(node_scores, complexes)
    elapsed = time.time() - t0
    print(f"  {elapsed:.3f}s")
    records.append(dict(
        sweep="complex_count", value=n,
        n_cells=2000, n_complexes=n,
        runtime_seconds=round(elapsed, 4),
    ))

pd.DataFrame(records).to_csv(
    os.path.join(RESULTS_DIR, f"runtime_netcomplex_rep{os.environ.get('RUNTIME_REPLICATE', '1')}.csv"), index=False
)
print("Done.")
