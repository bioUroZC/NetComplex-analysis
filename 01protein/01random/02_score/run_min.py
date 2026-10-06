"""
Minimum member gene expression (z-scored across samples).
Loops over all CPTAC cancers, saves min.csv per cancer.
"""

import os
import shutil
import sys

sys.path.insert(0, "/proj/c.zihao/work3/function")
from min import score_min

CANCERS = ["BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC"]
EXPR_DIR  = "/proj/c.zihao/work3/00data/cptacT/exprset"
CORUM     = os.environ["CORUM_PATH"]
OUT_BASE  = os.path.join(os.environ["RAND_OUT_ROOT"], "min")


def reset_out_dir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


reset_out_dir(OUT_BASE)

for cancer in CANCERS:
    expr_path = os.path.join(EXPR_DIR, f"{cancer}_exprSet_filtered.csv")
    out_dir   = os.path.join(OUT_BASE, cancer)
    os.makedirs(out_dir, exist_ok=True)
    print(f"[{cancer}] min ...", flush=True)
    score_min(expr_path, CORUM, out_dir)

print("Done.")
