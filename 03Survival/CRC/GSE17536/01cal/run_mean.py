"""
mean scoring for GSE17536.
"""
import os
import shutil
import sys

sys.path.insert(0, "/proj/c.zihao/work3/function")
from mean import score_mean

EXPR = "/proj/c.zihao/work3/03Survival/CRC/GSE17536/data/exprSet.csv"
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
OUT = "/proj/c.zihao/work3/03Survival/CRC/GSE17536/bench/mean"

if os.path.isdir(OUT):
    shutil.rmtree(OUT)
os.makedirs(OUT, exist_ok=True)

print(f"[GSE17536] mean ...", flush=True)
score_mean(EXPR, CORUM, OUT)
print("Done.")
