"""
mean scoring for GSE116918.
"""
import os
import shutil
import sys

sys.path.insert(0, "/proj/c.zihao/work3/function")
from mean import score_mean

EXPR = "/proj/c.zihao/work3/03Survival/PRAD/GSE116918/data/exprSet.csv"
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
OUT = "/proj/c.zihao/work3/03Survival/PRAD/GSE116918/bench/mean"

if os.path.isdir(OUT):
    shutil.rmtree(OUT)
os.makedirs(OUT, exist_ok=True)

print(f"[GSE116918] mean ...", flush=True)
score_mean(EXPR, CORUM, OUT)
print("Done.")
