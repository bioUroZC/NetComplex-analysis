"""
mean scoring for MTAB386.
"""
import os
import shutil
import sys

sys.path.insert(0, "/proj/c.zihao/work3/function")
from mean import score_mean

EXPR = "/proj/c.zihao/work3/03Survival/OV/MTAB386/data/exprSet.csv"
CORUM = "/proj/c.zihao/work3/00data/corum/complexData.csv"
OUT = "/proj/c.zihao/work3/03Survival/OV/MTAB386/bench/mean"

if os.path.isdir(OUT):
    shutil.rmtree(OUT)
os.makedirs(OUT, exist_ok=True)

print(f"[MTAB386] mean ...", flush=True)
score_mean(EXPR, CORUM, OUT)
print("Done.")
