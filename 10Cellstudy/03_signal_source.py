#!/usr/bin/env python3
"""Decompose selected NetComplex score changes into member and network inputs."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.append("/proj/c.zihao/work3/function")
import netComplex

EXPR_ROOT = Path("/proj/c.zihao/work3/00data/singleCell")
CORUM = Path("/proj/c.zihao/work3/00data/corum/complexData.csv")
LINKS = Path("/proj/c.zihao/work3/00data/string/links.csv")
RES_DIR = Path("/proj/c.zihao/work3/10Cellstudy/results")
OUT_DIR = RES_DIR / "signal_source"

RWR_ALPHA = 0.3
TOP_GENES = 20

CONDITIONS = {
    "GSE96583": {
        "GSM2560248": "control",
        "GSM2560249": "stim",
    },
    "GSE226572": {
        "GSM7079694": "control", "GSM7079695": "control", "GSM7079696": "stim",
        "GSM7079702": "control", "GSM7079703": "control", "GSM7079705": "stim",
        "GSM7079710": "control", "GSM7079711": "control", "GSM7079713": "stim",
    },
}

cases = pd.read_csv(RES_DIR / "case_complexes.csv")

OUT_DIR.mkdir(parents=True, exist_ok=True)
links = netComplex.read_links(str(LINKS))
complexes = netComplex.parse_complexes(pd.read_csv(CORUM))

targets = cases["complex"].drop_duplicates().tolist()
summary_rows = []
member_rows = []
for dataset, cond in CONDITIONS.items():

    contrib = {cx: {} for cx in targets}
    restart = {cx: {"stim": [], "control": []} for cx in targets}
    score = {cx: {"stim": [], "control": []} for cx in targets}
    raw_member = {cx: {g: {"stim": [], "control": []}
                       for g in complexes[cx]} for cx in targets}

    for sample, group in cond.items():
        expr_path = EXPR_ROOT / dataset / "out" / sample / "logNorm.txt"
        expr = pd.read_csv(expr_path, sep="\t", index_col=0)
        expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        nodes, smoothing, degree = netComplex.prepare_network(links, expr)
        ranks = netComplex.within_sample_rank(expr, nodes)
        node_scores, _, _ = netComplex.propagate_with_restart(
            ranks, smoothing, degree, alpha=RWR_ALPHA
        )
        pos = {gene: i for i, gene in enumerate(nodes)}
        smoothing = smoothing.tocsr()
        p_mean = node_scores.mean(axis=1)
        r_mean = ranks.mean(axis=1)
        print(f"[{dataset}/{sample}] {group}: {expr.shape[0]} genes x "
              f"{expr.shape[1]} cells", flush=True)

        for cx in targets:
            members = [g for g in complexes[cx] if g in pos]
            n_members = len(members)
            score[cx][group].append(
                float(np.mean([p_mean.iloc[pos[g]] for g in members])))
            restart[cx][group].append(
                float(np.sum([RWR_ALPHA * r_mean.iloc[pos[g]] for g in members])) / n_members)
            for gene in complexes[cx]:
                measured = gene in expr.index
                values = (expr.loc[gene].to_numpy(dtype=float)
                          if measured else np.zeros(expr.shape[1]))
                raw_member[cx][gene][group].append({
                    "mean_expression": float(np.mean(values)),
                    "detected_fraction": float(np.mean(values > 0)),
                    "measured": measured,
                })
            for gene in members:
                i = pos[gene]
                if degree[i] == 0:
                    continue
                start, end = smoothing.indptr[i], smoothing.indptr[i + 1]
                for k in range(start, end):
                    j = smoothing.indices[k]
                    add = ((1.0 - RWR_ALPHA) * float(smoothing.data[k])
                           * float(p_mean.iloc[j]) / n_members)
                    slot = contrib[cx].setdefault(nodes[j], {"stim": [], "control": []})
                    slot[group].append(add)

    for cx in targets:
        members = set(complexes[cx])
        case_members = []
        for gene, values in raw_member[cx].items():
            member = {
                "dataset": dataset, "complex": cx, "gene": gene,
                "mean_expression_control": float(np.mean(
                    [item["mean_expression"] for item in values["control"]]
                )),
                "mean_expression_stim": float(np.mean(
                    [item["mean_expression"] for item in values["stim"]]
                )),
                "detection_control": float(np.mean(
                    [item["detected_fraction"] for item in values["control"]]
                )),
                "detection_stim": float(np.mean(
                    [item["detected_fraction"] for item in values["stim"]]
                )),
                "measured_sample_fraction_control": float(np.mean(
                    [item["measured"] for item in values["control"]]
                )),
                "measured_sample_fraction_stim": float(np.mean(
                    [item["measured"] for item in values["stim"]]
                )),
            }
            member["expression_delta"] = member["mean_expression_stim"] - member["mean_expression_control"]
            member["detection_delta"] = member["detection_stim"] - member["detection_control"]
            member_rows.append(member)
            case_members.append(member)
        rows = []
        for gene, v in contrib[cx].items():
            s = float(np.mean(v["stim"])) if v["stim"] else 0.0
            c = float(np.mean(v["control"])) if v["control"] else 0.0
            rows.append({"gene": gene, "contribution_stim": s,
                         "contribution_control": c, "delta": s - c,
                         "is_member": gene in members})
        tab = pd.DataFrame(rows).sort_values("delta", ascending=False)
        safe = cx.replace("/", "_")
        tab.head(TOP_GENES).round(6).to_csv(OUT_DIR / f"{dataset}__{safe}.csv",
                                            index=False)

        s_score = float(np.mean(score[cx]["stim"]))
        c_score = float(np.mean(score[cx]["control"]))
        ext = tab.loc[~tab["is_member"], "delta"].sum()
        mem = tab.loc[tab["is_member"], "delta"].sum()
        rst = np.mean(restart[cx]["stim"]) - np.mean(restart[cx]["control"])
        total = ext + mem + rst
        detection_control = np.nanmedian([x["detection_control"] for x in case_members])
        detection_stim = np.nanmedian([x["detection_stim"] for x in case_members])
        expression_delta = np.nanmedian([x["expression_delta"] for x in case_members])
        summary_rows.append({
            "dataset": dataset, "complex": cx,
            "n_members_in_network": len([g for g in complexes[cx]]),
            "score_stim": s_score, "score_control": c_score,
            "score_delta": s_score - c_score,
            "delta_from_own_rank": rst,
            "delta_from_member_neighbours": mem,
            "delta_from_external_partners": ext,
            "external_share_of_delta": ext / total,
            "median_member_detection_control": detection_control,
            "median_member_detection_stim": detection_stim,
            "median_member_expression_delta": expression_delta,
            "top_external_drivers": ";".join(
                tab.loc[~tab["is_member"], "gene"].head(5).tolist()),
        })
        print(f"  {cx}: score {c_score:.4f} -> {s_score:.4f}, "
              f"external share of change {ext / total:.2f}", flush=True)

pd.DataFrame(summary_rows).round(6).to_csv(RES_DIR / "signal_source_summary.csv", index=False)
pd.DataFrame(member_rows).round(6).to_csv(RES_DIR / "member_signal_summary.csv", index=False)
print(f"\nWrote {RES_DIR / 'signal_source_summary.csv'} ({len(summary_rows)} rows)")
print(f"Wrote {RES_DIR / 'member_signal_summary.csv'} ({len(member_rows)} rows)")
