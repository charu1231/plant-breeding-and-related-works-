"""
06_summary_figure.py
====================
Combined summary (main) figure for the manuscript, assembling:

  (a) Unstructured GxE simulation — Task A (new environment), h2 = 0.4
  (b) Unstructured GxE simulation — Task B (new lines), h2 = 0.4
  (c) Real CIMMYT wheat — Task A & B model comparison
  (d) Structured (reaction-norm) GxE simulation — Task A (estimated components)
  (e) Unbalanced-data robustness check (30% MCAR missing)
  (f) Real-data estimated type-B genetic correlation (rG) heatmap

Reads output/02_sim_summary.csv, 03_real_validation.json,
04_structured_summary.csv, 05_unbalanced_summary.csv (the last two are optional
if present). Saves output/figures/06_summary.png.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
FIG = OUT / "figures"
sns.set_theme(style="whitegrid", context="notebook")

COL = {"naive": "#999999", "mt": "#2E86AB", "rn": "#A23B72", "single": "#C44536"}


def load_summary():
    s = pd.read_csv(OUT / "02_sim_summary.csv")
    real = json.load(open(OUT / "03_real_validation.json"))
    return s, real


def main():
    s, real = load_summary()
    fig, axes = plt.subplots(2, 3, figsize=(17, 9.5))

    # (a) unstructured Task A, h2=0.4
    ax = axes[0, 0]
    sub = s[s.h2 == 0.4]
    for m, lab in [("naive", "Naive mean"), ("mt", "MT-GBLUP"), ("rn", "RN-GBLUP")]:
        ax.errorbar(sub.rG, sub[f"taskA_{m}_mean"], yerr=sub[f"taskA_{m}_se"],
                    marker="o", label=lab, color=COL[m], capsize=3)
    ax.set_xlabel("Genetic correlation rG"); ax.set_ylabel("Accuracy (r)")
    ax.set_title("(a) Simulation (unstructured G×E): new environment")
    ax.set_ylim(0, 1); ax.legend(frameon=True, fontsize=8)

    # (b) unstructured Task B, h2=0.4
    ax = axes[0, 1]
    for m, lab in [("single", "Single-env GBLUP"), ("mt", "MT-GBLUP")]:
        ax.errorbar(sub.rG, sub[f"taskB_{m}_mean"], yerr=sub[f"taskB_{m}_se"],
                    marker="o", label=lab, color=COL[m], capsize=3)
    ax.set_xlabel("Genetic correlation rG")
    ax.set_title("(b) Simulation (unstructured G×E): new lines")
    ax.set_ylim(0, 1); ax.legend(frameon=True, fontsize=8)

    # (c) real data Task A & B
    ax = axes[0, 2]
    x = np.arange(3); w = 0.32
    ax.bar(x - w, [real["taskA"]["naive"], real["taskA"]["mt"], real["taskA"]["rn"]],
           width=w, color=[COL["naive"], COL["mt"], COL["rn"]], label="Task A: new env")
    ax.bar(x, [real["taskB"]["single"], real["taskB"]["mt"], np.nan],
           width=w, color=[COL["single"], COL["mt"], "#cccccc"], label="Task B: new lines")
    ax.bar(x + w, [real["taskA"]["naive"], real["taskA"]["mt"], real["taskA"]["rn"]],
           width=0, color="none")
    for xi, vals in zip(x - w, [real["taskA"]["naive"], real["taskA"]["mt"], real["taskA"]["rn"]]):
        ax.text(xi, vals + 0.01, f"{vals:.2f}", ha="center", fontsize=8)
    for xi, vals in zip(x, [real["taskB"]["single"], real["taskB"]["mt"], np.nan]):
        if not np.isnan(vals):
            ax.text(xi, vals + 0.01, f"{vals:.2f}", ha="center", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(["Naive /\nsingle-env", "MT-GBLUP", "RN-GBLUP"])
    ax.set_ylim(0, 1); ax.set_title("(c) Real CIMMYT wheat (estimated)")
    ax.legend(frameon=True, fontsize=8)

    # (d) structured G×E Task A
    ax = axes[1, 0]
    p = OUT / "04_structured_summary.csv"
    if p.exists():
        s4 = pd.read_csv(p)
        kappas = sorted(s4.kappa.unique()); h2s = sorted(s4.h2.unique())
        xx = np.arange(len(kappas)); w = 0.28
        for i, (m, lab) in enumerate([("naive", "Naive"), ("mt", "MT-GBLUP"), ("rn", "RN-GBLUP")]):
            vals = [s4[(s4.kappa == k) & (s4.h2 == 0.6)][f"taskA_{m}_mean"].iloc[0] for k in kappas]
            ax.bar(xx + (i - 1) * w, vals, width=w, label=lab, color=COL[m])
        ax.set_xticks(xx); ax.set_xticklabels([f"κ={k}" for k in kappas])
        ax.set_ylim(0, 1); ax.set_xlabel("Reaction-norm κ (slope/intercept)")
        ax.set_ylabel("Accuracy (r)")
        ax.set_title("(d) Structured G×E: new environment (h²=0.6)")
        ax.legend(frameon=True, fontsize=8)
    else:
        ax.text(0.5, 0.5, "pending", ha="center", va="center")
        ax.set_title("(d) Structured G×E")

    # (e) unbalanced
    ax = axes[1, 1]
    p = OUT / "05_unbalanced_summary.csv"
    if p.exists():
        s5 = pd.read_csv(p)
        xx = np.arange(len(s5)); w = 0.32
        for i, (m, lab, c) in enumerate([("naive", "Naive", COL["naive"]),
                                          ("mt", "MT-GBLUP", COL["mt"])]):
            ax.bar(xx + (i - 0.5) * w, s5[f"taskA_{m}_mean"], width=w, label=f"{lab} (A)", color=c)
        ax.bar(xx + 0.5 * w, s5["taskB_single_mean"], width=w, label="Single (B)", color=COL["single"])
        ax.bar(xx + 1.5 * w, s5["taskB_mt_mean"], width=w, label="MT (B)", color=COL["mt"])
        ax.set_xticks(xx + 0.5 * w); ax.set_xticklabels([f"rG={r}" for r in s5.rg])
        ax.set_ylim(0, 1); ax.set_xlabel("Genetic correlation rG (30% missing)")
        ax.set_title("(e) Unbalanced data (h²=0.4)")
        ax.legend(frameon=True, fontsize=7)
    else:
        ax.text(0.5, 0.5, "pending", ha="center", va="center")
        ax.set_title("(e) Unbalanced data")

    # (f) real type-B rG heatmap
    ax = axes[1, 2]
    rB = np.array(real["typeB_correlation"])
    env = real["env_names"]
    sns.heatmap(rB, annot=True, fmt=".2f", cmap="RdBu_r", vmin=-0.1, vmax=1,
                xticklabels=env, yticklabels=env, square=True, ax=ax,
                cbar_kws={"label": "rG"})
    ax.set_title("(f) Wheat type-B genetic correlation")

    fig.suptitle("Multi-environment genomic prediction under G×E — summary",
                 fontsize=15, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(FIG / "06_summary.png", dpi=150)
    plt.close("all")
    print("Saved:", FIG / "06_summary.png")


if __name__ == "__main__":
    main()
