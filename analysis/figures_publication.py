"""
figures_publication.py
======================
Regenerate publication-quality figures (300 DPI) directly from the saved
CSV/JSON outputs (no re-simulation), with a consistent journal style.

Produces (in output/figures_pub/):
  fig1_simulation.png        main benchmark (unstructured / compound symmetry)
  fig2_real_data.png         real-data validation (type-B rG + tasks)
  fig3_robustness.png        structured GxE + unbalanced + factor-analytic
  fig4_summary.png           compact 6-panel overview

Usage: .venv/bin/python analysis/figures_publication.py
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
PUB = OUT / "figures_pub"
PUB.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="paper")
plt.rcParams.update({
    "font.size": 8.5, "axes.titlesize": 9.5, "axes.labelsize": 8.5,
    "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": 0.6, "savefig.dpi": 300, "savefig.bbox": "tight",
    "figure.dpi": 100,
})

COL = {"naive": "#8c8c8c", "mt": "#2E86AB", "rn": "#A23B72",
       "single": "#C44536", "fa1": "#DDA0DD", "fa2": "#8FBC8F", "fa3": "#F4A460"}
MODELS_A = [("naive", "Naive mean"), ("mt", "MT-GBLUP"), ("rn", "RN-GBLUP")]
MODELS_B = [("single", "Single-env GBLUP"), ("mt", "MT-GBLUP")]


def fig1_simulation():
    s = pd.read_csv(OUT / "02_sim_summary.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6), sharey=True)
    for ax, h2 in zip(axes, [0.2, 0.4, 0.6]):
        sub = s[s.h2 == h2]
        for m, lab in MODELS_A:
            ax.errorbar(sub.rG, sub[f"taskA_{m}_mean"], yerr=sub[f"taskA_{m}_se"],
                        marker="o", ms=3, lw=1, label=lab, color=COL[m], capsize=2)
        ax.set_title(f"h\u00b2 = {h2}")
        ax.set_xlabel("Genetic correlation rG")
        ax.set_ylim(0, 0.9)
        ax.grid(alpha=0.3, lw=0.4)
    axes[0].set_ylabel("Prediction accuracy (r)")
    axes[0].legend(frameon=True, loc="upper left")
    fig.suptitle("Predicting breeding values in an untested environment "
                 "(compound-symmetry G\u00d7E, 15 replicates)", y=1.02)
    fig.tight_layout()
    fig.savefig(PUB / "fig1_simulation.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6), sharey=True)
    for ax, h2 in zip(axes, [0.2, 0.4, 0.6]):
        sub = s[s.h2 == h2]
        for m, lab in MODELS_B:
            ax.errorbar(sub.rG, sub[f"taskB_{m}_mean"], yerr=sub[f"taskB_{m}_se"],
                        marker="o", ms=3, lw=1, label=lab, color=COL[m], capsize=2)
        ax.set_title(f"h\u00b2 = {h2}")
        ax.set_xlabel("Genetic correlation rG")
        ax.set_ylim(0, 0.9)
        ax.grid(alpha=0.3, lw=0.4)
    axes[0].set_ylabel("Prediction accuracy (r)")
    axes[0].legend(frameon=True, loc="upper left")
    fig.suptitle("Predicting breeding values of untested lines "
                 "(5-fold CV, 15 replicates)", y=1.02)
    fig.tight_layout()
    fig.savefig(PUB / "fig2_simulation_lines.png")
    plt.close(fig)


def fig2_real_data():
    real = json.load(open(OUT / "03_real_validation.json"))
    rB = np.array(real["typeB_correlation"])
    env = real["env_names"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw={"width_ratios": [1, 1.2]})
    sns.heatmap(rB, annot=True, fmt=".2f", cmap="RdBu_r", vmin=-0.1, vmax=1,
                xticklabels=env, yticklabels=env, square=True, ax=axes[0],
                cbar_kws={"label": "type-B rG"}, annot_kws={"size": 8})
    axes[0].set_title("Estimated genetic correlation")
    x = np.arange(3); w = 0.32
    axes[1].bar(x - w, [real["taskA"]["naive"], real["taskA"]["mt"], real["taskA"]["rn"]],
                width=w, color=[COL["naive"], COL["mt"], COL["rn"]], label="Task A: new env")
    axes[1].bar(x, [real["taskB"]["single"], real["taskB"]["mt"], np.nan],
                width=w, color=[COL["single"], COL["mt"], "#cccccc"], label="Task B: new lines")
    for xi, v in zip(x - w, [real["taskA"]["naive"], real["taskA"]["mt"], real["taskA"]["rn"]]):
        axes[1].text(xi, v + 0.012, f"{v:.2f}", ha="center", fontsize=7.5)
    for xi, v in zip(x, [real["taskB"]["single"], real["taskB"]["mt"], np.nan]):
        if not np.isnan(v):
            axes[1].text(xi, v + 0.012, f"{v:.2f}", ha="center", fontsize=7.5)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(["Naive /\nsingle-env", "MT-GBLUP", "RN-GBLUP"])
    axes[1].set_ylim(0, 0.55)
    axes[1].set_title("Real-data validation (CIMMYT wheat)")
    axes[1].legend(frameon=True, fontsize=7)
    axes[1].grid(alpha=0.3, lw=0.4)
    fig.suptitle("Real-data validation (Crossa et al. 2010 wheat panel)", y=1.02)
    fig.tight_layout()
    fig.savefig(PUB / "fig3_real_data.png")
    plt.close(fig)


def fig3_robustness():
    s4 = pd.read_csv(OUT / "04_structured_summary.csv")
    s5 = pd.read_csv(OUT / "05_unbalanced_summary.csv")
    s7 = pd.read_csv(OUT / "07_fa_summary.csv")

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6))

    # structured
    ax = axes[0]
    xx = np.arange(2); w = 0.26
    for i, (m, lab) in enumerate(MODELS_A):
        vals = [s4[(s4.kappa == k) & (s4.h2 == 0.6)][f"taskA_{m}_mean"].iloc[0] for k in [0.3, 0.8]]
        ax.bar(xx + (i - 1) * w, vals, width=w, label=lab, color=COL[m])
    ax.set_xticks(xx); ax.set_xticklabels(["\u03ba = 0.3", "\u03ba = 0.8"])
    ax.set_ylim(0, 1); ax.set_ylabel("Accuracy (r)")
    ax.set_title("Structured (reaction-norm) G\u00d7E, h\u00b2 = 0.6")
    ax.grid(alpha=0.3, lw=0.4)

    # FA rank
    ax = axes[1]
    order = ["naive", "fa1", "fa2", "fa3", "rn", "mt"]
    labels = ["Naive", "FA1", "FA2", "FA3", "RN", "MT"]
    sub = s7[s7.scenario == "kappa0.8_h2_0.6"].iloc[0]
    vals = [sub[f"{m}_mean"] for m in order]
    ax.bar(np.arange(6), vals, color=[COL[m] for m in order])
    ax.set_xticks(np.arange(6)); ax.set_xticklabels(labels)
    ax.set_ylim(0, 1)
    ax.set_title("FA rank (reaction-norm G\u00d7E, \u03ba = 0.8)\n"
                 f"varExp FA1/2/3 = {sub.var_exp_1:.2f}/{sub.var_exp_2:.2f}/{sub.var_exp_3:.2f}")
    ax.grid(alpha=0.3, lw=0.4)

    # unbalanced
    ax = axes[2]
    xx = np.arange(2); w = 0.3
    ax.bar(xx - w / 2, s5.taskA_naive_mean, width=w, label="Task A: naive", color=COL["naive"])
    ax.bar(xx + w / 2, s5.taskA_mt_mean, width=w, label="Task A: MT", color=COL["mt"])
    ax.bar(xx - w / 2 + 0.5, s5.taskB_single_mean, width=w, label="Task B: single", color=COL["single"])
    ax.bar(xx + w / 2 + 0.5, s5.taskB_mt_mean, width=w, label="Task B: MT", color=COL["mt"])
    ax.set_xticks(xx + 0.5); ax.set_xticklabels([f"rG = {r}" for r in s5.rg])
    ax.set_ylim(0, 0.8)
    ax.set_title("30% missing data (h\u00b2 = 0.4)")
    ax.legend(frameon=True, fontsize=6)
    ax.grid(alpha=0.3, lw=0.4)

    fig.suptitle("Robustness analyses", y=1.02)
    fig.tight_layout()
    fig.savefig(PUB / "fig4_robustness.png")
    plt.close(fig)


def fig4_summary():
    """Compact overview reproducing 06_summary at publication quality."""
    s = pd.read_csv(OUT / "02_sim_summary.csv")
    real = json.load(open(OUT / "03_real_validation.json"))
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.0))
    sub = s[s.h2 == 0.4]
    ax = axes[0, 0]
    for m, lab in MODELS_A:
        ax.errorbar(sub.rG, sub[f"taskA_{m}_mean"], yerr=sub[f"taskA_{m}_se"],
                    marker="o", ms=2.5, lw=1, label=lab, color=COL[m], capsize=1.5)
    ax.set_ylim(0, 0.85); ax.set_title("(a) Sim: new environment"); ax.set_ylabel("r")
    ax.legend(frameon=True, fontsize=6)

    ax = axes[0, 1]
    for m, lab in MODELS_B:
        ax.errorbar(sub.rG, sub[f"taskB_{m}_mean"], yerr=sub[f"taskB_{m}_se"],
                    marker="o", ms=2.5, lw=1, label=lab, color=COL[m], capsize=1.5)
    ax.set_ylim(0, 0.85); ax.set_title("(b) Sim: new lines")
    ax.legend(frameon=True, fontsize=6)

    ax = axes[0, 2]
    x = np.arange(3); w = 0.3
    ax.bar(x - w, [real["taskA"]["naive"], real["taskA"]["mt"], real["taskA"]["rn"]],
           width=w, color=[COL["naive"], COL["mt"], COL["rn"]], label="A")
    ax.bar(x, [real["taskB"]["single"], real["taskB"]["mt"], np.nan],
           width=w, color=[COL["single"], COL["mt"], "#cccccc"], label="B")
    ax.set_xticks(x); ax.set_xticklabels(["Naive/\nsingle", "MT", "RN"])
    ax.set_ylim(0, 0.6); ax.set_title("(c) Real wheat")
    ax.legend(frameon=True, fontsize=6)

    ax = axes[1, 0]
    s4 = pd.read_csv(OUT / "04_structured_summary.csv")
    xx = np.arange(2); w = 0.26
    for i, (m, lab) in enumerate(MODELS_A):
        vals = [s4[(s4.kappa == k) & (s4.h2 == 0.6)][f"taskA_{m}_mean"].iloc[0] for k in [0.3, 0.8]]
        ax.bar(xx + (i - 1) * w, vals, width=w, label=lab, color=COL[m])
    ax.set_xticks(xx); ax.set_xticklabels(["\u03ba=0.3", "\u03ba=0.8"])
    ax.set_ylim(0, 1); ax.set_title("(d) Structured G\u00d7E")
    ax.legend(frameon=True, fontsize=6)

    ax = axes[1, 1]
    s7 = pd.read_csv(OUT / "07_fa_summary.csv")
    sub7 = s7[s7.scenario == "kappa0.8_h2_0.6"].iloc[0]
    order = ["naive", "fa1", "fa2", "fa3", "rn", "mt"]
    ax.bar(np.arange(6), [sub7[f"{m}_mean"] for m in order], color=[COL[m] for m in order])
    ax.set_xticks(np.arange(6)); ax.set_xticklabels(["N", "FA1", "FA2", "FA3", "RN", "MT"], fontsize=6.5)
    ax.set_ylim(0, 1); ax.set_title("(e) FA rank (RN G\u00d7E)")

    ax = axes[1, 2]
    rB = np.array(real["typeB_correlation"])
    sns.heatmap(rB, annot=True, fmt=".2f", cmap="RdBu_r", vmin=-0.1, vmax=1,
                xticklabels=real["env_names"], yticklabels=real["env_names"],
                square=True, ax=ax, cbar_kws={"label": "rG"}, annot_kws={"size": 7})
    ax.set_title("(f) Wheat type-B rG")

    fig.suptitle("Multi-environment genomic prediction under G\u00d7E", y=1.0)
    fig.tight_layout()
    fig.savefig(PUB / "fig5_overview.png")
    plt.close(fig)


if __name__ == "__main__":
    fig1_simulation()
    fig2_real_data()
    fig3_robustness()
    fig4_summary()
    print("Saved publication figures to", PUB)
    for p in sorted(PUB.glob("*.png")):
        print("  ", p.name)
