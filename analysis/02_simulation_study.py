"""
02_simulation_study.py
======================
Phase 2: multi-environment genomic-prediction benchmark under GxE.

Simulation design (calibrated to the CIMMYT wheat panel):
  - n=599 related inbred (DH) lines, m=1279 markers (0/1), 4 environments
  - marker allele frequencies resampled from the real wheat panel
  - GxE parameterized by the genetic correlation rG (compound-symmetry factor model)
  - scenario grid: rG in {0.2, 0.4, 0.6, 0.8}  x  h2 in {0.2, 0.4, 0.6}
  - 15 independent replicates per scenario

Models benchmarked (prediction target = true breeding value g_ij):
  Task A - predict a NEW environment (leave-one-env-out):
      naive   : phenotypic mean across training environments (baseline)
      MT-GBLUP: environments as correlated traits, oracle (rG, h2)
      RN-GBLUP: rank-2 linear reaction norm on the environment mean (z)
  Task B - predict NEW lines (5-fold CV):
      single  : single-environment GBLUP (per env, averaged)
      MT-GBLUP: multi-trait GBLUP (per env, averaged)

Metric: prediction accuracy = Pearson correlation(predicted BV, true BV).

Outputs:
  output/02_sim_results.json  (per-scenario means and SEs, per-replicate raw)
  output/figures/02_sim_*.png
"""
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import grm_vanraden
from simulator import simulate
from models import (gblup_single, mt_gblup, mt_predict_new_env,
                    mt_predict_new_lines, rn_covariance)

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="notebook")

RGS = [0.2, 0.4, 0.6, 0.8]
H2S = [0.2, 0.4, 0.6]
REPS = 15
E = 4
KFOLD = 5

# ---------------------------------------------------------------- calibration
res = pyreadr.read_r(str(ROOT / "data" / "wheat.RData"))
Xw = res["wheat.X"].values.astype(float)
MAF_EMPIRICAL = Xw.mean(axis=0)          # 0/1 markers -> mean = allele frequency


def task_a(d, G, eig):
    """Predict breeding values in each environment from the other three."""
    Y, Gt = d["Y"], d["G_true"]
    Sig_g, Sig_e = d["Sig_g"], d["Sig_e"]
    n, e = Y.shape
    acc = {"naive": [], "mt": [], "rn": []}
    z = d["mu"]                            # observable env covariate (env mean)
    Sig_rr = rn_covariance(Sig_g, z)
    for j in range(e):
        tr = [k for k in range(e) if k != j]
        Ytr = Y[:, tr] - Y[:, tr].mean(axis=0)
        # naive
        naive = Ytr.mean(axis=1)
        acc["naive"].append(np.corrcoef(naive, Gt[:, j])[0, 1])
        # MT
        Uh = mt_gblup(Ytr, G, Sig_g[np.ix_(tr, tr)], Sig_e[np.ix_(tr, tr)], eig=eig)
        mt = mt_predict_new_env(Uh, Sig_g, j, tr)
        acc["mt"].append(np.corrcoef(mt, Gt[:, j])[0, 1])
        # RN
        Uhr = mt_gblup(Ytr, G, Sig_rr[np.ix_(tr, tr)], Sig_e[np.ix_(tr, tr)], eig=eig)
        rn = mt_predict_new_env(Uhr, Sig_rr, j, tr)
        acc["rn"].append(np.corrcoef(rn, Gt[:, j])[0, 1])
    return {k: float(np.mean(v)) for k, v in acc.items()}


def task_b(d, G, seed):
    """Predict breeding values of new lines by 5-fold CV, per environment."""
    Y, Gt = d["Y"], d["G_true"]
    Sig_g, Sig_e = d["Sig_g"], d["Sig_e"]
    n, e = Y.shape
    rng = np.random.default_rng(seed + 999)
    idx = rng.permutation(n)
    folds = np.array_split(idx, KFOLD)
    acc_single = np.zeros(e)
    acc_mt = np.zeros(e)
    for tr in folds:
        va = np.ones(n, bool)
        va[tr] = False
        G_train = G[np.ix_(va, va)]
        G_cross = G[np.ix_(tr, va)]
        # MT on training lines only
        Uh = mt_gblup(Y[va] - Y[va].mean(axis=0), G_train, Sig_g, Sig_e)
        U_te = mt_predict_new_lines(G_cross, G_train, Uh)
        for j in range(e):
            yc = Y[:, j] - Y[:, j].mean()
            ps = gblup_single(G_train, yc[va], G_cross, d["sigma2_e"])
            acc_single[j] += np.corrcoef(ps, Gt[tr, j])[0, 1] / KFOLD
            acc_mt[j] += np.corrcoef(U_te[:, j], Gt[tr, j])[0, 1] / KFOLD
    return {
        "single": float(np.mean(acc_single)),
        "mt": float(np.mean(acc_mt)),
        "single_by_env": [float(x) for x in acc_single],
        "mt_by_env": [float(x) for x in acc_mt],
    }


def run_one(rg, h2, rep):
    d = simulate(n=599, m=1279, e=E, h2=h2, rg=rg, maf_empirical=MAF_EMPIRICAL,
                 nf=100, ngen=15, C=21, seed=rep)
    G = grm_vanraden(d["X"], coding="01")
    eig = np.linalg.eigh(G)
    a = task_a(d, G, eig)
    b = task_b(d, G, rep)
    return {"rep": rep, "taskA": a, "taskB": b}


def main():
    t0 = time.time()
    results = {}
    for rg in RGS:
        for h2 in H2S:
            key = f"rg{rg}_h2{h2}"
            reps = [run_one(rg, h2, rep) for rep in range(REPS)]
            results[key] = reps
            accA = {m: np.mean([r["taskA"][m] for r in reps]) for m in ["naive", "mt", "rn"]}
            accB = {m: np.mean([r["taskB"][m] for r in reps]) for m in ["single", "mt"]}
            print(f"[{key}] A: naive={accA['naive']:.3f} mt={accA['mt']:.3f} rn={accA['rn']:.3f}"
                  f"  |  B: single={accB['single']:.3f} mt={accB['mt']:.3f}"
                  f"  (elapsed {time.time()-t0:.0f}s)", flush=True)

    with open(OUT / "02_sim_results.json", "w") as f:
        json.dump({"RGS": RGS, "H2S": H2S, "REPS": REPS, "results": results}, f)

    # ------------------------------------------------------------ summary table
    rows = []
    for rg in RGS:
        for h2 in H2S:
            key = f"rg{rg}_h2{h2}"
            reps = results[key]
            row = {"rG": rg, "h2": h2}
            for task, models in [("taskA", ["naive", "mt", "rn"]),
                                 ("taskB", ["single", "mt"])]:
                for m in models:
                    vals = [r[task][m] for r in reps]
                    row[f"{task}_{m}_mean"] = round(float(np.mean(vals)), 4)
                    row[f"{task}_{m}_se"] = round(float(np.std(vals) / np.sqrt(REPS)), 4)
            rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "02_sim_summary.csv", index=False)

    # ------------------------------------------------------------ figures
    plot_results(summary, results)
    print(f"\nDone in {time.time()-t0:.0f}s")


def plot_results(summary, results):
    # Task A: line plots accuracy vs rG, one panel per h2
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
    colors = {"naive": "#999999", "mt": "#2E86AB", "rn": "#A23B72"}
    for ax, h2 in zip(axes, H2S):
        sub = summary[summary.h2 == h2]
        for m in ["naive", "mt", "rn"]:
            ax.errorbar(sub.rG, sub[f"taskA_{m}_mean"], yerr=sub[f"taskA_{m}_se"],
                        marker="o", label=m.upper() if m != "naive" else "Naive mean",
                        color=colors[m], capsize=3)
        ax.set_title(f"h² = {h2}")
        ax.set_xlabel("Genetic correlation rG (G×E strength ↓)")
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Prediction accuracy (r)")
    axes[0].legend(title="Task A: new environment", frameon=True)
    fig.suptitle("Predicting breeding values in an untested environment", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIG / "02_sim_taskA.png", dpi=160)

    # Task B
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
    colorsB = {"single": "#C44536", "mt": "#2E86AB"}
    for ax, h2 in zip(axes, H2S):
        sub = summary[summary.h2 == h2]
        for m in ["single", "mt"]:
            ax.errorbar(sub.rG, sub[f"taskB_{m}_mean"], yerr=sub[f"taskB_{m}_se"],
                        marker="o", label="Single-env GBLUP" if m == "single" else "MT-GBLUP",
                        color=colorsB[m], capsize=3)
        ax.set_title(f"h² = {h2}")
        ax.set_xlabel("Genetic correlation rG (G×E strength ↓)")
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Prediction accuracy (r)")
    axes[0].legend(title="Task B: new lines", frameon=True)
    fig.suptitle("Predicting breeding values of untested lines (5-fold CV)", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIG / "02_sim_taskB.png", dpi=160)
    plt.close("all")


if __name__ == "__main__":
    main()
