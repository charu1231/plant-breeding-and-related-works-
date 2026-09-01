"""
05_robustness_unbalanced.py
===========================
Preliminary robustness check: does the model ranking hold under MISSING /
unbalanced data? (Breeding multi-environment trials are frequently unbalanced.)

Design:
  - Use the unstructured (compound-symmetry) simulator at rG in {0.3, 0.7},
    h2 = 0.4, 20 replicates.
  - Mask 30% of the phenotype cells completely at random (MCAR).
  - Compare the SAME tasks as the main study, under a simple, transparent
    strategy: single-environment GBLUP and the naive mean use whatever cells
    are observed (they handle missingness natively), while MT-GBLUP is fitted
    on environment-mean imputed data with Haseman-Elston estimated components.

Interpretation: this is a screening check (not a full missing-data methodology);
it asks whether the qualitative conclusions (MT >= single-env; MT > naive) are
upset by 30% MCAR missingness under simple imputation.

Outputs:
  output/05_unbalanced_results.json, output/05_unbalanced_summary.csv
  output/figures/05_unbalanced.png
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
from utils import grm_vanraden, estimate_components_haseman_elston, make_psd
from simulator import simulate
from models import (gblup_single, mt_gblup, mt_predict_new_env,
                    mt_predict_new_lines)

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sns.set_theme(style="whitegrid", context="notebook")

RGS = [0.3, 0.7]
H2 = 0.4
REPS = 20
MISS = 0.30
E = 4
KFOLD = 5

res = pyreadr.read_r(str(ROOT / "data" / "wheat.RData"))
MAF_EMPIRICAL = res["wheat.X"].values.astype(float).mean(axis=0)


def run_one(rg, rep):
    d = simulate(n=599, m=1279, e=E, h2=H2, rg=rg, maf_empirical=MAF_EMPIRICAL, seed=rep)
    G = grm_vanraden(d["X"], coding="01")
    eig = np.linalg.eigh(G)
    Y, Gt = d["Y"].copy(), d["G_true"]
    n = Y.shape[0]

    # ---- mask 30% of cells MCAR ----
    rng = np.random.default_rng(rep + 777)
    mask = rng.random((n, E)) < MISS
    Y_obs = np.where(mask, np.nan, Y)
    # environment-mean imputation (column means of observed cells)
    Y_imp = Y_obs.copy()
    for j in range(E):
        col = Y_obs[:, j]
        mu_j = np.nanmean(col)
        Y_imp[np.isnan(col), j] = mu_j

    Sig_g_hat, Sig_e_hat = estimate_components_haseman_elston(Y_imp, G)
    Sig_g_hat = make_psd(Sig_g_hat)
    Sig_e_hat = np.diag(np.clip(np.diag(Sig_e_hat), 0.05, None))

    # ---- Task A (predict a new environment) ----
    accA = {"naive": [], "mt": []}
    for j in range(E):
        tr = [k for k in range(E) if k != j]
        # naive: row mean of OBSERVED cells in the training environments
        row_means = np.nanmean(Y_obs[:, tr], axis=1)
        valid = ~np.isnan(row_means)
        accA["naive"].append(np.corrcoef(row_means[valid], Gt[valid, j])[0, 1])
        # MT on imputed training environments
        Ytr = Y_imp[:, tr] - Y_imp[:, tr].mean(axis=0)
        Uh = mt_gblup(Ytr, G, Sig_g_hat[np.ix_(tr, tr)], Sig_e_hat[np.ix_(tr, tr)], eig=eig)
        mt = mt_predict_new_env(Uh, Sig_g_hat, j, tr)
        accA["mt"].append(np.corrcoef(mt, Gt[:, j])[0, 1])
    accA = {k: float(np.mean(v)) for k, v in accA.items()}

    # ---- Task B (predict new lines), 5-fold ----
    idx = rng.permutation(n)
    folds = np.array_split(idx, KFOLD)
    accB = {"single": np.zeros(E), "mt": np.zeros(E)}
    for tr in folds:
        va = np.ones(n, bool)
        va[tr] = False
        G_train = G[np.ix_(va, va)]
        G_cross = G[np.ix_(tr, va)]
        Uh = mt_gblup(Y_imp[va] - Y_imp[va].mean(axis=0), G_train, Sig_g_hat, Sig_e_hat)
        U_te = mt_predict_new_lines(G_cross, G_train, Uh)
        for j in range(E):
            # single-env on OBSERVED cells only
            obs = ~np.isnan(Y_obs[:, j])
            yc = Y_obs[:, j] - np.nanmean(Y_obs[:, j])
            yc = np.where(np.isnan(yc), 0.0, yc)
            lam = Sig_e_hat[j, j] / max(Sig_g_hat[j, j], 1e-6)
            tr_obs = va & obs
            ps = gblup_single(G[np.ix_(tr_obs, tr_obs)], yc[tr_obs],
                              G[np.ix_(tr, tr_obs)], lam)
            accB["single"][j] += np.corrcoef(ps, Gt[tr, j])[0, 1] / KFOLD
            accB["mt"][j] += np.corrcoef(U_te[:, j], Gt[tr, j])[0, 1] / KFOLD
    accB = {k: float(np.mean(v)) for k, v in accB.items()}
    return {"rep": rep, "taskA": accA, "taskB": accB}


def main():
    t0 = time.time()
    results = {}
    for rg in RGS:
        key = f"rg{rg}_h2{H2}"
        reps = [run_one(rg, rep) for rep in range(REPS)]
        results[key] = reps
        a = {m: np.mean([r["taskA"][m] for r in reps]) for m in reps[0]["taskA"]}
        b = {m: np.mean([r["taskB"][m] for r in reps]) for m in reps[0]["taskB"]}
        print(f"[{key}] A: naive={a['naive']:.3f} mt={a['mt']:.3f}"
              f"  |  B: single={b['single']:.3f} mt={b['mt']:.3f}  ({time.time()-t0:.0f}s)",
              flush=True)

    with open(OUT / "05_unbalanced_results.json", "w") as f:
        json.dump({"RGS": RGS, "H2": H2, "MISS": MISS, "REPS": REPS, "results": results}, f)

    rows = []
    for rg in RGS:
        key = f"rg{rg}_h2{H2}"
        reps = results[key]
        row = {"rg": rg, "miss": MISS}
        for task in ["taskA", "taskB"]:
            for m in reps[0][task]:
                vals = [r[task][m] for r in reps]
                row[f"{task}_{m}_mean"] = round(float(np.mean(vals)), 4)
                row[f"{task}_{m}_se"] = round(float(np.std(vals) / np.sqrt(REPS)), 4)
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "05_unbalanced_summary.csv", index=False)

    # ---- figure ----
    s = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0), sharey=True)
    x = np.arange(len(RGS)); w = 0.35
    for ax, task, models, colors, labels in [
        (axes[0], "taskA", ["naive", "mt"], ["#999999", "#2E86AB"], ["Naive mean", "MT-GBLUP (imputed)"]),
        (axes[1], "taskB", ["single", "mt"], ["#C44536", "#2E86AB"], ["Single-env GBLUP", "MT-GBLUP (imputed)"]),
    ]:
        for i, (m, c, lab) in enumerate(zip(models, colors, labels)):
            ax.bar(x + (i - 0.5) * w, s[f"{task}_{m}_mean"], yerr=s[f"{task}_{m}_se"],
                   width=w, color=c, label=lab, capsize=3)
        ax.set_xticks(x); ax.set_xticklabels([f"rG={rg}" for rg in RGS])
        ax.set_ylim(0, 1); ax.grid(alpha=0.3); ax.legend(frameon=True)
        ax.set_xlabel("Genetic correlation rG")
    axes[0].set_ylabel("Prediction accuracy (r)")
    axes[0].set_title("New environment (CV1)")
    axes[1].set_title("New lines (CV2)")
    fig.suptitle(f"30% random missing data (h² = {H2}) — ranking check", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIG / "05_unbalanced.png", dpi=160)
    plt.close("all")
    print(f"\nDone in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
