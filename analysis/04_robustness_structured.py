"""
04_robustness_structured.py
===========================
Robustness extension: STRUCTURED (reaction-norm / covariate-driven) GxE.

In contrast to the main simulation (compound-symmetry / unstructured GxE), here
the true generative model is a linear reaction norm g_ij = c_i + b_i z_j, i.e.
the genetic covariance is rank-2 and driven by an observable environmental
covariate z. This is the regime where reduced-rank reaction-norm models are
expected to be competitive.

Crucially, all models here use variance components ESTIMATED from the data
(multivariate Haseman-Elston), matching what a practitioner would do; the
reaction-norm model constrains the estimated covariance to the span of {1, z}
(via projection). For reference, the oracle-component versions are also
computed.

Scenarios: kappa (slope/intercept variance ratio) in {0.3, 0.8} x h2 in {0.3, 0.6},
20 replicates each.

Outputs:
  output/04_structured_results.json, output/04_structured_summary.csv
  output/figures/04_structured.png
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
from simulator import simulate_structured
from models import (gblup_single, mt_gblup, mt_predict_new_env,
                    mt_predict_new_lines, rn_covariance)

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sns.set_theme(style="whitegrid", context="notebook")

KAPPAS = [0.3, 0.8]
H2S = [0.3, 0.6]
REPS = 20
E = 4
KFOLD = 5

res = pyreadr.read_r(str(ROOT / "data" / "wheat.RData"))
MAF_EMPIRICAL = res["wheat.X"].values.astype(float).mean(axis=0)


def run_one(kappa, h2, rep):
    d = simulate_structured(n=599, m=1279, e=E, h2=h2, kappa=kappa,
                            maf_empirical=MAF_EMPIRICAL, seed=rep)
    G = grm_vanraden(d["X"], coding="01")
    eig = np.linalg.eigh(G)
    Sig_g_hat, Sig_e_hat = estimate_components_haseman_elston(d["Y"], G)
    Sig_g_hat = make_psd(Sig_g_hat)          # covariance 'bending' (finite-sample fix)
    Sig_e_hat = np.diag(np.clip(np.diag(Sig_e_hat), 0.05, None))
    # observable covariate for the reaction norm = environment mean phenotype
    z_obs = d["Y"].mean(axis=0)
    Sig_rr_hat = rn_covariance(Sig_g_hat, z_obs)

    Y, Gt = d["Y"], d["G_true"]
    n = Y.shape[0]

    # ---- Task A (predict a new environment), estimated components ----
    accA = {"naive": [], "mt": [], "rn": [], "mt_oracle": [], "rn_oracle": []}
    Sig_g_true, Sig_e_true = d["Sig_g"], d["Sig_e"]
    Sig_rr_true = rn_covariance(Sig_g_true, d["z"])
    for j in range(E):
        tr = [k for k in range(E) if k != j]
        Ytr = Y[:, tr] - Y[:, tr].mean(axis=0)
        accA["naive"].append(np.corrcoef(Ytr.mean(axis=1), Gt[:, j])[0, 1])
        # estimated
        Uh = mt_gblup(Ytr, G, Sig_g_hat[np.ix_(tr, tr)], Sig_e_hat[np.ix_(tr, tr)], eig=eig)
        accA["mt"].append(np.corrcoef(mt_predict_new_env(Uh, Sig_g_hat, j, tr), Gt[:, j])[0, 1])
        Uhr = mt_gblup(Ytr, G, Sig_rr_hat[np.ix_(tr, tr)], Sig_e_hat[np.ix_(tr, tr)], eig=eig)
        accA["rn"].append(np.corrcoef(mt_predict_new_env(Uhr, Sig_rr_hat, j, tr), Gt[:, j])[0, 1])
        # oracle
        Uo = mt_gblup(Ytr, G, Sig_g_true[np.ix_(tr, tr)], Sig_e_true[np.ix_(tr, tr)], eig=eig)
        accA["mt_oracle"].append(np.corrcoef(mt_predict_new_env(Uo, Sig_g_true, j, tr), Gt[:, j])[0, 1])
        Uor = mt_gblup(Ytr, G, Sig_rr_true[np.ix_(tr, tr)], Sig_e_true[np.ix_(tr, tr)], eig=eig)
        accA["rn_oracle"].append(np.corrcoef(mt_predict_new_env(Uor, Sig_rr_true, j, tr), Gt[:, j])[0, 1])
    accA = {k: float(np.mean(v)) for k, v in accA.items()}

    # ---- Task B (predict new lines), estimated components, 5-fold ----
    rng = np.random.default_rng(rep + 999)
    idx = rng.permutation(n)
    folds = np.array_split(idx, KFOLD)
    accB = {"single": np.zeros(E), "mt": np.zeros(E), "rn": np.zeros(E)}
    for tr in folds:
        va = np.ones(n, bool)
        va[tr] = False
        G_train = G[np.ix_(va, va)]
        G_cross = G[np.ix_(tr, va)]
        Uh_mt = mt_gblup(Y[va] - Y[va].mean(axis=0), G_train, Sig_g_hat, Sig_e_hat)
        U_mt = mt_predict_new_lines(G_cross, G_train, Uh_mt)
        Uh_rn = mt_gblup(Y[va] - Y[va].mean(axis=0), G_train, Sig_rr_hat, Sig_e_hat)
        U_rn = mt_predict_new_lines(G_cross, G_train, Uh_rn)
        for j in range(E):
            yc = Y[:, j] - Y[:, j].mean()
            lam = Sig_e_hat[j, j] / max(Sig_g_hat[j, j], 1e-6)
            ps = gblup_single(G_train, yc[va], G_cross, lam)
            accB["single"][j] += np.corrcoef(ps, Gt[tr, j])[0, 1] / KFOLD
            accB["mt"][j] += np.corrcoef(U_mt[:, j], Gt[tr, j])[0, 1] / KFOLD
            accB["rn"][j] += np.corrcoef(U_rn[:, j], Gt[tr, j])[0, 1] / KFOLD
    accB = {k: float(np.mean(v)) for k, v in accB.items()}

    # realized genetic-correlation range for reporting
    C = np.corrcoef(Gt.T)
    off = [C[j, k] for j in range(E) for k in range(j + 1, E)]
    return {"rep": rep, "taskA": accA, "taskB": accB,
            "realized_rG_min": float(np.min(off)), "realized_rG_max": float(np.max(off))}


def main():
    t0 = time.time()
    results = {}
    for kappa in KAPPAS:
        for h2 in H2S:
            key = f"kappa{kappa}_h2{h2}"
            reps = [run_one(kappa, h2, rep) for rep in range(REPS)]
            results[key] = reps
            a = {m: np.mean([r["taskA"][m] for r in reps]) for m in reps[0]["taskA"]}
            b = {m: np.mean([r["taskB"][m] for r in reps]) for m in reps[0]["taskB"]}
            rmin = np.mean([r["realized_rG_min"] for r in reps])
            rmax = np.mean([r["realized_rG_max"] for r in reps])
            print(f"[{key}] rG∈[{rmin:.2f},{rmax:.2f}]  A: naive={a['naive']:.3f} mt={a['mt']:.3f} rn={a['rn']:.3f}"
                  f" (oracle mt={a['mt_oracle']:.3f} rn={a['rn_oracle']:.3f})"
                  f"  |  B: single={b['single']:.3f} mt={b['mt']:.3f} rn={b['rn']:.3f}"
                  f"  ({time.time()-t0:.0f}s)", flush=True)

    with open(OUT / "04_structured_results.json", "w") as f:
        json.dump({"KAPPAS": KAPPAS, "H2S": H2S, "REPS": REPS, "results": results}, f)

    rows = []
    for kappa in KAPPAS:
        for h2 in H2S:
            key = f"kappa{kappa}_h2{h2}"
            reps = results[key]
            row = {"kappa": kappa, "h2": h2}
            for task in ["taskA", "taskB"]:
                for m in reps[0][task]:
                    vals = [r[task][m] for r in reps]
                    row[f"{task}_{m}_mean"] = round(float(np.mean(vals)), 4)
                    row[f"{task}_{m}_se"] = round(float(np.std(vals) / np.sqrt(REPS)), 4)
            row["realized_rG_min"] = round(float(np.mean([r["realized_rG_min"] for r in reps])), 3)
            row["realized_rG_max"] = round(float(np.mean([r["realized_rG_max"] for r in reps])), 3)
            rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "04_structured_summary.csv", index=False)

    plot_results(pd.DataFrame(rows))
    print(f"\nDone in {time.time()-t0:.0f}s")


def plot_results(summary):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
    for ax, h2 in zip(axes, H2S):
        sub = summary[summary.h2 == h2]
        x = np.arange(len(KAPPAS))
        w = 0.28
        models = [("naive", "Naive"), ("mt", "MT-GBLUP"), ("rn", "RN-GBLUP")]
        colors = {"naive": "#999999", "mt": "#2E86AB", "rn": "#A23B72"}
        for i, (m, lab) in enumerate(models):
            ax.bar(x + (i - 1) * w, sub[f"taskA_{m}_mean"],
                   yerr=sub[f"taskA_{m}_se"], width=w, label=lab, color=colors[m], capsize=2)
        ax.set_xticks(x)
        ax.set_xticklabels([f"κ={k}" for k in KAPPAS])
        ax.set_title(f"h² = {h2}")
        ax.set_xlabel("Reaction-norm slope/intercept ratio κ")
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Prediction accuracy (r), new environment")
    axes[0].legend(frameon=True)
    fig.suptitle("Structured (reaction-norm) G×E — estimated variance components", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIG / "04_structured.png", dpi=160)
    plt.close("all")


if __name__ == "__main__":
    main()
