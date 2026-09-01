"""
07_robustness_factor_analytic.py
================================
Factor-analytic (FA) rank analysis: how much of the environment covariance
structure is needed for multi-environment prediction?

Motivation: the main study showed MT-GBLUP (full covariance) is best under
unstructured GxE, while the reaction-norm model suffices under low-rank
structured GxE. Here we quantify "how much rank do you need" by fitting the
family of FA(k) models (Burgueno et al. 2012) to the Haseman-Elston estimated
covariance, for k = 1..4 (k = 4 = full MT), and reporting:

  (i) the proportion of genetic variance explained by the leading k factors
      (the standard rank diagnostic), and
  (ii) Task A (new environment) prediction accuracy of FA(k) vs naive vs RN.

Scenarios: unstructured GxE (rG = 0.3 / 0.6, h2 = 0.4) and structured
reaction-norm GxE (kappa = 0.3 / 0.8, h2 = 0.6); 10 replicates each.

Outputs: output/07_fa_results.json, output/07_fa_summary.csv,
         output/figures/07_factor_analytic.png
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
from simulator import simulate, simulate_structured
from models import (mt_gblup, mt_predict_new_env, rn_covariance, fa_covariance)

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sns.set_theme(style="whitegrid", context="notebook")

REPS = 10
E = 4
res = pyreadr.read_r(str(ROOT / "data" / "wheat.RData"))
MAF_EMPIRICAL = res["wheat.X"].values.astype(float).mean(axis=0)


def task_a_models(d, G, eig, Sig_g, Sig_e, z_obs):
    Y, Gt = d["Y"], d["G_true"]
    e = Y.shape[1]
    Sig_rr = rn_covariance(Sig_g, z_obs)
    out = {"naive": [], "mt": [], "rn": [],
           "fa1": [], "fa2": [], "fa3": []}
    for j in range(e):
        tr = [k for k in range(e) if k != j]
        Ytr = Y[:, tr] - Y[:, tr].mean(axis=0)
        out["naive"].append(np.corrcoef(Ytr.mean(axis=1), Gt[:, j])[0, 1])
        for name, Sig in [("mt", Sig_g), ("rn", Sig_rr),
                          ("fa1", fa_covariance(Sig_g, 1)),
                          ("fa2", fa_covariance(Sig_g, 2)),
                          ("fa3", fa_covariance(Sig_g, 3))]:
            Uh = mt_gblup(Ytr, G, Sig[np.ix_(tr, tr)], Sig_e[np.ix_(tr, tr)], eig=eig)
            out[name].append(np.corrcoef(mt_predict_new_env(Uh, Sig, j, tr), Gt[:, j])[0, 1])
    return {k: float(np.mean(v)) for k, v in out.items()}


def factor_explained(Sig_g):
    w = np.linalg.eigvalsh((Sig_g + Sig_g.T) / 2.0)
    w = np.clip(w, 0.0, None)
    total = w.sum()
    if total <= 0:
        return [0.0, 0.0, 0.0]
    order = np.sort(w)[::-1]
    return [float(order[:k].sum() / total) for k in (1, 2, 3)]


def run_one(kind, p1, p2, rep):
    if kind == "unstructured":
        d = simulate(n=599, m=1279, e=E, h2=p2, rg=p1, maf_empirical=MAF_EMPIRICAL, seed=rep)
    else:
        d = simulate_structured(n=599, m=1279, e=E, h2=p2, kappa=p1,
                                maf_empirical=MAF_EMPIRICAL, seed=rep)
    G = grm_vanraden(d["X"], coding="01")
    eig = np.linalg.eigh(G)
    Sig_g, Sig_e = estimate_components_haseman_elston(d["Y"], G)
    Sig_g = make_psd(Sig_g)
    Sig_e = np.diag(np.clip(np.diag(Sig_e), 0.05, None))
    z_obs = d["Y"].mean(axis=0)
    acc = task_a_models(d, G, eig, Sig_g, Sig_e, z_obs)
    acc["var_exp"] = factor_explained(Sig_g)
    return acc


def main():
    t0 = time.time()
    results = {}
    scenarios = [
        ("unstructured", "rg0.3_h2_0.4", 0.3, 0.4),
        ("unstructured", "rg0.6_h2_0.4", 0.6, 0.4),
        ("structured", "kappa0.3_h2_0.6", 0.3, 0.6),
        ("structured", "kappa0.8_h2_0.6", 0.8, 0.6),
    ]
    for kind, key, p1, p2 in scenarios:
        reps = [run_one(kind, p1, p2, rep) for rep in range(REPS)]
        results[key] = reps
        m = {k: np.mean([r[k] for r in reps]) for k in reps[0] if k != "var_exp"}
        ve = np.mean([r["var_exp"] for r in reps], axis=0)
        print(f"[{key:18s}] varExp(1,2,3)={ve[0]:.2f},{ve[1]:.2f},{ve[2]:.2f}  "
              f"acc: naive={m['naive']:.3f} mt={m['mt']:.3f} rn={m['rn']:.3f} "
              f"fa1={m['fa1']:.3f} fa2={m['fa2']:.3f} fa3={m['fa3']:.3f}  ({time.time()-t0:.0f}s)",
              flush=True)

    with open(OUT / "07_fa_results.json", "w") as f:
        json.dump({"REPS": REPS, "results": results}, f)

    rows = []
    for key, reps in results.items():
        row = {"scenario": key}
        for m in ["naive", "mt", "rn", "fa1", "fa2", "fa3"]:
            vals = [r[m] for r in reps]
            row[f"{m}_mean"] = round(float(np.mean(vals)), 4)
            row[f"{m}_se"] = round(float(np.std(vals) / np.sqrt(REPS)), 4)
        ve = np.mean([r["var_exp"] for r in reps], axis=0)
        for k in (1, 2, 3):
            row[f"var_exp_{k}"] = round(float(ve[k - 1]), 3)
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "07_fa_summary.csv", index=False)

    plot_fa(pd.DataFrame(rows))
    print(f"\nDone in {time.time()-t0:.0f}s")


def plot_fa(s):
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.0), sharey=True)
    labels = {"naive": "Naive", "mt": "MT", "rn": "RN", "fa1": "FA1", "fa2": "FA2", "fa3": "FA3"}
    colors = {"naive": "#999999", "mt": "#2E86AB", "rn": "#A23B72",
              "fa1": "#DDA0DD", "fa2": "#8FBC8F", "fa3": "#F4A460"}
    order = ["naive", "fa1", "fa2", "fa3", "rn", "mt"]
    for ax, key in zip(axes, s.scenario):
        sub = s[s.scenario == key].iloc[0]
        xs = np.arange(len(order))
        vals = [sub[f"{m}_mean"] for m in order]
        errs = [sub[f"{m}_se"] for m in order]
        ax.bar(xs, vals, yerr=errs, color=[colors[m] for m in order], capsize=2)
        ax.set_xticks(xs)
        ax.set_xticklabels([labels[m] for m in order], rotation=45, fontsize=8)
        ax.set_ylim(0, 1)
        ax.set_title(f"{key}\n(varExp FA1/2/3 = {sub.var_exp_1:.2f}/{sub.var_exp_2:.2f}/{sub.var_exp_3:.2f})",
                     fontsize=9)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Prediction accuracy (r), new environment")
    fig.suptitle("Factor-analytic rank analysis of the estimated environment covariance",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(FIG / "07_factor_analytic.png", dpi=160)
    plt.close("all")


if __name__ == "__main__":
    main()
