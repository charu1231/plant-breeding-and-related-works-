"""
03_real_data_validation.py
==========================
Phase 3: apply the benchmark models to the REAL CIMMYT wheat data
(599 lines x 1279 markers x 4 environments), using variance components
ESTIMATED from the data (multivariate Haseman-Elston / GCTA-style moments).

Models (mirror of the simulation study):
  Task A - predict a NEW environment (leave-one-env-out, "CV1"):
      naive   : phenotypic mean across training environments
      MT-GBLUP: multi-trait GBLUP with estimated Sig_g / Sig_e
      RN-GBLUP: rank-2 reaction norm on environment mean (estimated Sig)
  Task B - predict NEW lines (10-fold CV using the wheat.sets partition, "CV2"):
      single  : single-environment GBLUP (per-env lambda = se2/sg2)
      MT-GBLUP: multi-trait GBLUP

Outputs:
  output/03_real_validation.json
  output/figures/03_real_validation.png
"""
import json
import sys
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
from utils import (grm_vanraden, estimate_components_haseman_elston,
                   two_way_anova_components, make_psd)
from models import (gblup_single, mt_gblup, mt_predict_new_env,
                    mt_predict_new_lines, rn_covariance)

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sns.set_theme(style="whitegrid", context="notebook")

# ---------------------------------------------------------------- load
res = pyreadr.read_r(str(ROOT / "data" / "wheat.RData"))
Y = res["wheat.Y"].values.astype(float)
X = res["wheat.X"].values.astype(float)
sets = res["wheat.sets"].values.ravel().astype(int)
env_names = list(res["wheat.Y"].columns)
n, m = X.shape
e = Y.shape[1]

G = grm_vanraden(X, coding="01")

# ---------------------------------------------------------------- estimate
Sig_g, Sig_e = estimate_components_haseman_elston(Y, G)
Sig_g = make_psd(Sig_g)                       # covariance bending (finite sample)
Sig_e = np.diag(np.clip(np.diag(Sig_e), 0.05, None))
rB = Sig_g / np.sqrt(np.outer(np.diag(Sig_g), np.diag(Sig_g)))
print("Estimated genetic covariance Sig_g:")
print(pd.DataFrame(Sig_g, index=env_names, columns=env_names).round(3))
print("\nEstimated type-B genetic correlation (rG):")
rB_df = pd.DataFrame(rB, index=env_names, columns=env_names)
print(rB_df.round(2))
print("\nEstimated residual variances (diag Sig_e):", np.diag(Sig_e).round(3))
vc = two_way_anova_components(Y)
print("Two-way H2 (entry means):", round(vc["H2"], 3))

# ---------------------------------------------------------------- Task A (CV1)
accA = {"naive": [], "mt": [], "rn": []}
z = Y.mean(axis=0)                     # observable env covariate = env mean yield
Sig_rr = rn_covariance(Sig_g, z)
eig = np.linalg.eigh(G)
for j in range(e):
    tr = [k for k in range(e) if k != j]
    Ytr = Y[:, tr] - Y[:, tr].mean(axis=0)
    target = Y[:, j]                   # realised phenotype as proxy for BV (r=1)
    accA["naive"].append(np.corrcoef(Ytr.mean(axis=1), target)[0, 1])
    Uh = mt_gblup(Ytr, G, Sig_g[np.ix_(tr, tr)], Sig_e[np.ix_(tr, tr)], eig=eig)
    accA["mt"].append(np.corrcoef(mt_predict_new_env(Uh, Sig_g, j, tr), target)[0, 1])
    Uhr = mt_gblup(Ytr, G, Sig_rr[np.ix_(tr, tr)], Sig_e[np.ix_(tr, tr)], eig=eig)
    accA["rn"].append(np.corrcoef(mt_predict_new_env(Uhr, Sig_rr, j, tr), target)[0, 1])
accA = {k: float(np.mean(v)) for k, v in accA.items()}
print("\nTask A (predict new env, realised phenotype target):", accA)

# ---------------------------------------------------------------- Task B (CV2)
folds = sorted(set(sets.tolist()))
accB = {"single": np.zeros(e), "mt": np.zeros(e)}
for f in folds:
    tr = np.flatnonzero(sets != f)
    va = np.flatnonzero(sets == f)
    G_train = G[np.ix_(tr, tr)]
    G_cross = G[np.ix_(va, tr)]
    Uh = mt_gblup(Y[tr] - Y[tr].mean(axis=0), G_train, Sig_g, Sig_e)
    U_te = mt_predict_new_lines(G_cross, G_train, Uh)
    for j in range(e):
        yc = Y[:, j] - Y[:, j].mean()
        lam = Sig_e[j, j] / max(Sig_g[j, j], 1e-6)
        ps = gblup_single(G_train, yc[tr], G_cross, lam)
        accB["single"][j] += np.corrcoef(ps, Y[va, j])[0, 1] / len(folds)
        accB["mt"][j] += np.corrcoef(U_te[:, j], Y[va, j])[0, 1] / len(folds)
accB_mean = {"single": float(np.mean(accB["single"])),
             "mt": float(np.mean(accB["mt"])),
             "single_by_env": [round(float(x), 4) for x in accB["single"]],
             "mt_by_env": [round(float(x), 4) for x in accB["mt"]]}
print("Task B (predict new lines, 10-fold):", accB_mean)

# ---------------------------------------------------------------- save
out = {
    "dataset": "CIMMYT wheat (Crossa et al. 2010, Genetics 186:713-724)",
    "Sig_g": Sig_g.tolist(),
    "Sig_e": Sig_e.tolist(),
    "typeB_correlation": rB.round(3).tolist(),
    "H2_entry_means": round(vc["H2"], 3),
    "taskA": {k: round(v, 4) for k, v in accA.items()},
    "taskB": accB_mean,
    "env_names": env_names,
}
with open(OUT / "03_real_validation.json", "w") as f:
    json.dump(out, f, indent=2)

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
sns.heatmap(rB_df, annot=True, fmt=".2f", cmap="RdBu_r", vmin=-1, vmax=1,
            square=True, ax=axes[0], cbar_kws={"label": "rG"})
axes[0].set_title("Estimated type-B genetic correlation")

labels = ["Naive\nmean", "MT-GBLUP", "RN-GBLUP"]
valsA = [accA["naive"], accA["mt"], accA["rn"]]
x = np.arange(len(labels))
axes[1].bar(x - 0.2, valsA, width=0.4, label="Task A: new env", color="#2E86AB")
valsB = [accB_mean["single"], accB_mean["mt"], np.nan]
axes[1].bar(x + 0.2, [accB_mean["single"], accB_mean["mt"], np.nan],
            width=0.4, label="Task B: new lines", color="#C44536")
for xi, vi in zip(x - 0.2, valsA):
    axes[1].text(xi, vi + 0.01, f"{vi:.3f}", ha="center", fontsize=9)
for xi, vi in zip(x + 0.2, [accB_mean["single"], accB_mean["mt"], np.nan]):
    if not np.isnan(vi):
        axes[1].text(xi, vi + 0.01, f"{vi:.3f}", ha="center", fontsize=9)
axes[1].set_xticks(x)
axes[1].set_xticklabels(["Naive / single-env", "MT-GBLUP", "RN"])
axes[1].set_ylim(0, 1)
axes[1].set_ylabel("Prediction accuracy (r)")
axes[1].legend()
axes[1].set_title("Real-data validation (CIMMYT wheat)")
fig.tight_layout()
fig.savefig(FIG / "03_real_validation.png", dpi=160)
print("\nSaved:", OUT / "03_real_validation.json")
print("Saved:", FIG / "03_real_validation.png")
