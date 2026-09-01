"""
01_real_data_baseline.py
========================
Baseline analyses on the CIMMYT wheat dataset (BGLR `wheat`; Crossa et al. 2010,
Genetics 186:713-724):
  599 historical CIMMYT wheat lines, 1279 DArT markers, grain yield in 4 environments.

Outputs (saved to output/):
  1. environment correlation matrix + GxE variance components (JSON + heatmap PNG)
  2. GBLUP (VanRaden method-1) 5-fold CV prediction accuracy per environment (JSON + bar PNG)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import grm_vanraden, gblup_cv, two_way_anova_components

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "wheat.RData"
OUT = ROOT / "output"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="notebook")

# ---------------------------------------------------------------- load
res = pyreadr.read_r(str(DATA))
Y = res["wheat.Y"].values.astype(float)          # 599 x 4
X = res["wheat.X"].values.astype(float)          # 599 x 1279
A = res["wheat.A"].values.astype(float)
sets = res["wheat.sets"].values.ravel().astype(int)
env_names = list(res["wheat.Y"].columns)

n, m = X.shape
print(f"Lines={n}, markers={m}, environments={len(env_names)}")

# ---------------------------------------------------------------- 1. GxE structure
vc = two_way_anova_components(Y)
env_corr = np.corrcoef(Y.T)
env_corr_df = pd.DataFrame(env_corr, index=env_names, columns=env_names)

G = grm_vanraden(X, coding="01")
np.fill_diagonal(G, np.diag(G))  # keep as-is
print("\n=== Variance components (two-way, r=1) ===")
for k, v in vc.items():
    print(f"  {k:12s} {v:10.4f}")
print("\nEnvironment correlation matrix:")
print(env_corr_df.round(2))

# ---------------------------------------------------------------- 2. GBLUP CV per env
acc = {}
for j, env in enumerate(env_names):
    r = gblup_cv(G, Y[:, j], k=5, seed=2026)
    acc[env] = {"corr": round(float(r["corr"]), 4),
                "mean_best_lambda": round(float(np.mean(r["best_lams"])), 3)}
    print(f"  env {env}: 5-fold CV accuracy (r) = {acc[env]['corr']:.4f}")

# ---------------------------------------------------------------- save results
summary = {
    "dataset": "CIMMYT wheat (BGLR 'wheat', Crossa et al. 2010 Genetics 186:713-724)",
    "n_lines": n, "n_markers": m, "n_env": len(env_names),
    "variance_components": {k: round(float(v), 6) for k, v in vc.items()},
    "env_correlation": env_corr_df.round(4).to_dict(),
    "gblup_cv_accuracy": acc,
    "overall_mean_accuracy": round(float(np.mean([a["corr"] for a in acc.values()])), 4),
}
with open(OUT / "01_baseline_results.json", "w") as f:
    json.dump(summary, f, indent=2)

# ---------------------------------------------------------------- figures
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

sns.heatmap(env_corr_df, annot=True, fmt=".2f", cmap="RdBu_r",
            vmin=-1, vmax=1, square=True, cbar_kws={"label": "Pearson r"},
            ax=axes[0])
axes[0].set_title("Grain yield: environment correlation\n(low r = strong G×E)")

envs = list(acc.keys())
vals = [acc[e]["corr"] for e in envs]
bars = axes[1].bar(envs, vals, color="#4C72B0")
axes[1].axhline(summary["overall_mean_accuracy"], ls="--", color="grey",
                label=f"mean = {summary['overall_mean_accuracy']:.3f}")
axes[1].set_ylim(0, 1)
axes[1].set_xlabel("Environment")
axes[1].set_ylabel("Prediction accuracy (r), 5-fold CV")
axes[1].set_title("GBLUP per-environment prediction")
for b, v in zip(bars, vals):
    axes[1].text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.3f}",
                 ha="center", fontsize=9)
axes[1].legend()

fig.tight_layout()
fig.savefig(FIG / "01_baseline.png", dpi=160)
print("\nSaved:", OUT / "01_baseline_results.json")
print("Saved:", FIG / "01_baseline.png")
