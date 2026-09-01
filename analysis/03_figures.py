#!/usr/bin/env python3
"""03_figures.py — render all paper figures from processed data + result tables."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRO = ROOT / "data" / "processed"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

CB = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#000000"]
plt.rcParams.update({"font.size": 9.5, "axes.titlesize": 10.5, "axes.labelsize": 9.5,
                     "figure.dpi": 300, "savefig.dpi": 300})

df = pd.read_csv(PRO / "trial_clean.csv")
gm = pd.read_csv(PRO / "gen_env_means.csv", dtype={"gen": str})
feat = pd.read_csv(PRO / "env_climate_features.csv")
lt = pd.read_csv(PRO / "site_longterm_monthly.csv")
h2 = pd.read_csv(TAB / "per_env_anova_heritability.csv")
fw = pd.read_csv(TAB / "finlay_wilkinson_stability.csv", dtype={"gen": str})
panel = pd.read_csv(TAB / "resilient_donor_panel.csv", dtype={"gen": str})
ammi_g = pd.read_csv(TAB / "ammi_genotype_scores.csv", dtype={"gen": str})
ammi_e = pd.read_csv(TAB / "ammi_environment_scores.csv")
spear = pd.read_csv(TAB / "env_spearman_rankcorr.csv", index_col=0)
cdist = pd.read_csv(TAB / "env_climate_distance.csv", index_col=0)

wide = gm.pivot_table(index="gen", columns="env", values="yield_mean")
MONTHS = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]

# ---------------------------------------------------------------- fig 1: climate context
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
colors = {"Bako": CB[0], "Assosa": CB[1]}
for site, d in lt.groupby("site"):
    ann = d.groupby("year")["T2M"].mean()
    ax[0].plot(ann.index, ann.values, ".", ms=3, color=colors[site], alpha=0.45)
    z = np.polyfit(ann.index, ann.values, 1)
    ax[0].plot(ann.index, np.polyval(z, ann.index), color=colors[site], lw=1.6,
               label=f"{site} ({z[0]*10:+.2f} °C/decade)")
ax[0].axvspan(2018, 2021, color="grey", alpha=0.15)
ax[0].text(2019.5, ax[0].get_ylim()[0], "trial\nperiod", ha="center", va="bottom", fontsize=8)
ax[0].set_xlabel("Year"); ax[0].set_ylabel("Mean annual temperature (°C)")
ax[0].set_title("(a) Annual temperature 1991–2024 (NASA POWER)")
ax[0].legend(frameon=False, fontsize=8.5)

seas = feat.sort_values(["site", "year"])
xloc = np.arange(len(seas))
ax[1].bar(xloc, seas["P_JunSep"], color=[colors[s] for s in seas["site"]], width=0.65)
for i, (s, n) in enumerate(zip(seas["site"], seas["P_JunSep_normal"])):
    ax[1].plot([i - 0.33, i + 0.33], [n, n], color="black", lw=1.2)
ax[1].set_xticks(xloc); ax[1].set_xticklabels(seas["env"], rotation=35, ha="right", fontsize=8)
ax[1].set_ylabel("Jun–Sep rainfall (mm)")
ax[1].set_title("(b) Main-season rainfall in trial years vs 1991–2020 normal")
ax[1].legend(handles=[Line2D([0], [0], color="black", lw=1.2, label="1991–2020 normal")],
             frameon=False, fontsize=8.5)

for k, site in enumerate(["Bako", "Assosa"]):
    d = lt[lt.site == site].copy()
    norm = d[(d.year >= 1991) & (d.year <= 2020)].groupby("month")["P_mm"].mean()
    axx = ax[2]
    axx.plot(range(1, 13), norm.values, color=colors[site], lw=2,
             ls="-" if site == "Bako" else "--", label=site)
for e, c in zip(["2018", "2021"], [CB[2], CB[3]]):
    for site in ["Bako", "Assosa"]:
        d = lt[(lt.site == site) & (lt.year == int(e))]
        ax[2].plot(d.month, d.P_mm, ".", ms=4, alpha=0.55,
                   color=colors[site], marker="x" if e == "2018" else "o")
ax[2].axvspan(6, 9, color="grey", alpha=0.12)
ax[2].set_xticks(range(1, 13)); ax[2].set_xticklabels(MONTHS, fontsize=7.5, rotation=45)
ax[2].set_xlabel("Month"); ax[2].set_ylabel("Rainfall (mm/month)")
ax[2].set_title("(c) Monthly rainfall climatology (lines); 2018 = ×, 2021 = ○")
ax[2].legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig(FIG / "fig1_climate_context.png"); plt.close(fig)

# ---------------------------------------------------------------- fig 2: env performance + H2 + rainfall
envs = list(h2["env"])
fig, ax1 = plt.subplots(figsize=(8.2, 4.6))
cols = [colors[s] for s in h2["site"]]
bars = ax1.bar(range(len(envs)), h2["mean_yield"], color=cols, width=0.66)
for i, (v, hh, p, pv) in enumerate(zip(h2["mean_yield"], h2["H2_entry_mean"], h2["P_JunSep"] if "P_JunSep" in h2 else [0]*len(h2), h2["p_gen"])):
    star = "***" if pv < 0.001 else "**" if pv < 0.01 else "*" if pv < 0.05 else "ns"
    ax1.text(i, v + 40, f"H²={hh:.2f}\n{star}", ha="center", fontsize=8)
ax1.set_xticks(range(len(envs))); ax1.set_xticklabels(envs, rotation=35, ha="right")
ax1.set_ylabel("Mean grain yield (kg/ha)"); ax1.set_ylim(0, max(h2["mean_yield"]) * 1.32)
f2 = feat.set_index("env").loc[envs]
ax2 = ax1.twinx()
ax2.plot(range(len(envs)), f2["P_JunSep"], "o-", color="black", lw=1.4, ms=5)
ax2.set_ylabel("Jun–Sep rainfall (mm)", color="black")
ax2.tick_params(axis="y", labelcolor="black")
ax1.set_title("Environment mean yield, genotype F-test significance, heritability, and seasonal rainfall")
fig.tight_layout(); fig.savefig(FIG / "fig2_environment_performance.png"); plt.close(fig)

# ---------------------------------------------------------------- fig 3: AMMI1-style biplot
fig, ax = plt.subplots(figsize=(7.6, 5.6))
ax.axhline(0, color="grey", lw=0.7); ax.axvline(wide.stack().mean(), color="grey", lw=0.7, ls=":")
gmm = wide.stack().mean()
ax.scatter(ammi_g["mean"], ammi_g["IPCA1"], s=26, color="#4C72B0", alpha=0.85, zorder=3, label="genotypes")
for r in ammi_g.itertuples():
    if abs(r.IPCA1) > ammi_g.IPCA1.abs().quantile(0.72) or r.mean > gmm:
        ax.annotate(str(r.gen), (r.mean, r.IPCA1), fontsize=7.5, xytext=(3, 3),
                    textcoords="offset points")
ax.scatter(ammi_e["mean"], ammi_e["IPCA1"], s=95, color="#C44E52", marker="s", zorder=4, label="environments")
for r in ammi_e.itertuples():
    ax.annotate(r.env, (r.mean, r.IPCA1), fontsize=8.5, weight="bold", color="#C44E52",
                xytext=(5, 5), textcoords="offset points")
ax.set_xlabel("Main effect: mean yield (kg/ha)")
ax.set_ylabel("IPCA1 score (56.3% of G×E variation)")
ax.set_title("AMMI biplot (complete 11-genotype × 7-environment subset)")
ax.legend(frameon=False, loc="lower right")
fig.tight_layout(); fig.savefig(FIG / "fig3_ammi_biplot.png"); plt.close(fig)

# ---------------------------------------------------------------- fig 4: Finlay-Wilkinson stability map
p5 = panel.copy()
fig, ax = plt.subplots(figsize=(7.6, 5.6))
gmm = p5["mean"].mean() * 0 + wide.stack().mean()
ax.axhline(1.0, color="grey", ls="--", lw=1); ax.axvline(gmm, color="grey", ls="--", lw=1)
don = p5["resilient_donor"] == "yes"
ax.scatter(p5.loc[~don, "mean"], p5.loc[~don, "fw_slope_b"], s=45, color="#888888", zorder=3)
ax.scatter(p5.loc[don, "mean"], p5.loc[don, "fw_slope_b"], s=80, color=CB[2], marker="*", zorder=4)
for r in p5.itertuples():
    ax.annotate(str(r.gen), (r.mean, r.fw_slope_b), fontsize=8, xytext=(4, 4), textcoords="offset points")
ax.set_xlabel("Genotype mean yield across environments (kg/ha)")
ax.set_ylabel("Finlay–Wilkinson slope b")
ax.set_title("Yield-stability map of genotypes tested in ≥5 environments")
ax.text(0.02, 0.03, "stable & high-yielding (donors, ★)\nbelow b = 1, right of grand mean",
        transform=ax.transAxes, fontsize=8, color=CB[2])
fig.tight_layout(); fig.savefig(FIG / "fig4_stability_finlay_wilkinson.png"); plt.close(fig)

# ---------------------------------------------------------------- fig 5: env similarity vs climate distance matrics
S = (spear.to_numpy() + spear.to_numpy().T) / 2
diss = 1 - S
np.fill_diagonal(diss, 0.0)
cond = squareform(diss, checks=False)
L = linkage(cond, method="average")
order = np.asarray(leaves_list(L), dtype=int)
lab = [str(spear.columns[i]) for i in order]
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6))
for axx, M, ttl, cmap, vmin, vmax in [(ax[0], spear.loc[lab, lab].to_numpy(),
                                       "(a) G×E structure: genotype-rank\nSpearman correlation between environments", "viridis", -0.6, 1.0),
                                      (ax[1], cdist.loc[lab, lab].to_numpy(),
                                       "(b) NASA POWER standardized climate\ndistance between environments", "magma_r", 0, 7)]:
    im = axx.imshow(M, cmap=cmap, vmin=vmin, vmax=vmax)
    axx.set_xticks(range(len(lab))); axx.set_xticklabels(lab, rotation=35, ha="right", fontsize=7.5)
    axx.set_yticks(range(len(lab))); axx.set_yticklabels(lab, fontsize=7.5)
    axx.set_title(ttl)
    for i in range(len(lab)):
        for j in range(len(lab)):
            axx.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=6,
                     color="white" if (M[i, j] - vmin) / (vmax - vmin) < 0.55 else "black")
    fig.colorbar(im, ax=axx, fraction=0.045, pad=0.03)
fig.suptitle("Phenotypic vs climatic similarity of trial environments (Mantel r = 0.06, p = 0.39)")
fig.tight_layout(rect=[0, 0, 1, 0.93]); fig.savefig(FIG / "fig5_env_similarity_vs_climate.png"); plt.close(fig)

# ---------------------------------------------------------------- fig 6: rank crossover heatmap
ranks = wide.loc[p5["gen"]].rank(ascending=False).loc[:, feat.sort_values("P_JunSep")["env"]]
fig, ax = plt.subplots(figsize=(8.6, 5.2))
M = ranks.to_numpy().astype(float)
im = ax.imshow(M, cmap="RdYlGn_r", vmin=1, vmax=ranks.shape[0], aspect="auto")
ax.set_xticks(range(len(ranks.columns))); ax.set_xticklabels(ranks.columns, rotation=35, ha="right")
ax.set_yticks(range(len(ranks.index)))
ax.set_yticklabels(["★ " + g if g in set(panel.loc[panel.resilient_donor == "yes", "gen"]) else g
                    for g in ranks.index], fontsize=8)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j, i, f"{M[i, j]:.0f}", ha="center", va="center", fontsize=7,
                color="white" if M[i, j] < 4 or M[i, j] > 10 else "black")
ax.set_title("Genotype yield ranks across environments (ordered by Jun–Sep rainfall, wettest on right)")
fig.colorbar(im, ax=ax, label="rank (1 = best)", fraction=0.045, pad=0.03)
fig.tight_layout(); fig.savefig(FIG / "fig6_rank_crossover_heatmap.png"); plt.close(fig)

print("figures written:", sorted(p.name for p in FIG.glob("*.png")))
