#!/usr/bin/env python3
"""
02_gxe_stability.py — variance partitioning, stability analyses, AMMI, and the topicandduration
continuously from a normal to prevent climate-GxE linkage analysis.

All inputs come from data/processed (produced by 01_preprocess.py).
Outputs are written to results/tables/.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats as st
import statsmodels.api as sm
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[1]
PRO = ROOT / "data" / "processed"
TAB = ROOT / "results" / "tables"
TAB.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(7)

df = pd.read_csv(PRO / "trial_clean.csv")
gm = pd.read_csv(PRO / "gen_env_means.csv")
feat = pd.read_csv(PRO / "env_climate_features.csv")
lt = pd.read_csv(PRO / "site_longterm_monthly.csv")
df["gen"] = df["gen"].astype(str); gm["gen"] = gm["gen"].astype(str)

summary = {}

# ================================================================ 1. ANOVA + variance components
# pooled ANOVA via nested-model comparisons (robust for the rank-deficient,
# unbalanced GxE design): SS(term) = SSE(reduced) - SSE(full), df = rank difference
df["yld"] = df["yield"]   # 'yield' is a Python keyword -> invalid in patsy formulas

def fit_sse(formula):
    m = smf.ols(formula, data=df).fit()
    r = np.linalg.matrix_rank(m.model.exog)
    return float(m.ssr), int(r), float(m.df_resid)

TERMS = [("Environment (E)", "yld ~ 1", "yld ~ C(env)"),
         ("Genotype (G | E)", "yld ~ C(env)", "yld ~ C(env) + C(gen)"),
         ("Rep within E", "yld ~ C(env) + C(gen)", "yld ~ C(env) + C(gen) + C(env):C(rep)"),
         ("GxE (G:E | E+G)", "yld ~ C(env) + C(gen) + C(env):C(rep)",
          "yld ~ C(env) + C(gen) + C(env):C(rep) + C(env):C(gen)")]
m_full = smf.ols("yld ~ C(env) + C(gen) + C(env):C(rep) + C(env):C(gen)", data=df).fit()
sse_full, df_full = float(m_full.ssr), int(m_full.df_resid)
mse_full = sse_full / df_full
rows = []
for name, f0, f1 in TERMS:
    s0, r0, _ = fit_sse(f0)
    s1, r1, _ = fit_sse(f1)
    ss, dft = s0 - s1, r1 - r0
    F = (ss / dft) / mse_full
    rows.append({"term": name, "df": dft, "SS": ss, "MS": ss / dft,
                 "F": F, "p": 1 - st.f.cdf(F, dft, df_full)})
rows.append({"term": "Residual", "df": df_full, "SS": sse_full, "MS": mse_full,
             "F": np.nan, "p": np.nan})
aov = pd.DataFrame(rows)
aov.to_csv(TAB / "anova_model_comparison.csv", index=False)
ss_tot = aov["SS"].sum()
summary["pct_SS_gxe"] = round(100 * aov.loc[3, "SS"] / ss_tot, 2)
summary["pct_SS_env"] = round(100 * aov.loc[0, "SS"] / ss_tot, 2)
summary["pct_SS_gen"] = round(100 * aov.loc[1, "SS"] / ss_tot, 2)
summary["anova_F"] = {r.term: round(r.F, 2) for r in aov.itertuples() if r.term != "Residual"}
summary["anova_p"] = {r.term: float(r.p) for r in aov.itertuples() if r.term != "Residual"}

# variance components: yield ~ 1 + (gen) + (env) + (gen:env) + (env:rep), REML
df["gxe"] = df["gen"] + "|" + df["env"]
df["rep_e"] = df["rep"].astype(str) + "|" + df["env"]
try:
    md = smf.mixedlm("yld ~ 1", data=df, groups=df.index.to_series(),
                     vc_formula={"gen": "0 + C(gen)", "env": "0 + C(env)",
                                 "gxe": "0 + C(gxe)", "rep_in_env": "0 + C(rep_e)"})
    m_vc = md.fit(reml=True)
    vc = m_vc.vcomp
    s2 = {k: float(v) for k, v in zip(["gen", "env", "gxe", "rep_in_env"], vc)}
    s2e = float(m_vc.scale)
    vc_method = "REML MixedLM"
except Exception as e:   # pragma: no cover - fallback to method of moments on balanced subset
    vc_method = f"method-of-moments fallback ({type(e).__name__}), complete 11 x 7 subset"
    wide_tmp = gm.pivot_table(index="gen", columns="env", values="yield_mean")
    full_g = wide_tmp.dropna().index
    dsb = df[df.gen.isin(full_g) & df.rep.astype(str).isin(["R1", "R2", "R3"])].copy()
    a2 = sm.stats.anova_lm(smf.ols("yld ~ C(env) + C(gen) + C(env):C(gen)", data=dsb).fit(), typ=2)
    # actual mean number of replicate plots per genotype-environment cell
    n_rep_cell: float = dsb.groupby(["gen", "env"]).size().mean()
    E_ = dsb.env.nunique(); G_ = dsb.gen.nunique()
    MSg = a2.loc["C(gen)", "sum_sq"] / a2.loc["C(gen)", "df"]
    MSge = a2.loc["C(env):C(gen)", "sum_sq"] / a2.loc["C(env):C(gen)", "df"]
    MSe = a2.loc["Residual", "sum_sq"] / a2.loc["Residual", "df"]
    s2 = {"gen": max((MSg - MSge) / (E_ * n_rep_cell), 0.0),
          "env": np.nan, "gxe": max((MSge - MSe) / n_rep_cell, 0.0),
          "rep_in_env": np.nan}
    s2e = float(MSe)
summary["vc_method"] = vc_method
summary["variance_components_kgha2"] = {**s2, "residual": round(s2e, 1)}
e_bar = df.groupby("gen")["env"].nunique().mean()          # avg envs per genotype
r_bar = df.groupby(["gen", "env"])["rep"].count().mean()   # avg reps per cell
H2 = s2["gen"] / (s2["gen"] + s2["gxe"] / e_bar + s2e / (e_bar * r_bar))
summary["H2_entry_mean_mixedmodel"] = round(H2, 3)

# per-environment ANOVA + heritability (entry-mean basis)
h2_rows = []
for env, d in df.groupby("env"):
    a = sm.stats.anova_lm(smf.ols("yld ~ C(rep) + C(gen)", data=d).fit(), typ=2)
    msg, mse = a.loc["C(gen)", "sum_sq"] / a.loc["C(gen)", "df"], a.loc["Residual", "sum_sq"] / a.loc["Residual", "df"]
    F = msg / mse
    df_envir = {"env": env, "site": d.site.iloc[0], "year": int(d.year.iloc[0]),
                "n_gen": d.gen.nunique(), "n_rep": d.rep.nunique(),
                "MS_gen": msg, "MS_error": mse, "F_gen": F,
                "p_gen": 1 - st.f.cdf(F, a.loc["C(gen)", "df"], a.loc["Residual", "df"]),
                "H2_entry_mean": (msg - mse) / msg,
                "mean_yield": d["yield"].mean()}
    h2_rows.append(df_envir)
h2 = pd.DataFrame(h2_rows).sort_values(["site", "year"])
h2.to_csv(TAB / "per_env_anova_heritability.csv", index=False)

# ================================================================ 2. Finlay-Wilkinson stability
wide = gm.pivot_table(index="gen", columns="env", values="yield_mean")
env_index = wide.mean(axis=0)                      # environmental index (kg/ha)
fw_rows = []
for gname, row in wide.iterrows():
    y = row.dropna()
    if len(y) < 3:
        continue
    x = env_index[y.index].to_numpy(); yy = y.to_numpy()
    sl, ic, r, p, se = st.linregress(x, yy)
    yhat = ic + sl * x
    s2d = np.sum((yy - yhat) ** 2) / (len(yy) - 2)
    fw_rows.append({"gen": gname, "n_env": len(y), "mean": yy.mean(), "min_env_tested": ",".join(y.index),
                    "fw_slope_b": sl, "fw_intercept": ic, "fw_R2": r ** 2,
                    "s2d_deviation": s2d, "b_pvalue": p})
fw = pd.DataFrame(fw_rows).sort_values("mean", ascending=False)
fw.to_csv(TAB / "finlay_wilkinson_stability.csv", index=False)
summary["n_genotypes_FW"] = len(fw)
summary["grand_mean_kgha"] = round(float(wide.stack().mean()), 1)

# ================================================================ 3. Shukla + Wricke on complete 11x7 subset
n_env_per_gen = wide.notna().sum(axis=1)
complete_gens = n_env_per_gen[n_env_per_gen == wide.shape[1]].index
sub = wide.loc[complete_gens]
G, E = sub.shape
grand = sub.to_numpy().mean()
g_eff = sub.mean(axis=1) - grand
e_eff = sub.mean(axis=0) - grand
resid = sub.subtract(sub.mean(axis=1), axis=0).subtract(sub.mean(axis=0), axis=1) + grand
shukla = (resid ** 2).sum(axis=1) * G / ((G - 1) * (E - 2) / (G - 1))
# standard Shukla variance: sigma_i^2 = G/((G-2)(E-1)) * sum_j resid_ij^2 - SS_GE/((G-1)(G-2)(E-1))
SSge = (resid ** 2).sum().sum()
shukla_raw = (resid ** 2).sum(axis=1)
shukla_var = G / ((G - 2) * (E - 1)) * shukla_raw - SSge / ((G - 1) * (E - 1) * (G - 2))
ecov = shukla_raw / SSge * 100
stab_complete = pd.DataFrame({"gen": complete_gens, "mean": sub.mean(axis=1).values,
                              "shukla_var": np.maximum(shukla_var, 0), "wricke_pct": ecov.values})
stab_complete = stab_complete.sort_values("mean", ascending=False)
stab_complete.to_csv(TAB / "stability_complete_subset_11gen.csv", index=False)

# ================================================================ 4. AMMI on complete subset
Z = resid.to_numpy()                                # double-centered GxE
U, S, Vt = np.linalg.svd(Z, full_matrices=False)
lam = S ** 2
expl = 100 * lam / lam.sum()
gen_scores = U[:, :3] * S[:3]                       # principal coordinate for genotypes
env_scores = Vt.T[:, :3] * S[:3]
ammi_gen = pd.DataFrame({"gen": sub.index, "IPCA1": gen_scores[:, 0],
                         "IPCA2": gen_scores[:, 1], "IPCA3": gen_scores[:, 2],
                         "mean": sub.mean(axis=1).values})
ammi_env = pd.DataFrame({"env": sub.columns, "IPCA1": env_scores[:, 0],
                         "IPCA2": env_scores[:, 1], "IPCA3": env_scores[:, 2],
                         "mean": sub.mean(axis=0).values})
ammi_gen.to_csv(TAB / "ammi_genotype_scores.csv", index=False)
ammi_env.to_csv(TAB / "ammi_environment_scores.csv", index=False)
summary["AMMI_IPCA_var_explained_pct"] = [round(float(x), 1) for x in expl[:3]]

# ================================================================ 5. environment similarity vs climate distance + Mantel
spear = pd.DataFrame(np.full((7, 7), np.nan), columns=sub.columns, index=sub.columns)
for e1 in sub.columns:
    for e2 in sub.columns:
        if e1 == e2:
            spear.loc[e1, e2] = 1.0
            continue
        pair = wide[[e1, e2]].dropna()
        spear.loc[e1, e2] = st.spearmanr(pair[e1], pair[e2]).statistic
spear.to_csv(TAB / "env_spearman_rankcorr.csv")

CLIM = ["P_MarMay", "P_JunSep", "P_MayOct", "T2M_JunSep", "T2M_MAX_JunSep", "RH2M_JunSep"]
X = feat.set_index("env")[CLIM]
Xz = (X - X.mean()) / X.std(ddof=0)
cdist = pd.DataFrame([[np.linalg.norm(Xz.loc[a] - Xz.loc[b]) for b in Xz.index] for a in Xz.index],
                     index=Xz.index, columns=Xz.index)
cdist.to_csv(TAB / "env_climate_distance.csv")

iu = np.triu_indices(7, 1)
ph_dist = (1 - spear).to_numpy()[iu]                # phenotypic dissimilarity
cl_dist = cdist.to_numpy()[iu]
r_obs = st.pearsonr(ph_dist, cl_dist).statistic
perms = []
M = (1 - spear).to_numpy()
for _ in range(4999):
    p = rng.permutation(7)
    perms.append(st.pearsonr(M[np.ix_(p, p)][iu], cl_dist).statistic)
p_mantel = (np.sum(np.array(perms) >= r_obs) + 1) / (len(perms) + 1)
summary["mantel_r_phendist_climdist"] = round(float(r_obs), 3)
summary["mantel_p"] = round(float(p_mantel), 4)

# ================================================================ 6. climate -> environment means & GxE axis
envs_all = env_index.index
link = feat.set_index("env").loc[envs_all]
rows = []
for c in CLIM + ["P_JunSep_anom_pct", "P_MarMay", "T2M_JunSep_anom"]:
    rows.append({"climate_feature": c,
                 "pearson_r_envmean": st.pearsonr(link[c], env_index).statistic,
                 "pearson_p_envmean": st.pearsonr(link[c], env_index).pvalue,
                 "pearson_r_IPCA1": st.pearsonr(link[c], ammi_env.set_index("env").loc[envs_all, "IPCA1"]).statistic,
                 "spearman_r_IPCA1": st.spearmanr(link[c], ammi_env.set_index("env").loc[envs_all, "IPCA1"]).statistic})
link_tab = pd.DataFrame(rows)
link_tab.to_csv(TAB / "climate_links_envmean_IPCA1.csv", index=False)

# ================================================================ 7. long-term climate trends (1991-2024)
DAYS = {m: [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31][m - 1] for m in range(1, 13)}
trend_rows = []
for site, d in lt.groupby("site"):
    d = d.copy()
    d["year"] = d["year"].astype(int)
    ann = d.groupby("year").agg(T2M=("T2M", "mean"))
    p_seas = d[d.month.isin([6, 7, 8, 9])].groupby("year").apply(
        lambda x: np.sum(x["PRECTOTCORR"] * [DAYS[m] for m in x.month]), include_groups=False)
    t_row = {"site": site}
    for nm, ser in [("T2M_annual", ann.T2M), ("P_JunSep", p_seas)]:
        sl, ic, r, p, se = st.linregress(ser.index.to_numpy(), ser.to_numpy())
        ts = st.theilslopes(ser.to_numpy(), ser.index.to_numpy())
        tau, p_kt = st.kendalltau(ser.index.to_numpy(), ser.to_numpy())
        trend_rows.append({"site": site, "series": nm, "period": "1991-2024",
                           "slope_per_decade": sl * 10, "p_ols": p,
                           "theil_sen_per_decade": ts.slope * 10, "p_kendall": p_kt})
trend = pd.DataFrame(trend_rows)
trend.to_csv(TAB / "longterm_climate_trends.csv", index=False)
summary["warming_1991_2024_C_per_decade"] = {r.site: round(r.slope_per_decade, 3)
                                             for r in trend.itertuples() if r.series == "T2M_annual"}

# ================================================================ 8. climate-resilient donor panel
driest = feat.loc[feat["P_JunSep"].idxmin(), "env"]
summary["driest_environment"] = str(driest)
panel = fw[fw.n_env >= 5].copy()
panel["yield_driest"] = panel["gen"].map(wide[driest] if driest in wide.columns else pd.Series(dtype=float))
panel["rank_driest"] = panel["yield_driest"].rank(ascending=False)
stab_map = stab_complete.set_index("gen")[["shukla_var", "wricke_pct"]]
panel = panel.join(stab_map, on="gen")
gm_all = float(wide.stack().mean())
panel["resilient_donor"] = ((panel["mean"] >= gm_all) & (panel["fw_slope_b"] <= 1.05) &
                            (panel["rank_driest"] <= 5)).map({True: "yes", False: "no"})
panel = panel.sort_values("mean", ascending=False)
panel.to_csv(TAB / "resilient_donor_panel.csv", index=False)
summary["resilient_donors"] = panel.loc[panel.resilient_donor == "yes", "gen"].tolist()

# top performers per environment (winners)
winners = wide.idxmax()
summary["env_winners"] = {e: str(winners[e]) for e in wide.columns}

with open(TAB / "summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print("== nested-model-comparison ANOVA =="); print(aov.round(3).to_string(index=False))
print("\n== variance components (kg/ha)^2 ==");
for k, v in s2.items(): print(f"  {k:>10}: {v:10.1f}")
print(f"  {'residual':>10}: {s2e:10.1f}")
print(f"  H2 (entry mean): {H2:.3f}")
print("\n== per-environment heritability =="); print(h2.round(3).to_string(index=False))
print("\n== AMMI %var explained:", [round(float(x), 1) for x in expl[:4]])
print("== environment Spearman rank correlations =="); print(spear.round(2).to_string())
print("== climate distance (z, May-Oct / Jun-Sep window features) =="); print(cdist.round(2).to_string())
print(f"\n== Mantel: r = {r_obs:.3f}, p = {p_mantel:.4f}")
print("\n== climate links =="); print(link_tab.round(3).to_string(index=False))
print("\n== long-term trends =="); print(trend.round(4).to_string(index=False))
print("\n== donor panel (genotypes tested in >=5 envs) ==")
print(panel[["gen", "n_env", "mean", "fw_slope_b", "yield_driest", "rank_driest",
             "wricke_pct", "resilient_donor"]].round(2).to_string(index=False))
print("\n== environment winners ==", summary["env_winners"])
