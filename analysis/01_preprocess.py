#!/usr/bin/env python3
"""
01_preprocess.py — clean the raw trial + NASA POWER data and build analysis tables.

Inputs  (data/raw):
  tesfaye_millet_raw.csv                       plot-level MET data (agridat `tesfaye.millet`)
  nasa_power/{site}_monthly_2018_2021.csv      monthly ag-level NASA POWER, trial years
  nasa_power/{site}_monthly_1991_2024_t2m_prect.csv  monthly NASA POWER long term (T2M, PRECTOTCORR)

Outputs (data/processed):
  trial_clean.csv            cleaned plot-level data with parsed env/site/year
  trial_structure.csv        per-environment design summary (genotypes, reps, plots)
  gen_env_means.csv          genotype x environment adjusted (= arithmetic, design is
                             complete factorial) mean yield, kg/ha  (long + wide written)
  env_climate_features.csv   per-environment agroclimatic indices + 1991-2020 anomalies
  site_longterm_monthly.csv  tidy long-term monthly series per site (for trends/figures)

All derived numbers are computed from the raw CSVs only.
"""
import calendar
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

SITE_COORDS = {"Bako": (9.10, 37.15), "Assosa": (10.033, 34.567)}  # from Tesfaye et al. 2023, M&M

# ---------------------------------------------------------------- trial data
df = pd.read_csv(RAW / "tesfaye_millet_raw.csv")
df.columns = [c.strip() for c in df.columns]
for c in ["site", "rep", "gen"]:
    df[c] = df[c].astype(str).str.strip()
df["year"] = df["year"].astype(int)
df["yield"] = pd.to_numeric(df["yield"], errors="raise")
df["env"] = df["site"] + " " + df["year"].astype(str)
df = df.sort_values(["env", "gen", "rep"]).reset_index(drop=True)
assert df["yield"].notna().all() and (df["yield"] > 0).all()
df.to_csv(OUT / "trial_clean.csv", index=False)

structure = (df.groupby(["site", "year", "env"])
               .agg(n_plots=("yield", "size"),
                    n_gen=("gen", "nunique"),
                    n_rep=("rep", "nunique"),
                    mean_yield=("yield", "mean"),
                    sd_yield=("yield", "std"))
               .reset_index())
structure["plots_per_gen"] = structure["n_plots"] / structure["n_gen"]
structure.to_csv(OUT / "trial_structure.csv", index=False)
print("== trial structure ==")
print(structure.to_string(index=False, float_format=lambda x: f"{x:.1f}"))

gens_per_env = df.groupby("env")["gen"].apply(lambda s: set(s)).to_dict()
envs = sorted(gens_per_env)
print("\nenvironments:", envs)
n_envs_per_gen = df.groupby("gen")["env"].nunique().sort_values(ascending=False)
print("\ngenotypes per #environments tested:")
print(n_envs_per_gen.value_counts().sort_index().to_string())

# genotype x environment means (design is a complete factorial within each env;
# verified by plots_per_gen == n_rep for every env)
check = df.groupby(["env", "gen"])["rep"].nunique().reset_index(name="r")
rep_n = df.groupby("env")["rep"].nunique().to_dict()
complete = all(check[check.env == e]["r"].eq(rep_n[e]).all() for e in envs)
print("\ncomplete-factorial design within every environment:", complete)

gm_long = (df.groupby(["env", "site", "year", "gen"], as_index=False)
             .agg(yield_mean=("yield", "mean"), n_reps=("yield", "size")))
gm_long.to_csv(OUT / "gen_env_means.csv", index=False)
gm_wide = gm_long.pivot_table(index="gen", columns="env", values="yield_mean")
gm_wide.to_csv(OUT / "gen_env_means_wide.csv")
print("\ngenotype x environment LS-mean matrix:", gm_wide.shape,
      f"({gm_wide.notna().mean().mean()*100:.0f}% filled)")

# ---------------------------------------------------------------- climate data
MONTHS = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]

def load_monthly(path):
    d = pd.read_csv(path)
    long = d.melt(id_vars=["PARAMETER", "YEAR"], value_vars=MONTHS,
                  var_name="month", value_name="value")
    long["month"] = long["month"].map({m: i + 1 for i, m in enumerate(MONTHS)})
    wide = long.pivot_table(index=["YEAR", "month"], columns="PARAMETER",
                            values="value").reset_index()
    wide = wide.rename(columns={"YEAR": "year"})
    return wide

def add_precip_totals(w):
    days = w.apply(lambda r: calendar.monthrange(int(r["year"]), int(r["month"]))[1], axis=1)
    w = w.copy()
    w["P_mm"] = w["PRECTOTCORR"] * days          # mm/day -> mm/month
    return w

site_monthly_trial, site_monthly_long = {}, {}
for site, fname in [("Bako", "bako_monthly_2018_2021.csv"),
                    ("Assosa", "assosa_monthly_2018_2021.csv")]:
    w = add_precip_totals(load_monthly(RAW / "nasa_power" / fname))
    w["site"] = site
    site_monthly_trial[site] = w
for site, fname in [("Bako", "bako_monthly_1991_2024_t2m_prect.csv"),
                    ("Assosa", "assosa_monthly_1991_2024_t2m_prect.csv")]:
    w = add_precip_totals(load_monthly(RAW / "nasa_power" / fname))
    w["site"] = site
    site_monthly_long[site] = w

lt = pd.concat(site_monthly_long, ignore_index=True)
lt.to_csv(OUT / "site_longterm_monthly.csv", index=False)

# growing-season windows follow the documented rainfall regimes (Tesfaye et al.
# 2023): short rains Mar-May, main rainy season May-Oct (peak Jul-Aug);
# core crop window Jun-Sep.
def seasonal(w, months):
    sub = w[w["month"].isin(months)]
    g = sub.groupby("year")
    out = pd.DataFrame({
        f"P_{months[0]}_{months[-1]}": g["P_mm"].sum(),
        f"T2M_{months[0]}_{months[-1]}": g["T2M"].mean() if "T2M" in sub else np.nan,
    })
    for p in ["T2M_MAX", "T2M_MIN", "RH2M"]:
        if p in sub:
            out[f"{p}_{months[0]}_{months[-1]}"] = g[p].mean()
    return out

WINDOWS = {"MarMay": [3, 4, 5], "MayOct": [5, 6, 7, 8, 9, 10],
           "JunSep": [6, 7, 8, 9], "JulAug": [7, 8]}

feat_rows = []
for site, w34 in site_monthly_long.items():
    wt = site_monthly_trial[site]
    # 1991-2020 normals (from long series)
    norm = w34[(w34.year >= 1991) & (w34.year <= 2020)]
    norm_seas = {k: seasonal(norm, m) for k, m in WINDOWS.items()}
    for year in [2018, 2019, 2020, 2021]:
        if site == "Bako" and year == 2019:      # no Bako 2019 trial
            continue
        row = {"env": f"{site} {year}", "site": site, "year": year}
        cur_t = wt[wt.year == year]
        cur_l = w34[w34.year == year]
        for k, m in WINDOWS.items():
            s_t = seasonal(cur_t, m)     # trial file has all 5 parameters
            row[f"P_{k}"] = float(s_t[f"P_{m[0]}_{m[-1]}"].iloc[0])
            for p in ["T2M", "T2M_MAX", "T2M_MIN", "RH2M"]:
                row[f"{p}_{k}"] = float(s_t[f"{p}_{m[0]}_{m[-1]}"].iloc[0])
            n = norm_seas[k]
            row[f"P_{k}_normal"] = float(n[f"P_{m[0]}_{m[-1]}"].mean())
            row[f"T2M_{k}_normal"] = float(n[f"T2M_{m[0]}_{m[-1]}"].mean())
            row[f"P_{k}_anom"] = row[f"P_{k}"] - row[f"P_{k}_normal"]
            row[f"P_{k}_anom_pct"] = 100 * row[f"P_{k}_anom"] / row[f"P_{k}_normal"]
            row[f"T2M_{k}_anom"] = row[f"T2M_{k}"] - row[f"T2M_{k}_normal"]
        row["P_annual"] = float(cur_l["P_mm"].sum())
        row["T2M_annual"] = float(cur_l["T2M"].mean())
        feat_rows.append(row)

feat = pd.DataFrame(feat_rows).sort_values(["site", "year"]).reset_index(drop=True)
feat.to_csv(OUT / "env_climate_features.csv", index=False)
cols = ["env", "P_MarMay", "P_MayOct", "P_JunSep", "P_JunSep_anom_pct",
        "T2M_JunSep", "T2M_MAX_JunSep", "T2M_JunSep_anom", "RH2M_JunSep"]
print("\n== environment climate features (NASA POWER) ==")
print(feat[cols].to_string(index=False, float_format=lambda x: f"{x:.1f}"))

dry = feat.loc[feat["P_MayOct"].idxmin(), "env"]
wet = feat.loc[feat["P_MayOct"].idxmax(), "env"]
print(f"\ndriest trial environment (May-Oct rain): {dry}; wettest: {wet}")
print("\nwrote:", [p.name for p in sorted(OUT.glob('*.csv'))])
