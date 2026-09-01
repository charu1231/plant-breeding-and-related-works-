# Climate-driven re-analysis of a rarely used finger millet MET dataset

This repository contains a complete, reproducible secondary analysis of an **open,
rarely reused plant-breeding dataset** — 7-environment finger millet
(*Eleusine coracana*) multi-environment trial (MET) data from Ethiopia — joined
with **NASA POWER satellite/reanalysis agro-climate data**. **Every number in
this project is derived from open data pulled programmatically as CSV; nothing
was fabricated or typed in by hand.**

Full provenance: [`data/raw/SOURCES.md`](data/raw/SOURCES.md).
Manuscript draft (PDF/DOCX/Markdown): **held outside the public repository**
(unpublished work; available from the authors on request).

## The dataset

| Layer | Source | Content |
|---|---|---|
| Phenotypes | Tesfaye et al. (2023), *PLoS ONE* 18(2):e0277499 (CC-BY), distributed via the R package `agridat` v1.27 (`tesfaye.millet`) | 415 plot records, 47 genotypes, 2 sites (Bako, Assosa), 2018–2021, grain yield |
| Climate | NASA POWER Data Access API v2, AG community | monthly T2M / T2M_MAX / T2M_MIN / RH2M / PRECTOTCORR for both trial sites (trial years + 1991–2024) |

## Reproduce

```bash
python3 -m pip install pandas numpy matplotlib scipy statsmodels
python3 analysis/01_preprocess.py     # raw CSV -> processed tables
python3 analysis/02_gxe_stability.py  # ANOVA, variance components, stability, AMMI, Mantel
python3 analysis/03_figures.py        # figures 1-6
```

Re-pulling the raw trial CSV requires the `agridat` GitHub mirror
(`data/tesfaye.millet.txt`); the NASA POWER files can be re-downloaded with the
URLs in `data/raw/SOURCES.md`.

## Headline findings

* 2018 was a marked drought year at both stations (Jun–Sep rainfall anomalies
  −16 % at Assosa, −25 % at Bako vs 1991–2020), spanning a useful stress →
  favorable gradient for G×E analysis.
* Environment, genotype and G×E accounted for 60.6 %, 12.2 % and 11.4 % of the
  total yield variation; entry-mean heritability ≈ 0.62.
* The driest Assosa environment had the **highest** heritability in the whole
  series (H² = 0.93) whereas the drought year at Bako erased genetic signal
  (H² = 0.16, n.s.) — identical water stress, opposite selection value.
* Genotype rankings were consistent across Assosa seasons (up to ρ = 0.88) but
  broke down at Bako 2018 (as low as ρ = −0.45 vs Assosa 2020).
* Seasonal rainfall + humidity explain environment productivity levels well,
  but overall G×E *structure* was NOT significantly associated with the climate
  distance between environments (Mantel r = 0.06, p = 0.39) — an honest null
  result with methodological implications for characterizing METs.
* Stability screening (Finlay–Wilkinson + rank persistence in the driest
  environment) nominates **203347, 203364 and 203263** as climate-resilient
  donor parents; the original study's top pick "Bako-09" is shown to be highly
  unstable (b = 1.66) and suited only to favorable environments.
* NASA POWER data also show a significant warming trend at Bako
  (+0.15 °C/decade, Theil–Sen, p = 0.03 over 1991–2024).

## Layout

```
data/raw/               pristine CSVs + SOURCES.md (provenance)
data/processed/         cleaned / derived tables
analysis/               01-03 python scripts (this is the full pipeline)
results/tables/         every number quoted in the paper
results/figures/        figures 1-6
paper/PAPER.md          manuscript draft referencing the artifacts above
```

*This repository was assembled with the help of an AI coding agent; all
statistical outputs are produced by the scripts above from the cited open data
and can be audited end-to-end.*
