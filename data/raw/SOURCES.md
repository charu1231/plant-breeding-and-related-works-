# Data sources and provenance

All data in this project were pulled from public, open sources on **2026-08-31** (UTC).
No observations were simulated, imputed at source, or edited by hand.

## 1. Multi-environment trial data (phenotypes)

`tesfaye_millet_raw.csv` (415 plot-level records)

- Original study: Tesfaye K, Alemu T, Argaw T, de Villiers S, Assefa E (2023).
  *Evaluation of finger millet (Eleusine coracana (L.) Gaertn.) in multi-environment
  trials using enhanced statistical models.* PLoS ONE 18(2): e0277499.
  DOI: https://doi.org/10.1371/journal.pone.0277499
  License: Creative Commons Attribution (CC BY 4.0). Data are within the paper and
  its Supporting Information files.
- Distribution copy used here: R package **agridat** (v1.27, MIT license),
  dataset `tesfaye.millet`, source repository https://github.com/kwstat/agridat
  (file `data/tesfaye.millet.txt`, converted verbatim to CSV).
- Trials: Bako Agricultural Research Center (9°6' N, 37°9' E) and Assosa
  Agricultural Research Center (10°02' N, 34°34' E), Ethiopia, 2018–2021,
  7 site-year environments, RCBD, grain yield (kg/ha).
- Site coordinates and site climate descriptions above are quoted from the
  Materials and Methods section of Tesfaye et al. (2023).

## 2. Agro-climate data

NASA POWER (Prediction Of Worldwide Energy Resources), Agroclimatology (AG)
community, Data Access API v2. Files were saved from API responses exactly as
returned (`format=CSV`).

| File | Endpoint (template) | Window | Parameters |
|---|---|---|---|
| `nasa_power/bako_monthly_2018_2021.csv` | `https://power.larc.nasa.gov/api/temporal/monthly/point?longitude=37.15&latitude=9.10&community=AG&temporal=monthly&start=2018&end=2021&format=CSV&header=false&time-standard=LST` | monthly 2018–2021 | T2M, T2M_MAX, T2M_MIN, PRECTOTCORR, RH2M |
| `nasa_power/assosa_monthly_2018_2021.csv` | same, `longitude=34.567&latitude=10.033` | monthly 2018–2021 | T2M, T2M_MAX, T2M_MIN, PRECTOTCORR, RH2M |
| `nasa_power/bako_monthly_1991_2024_t2m_prect.csv` | same, `start=1991&end=2024` | monthly 1991–2024 | T2M, PRECTOTCORR |
| `nasa_power/assosa_monthly_1991_2024_t2m_prect.csv` | same, `start=1991&end=2024` | monthly 1991–2024 | T2M, PRECTOTCORR |

Units (NASA POWER AG conventions): temperatures in °C (monthly means of daily
mean / daily max / daily min), RH2M in %, PRECTOTCORR in mm/day (monthly mean
of corrected daily precipitation; multiplied by days per month to obtain
monthly totals). Grid cell ~0.5° × 0.625° (MERRA-2 / GEOS-IT assimilation
products). Citation: NASA Langley Research Center POWER Project,
https://power.larc.nasa.gov (accessed 2026-08-31).
