# Multi-environment Genomic Prediction under G×E — Simulation + Real-data validation

**Working title:** *Multi-environment genomic prediction under genotype × environment
interaction: a simulation study validated in CIMMYT wheat.*

A fully reproducible, pure-Python quantitative-genetics pipeline that benchmarks
multi-environment genomic-prediction models across a controlled G×E gradient and
validates the conclusions on the public CIMMYT wheat panel (Crossa et al. 2010,
*Genetics* 186:713–724).

## Data
- `data/wheat.RData` — CIMMYT wheat: 599 lines × 1,279 DArT markers × 4 environments
  (Crossa et al. 2010, via the public `gdlc/BGLR` GitHub repository).

## Setup (sandbox: only PyPI + GitHub reachable; no R needed)
```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Pipeline
```bash
.venv/bin/python analysis/01_real_data_baseline.py      # G×E structure + GBLUP CV
.venv/bin/python analysis/02_simulation_study.py        # rG × h² benchmark (15 reps)
.venv/bin/python analysis/03_real_data_validation.py    # MT vs baselines on wheat
.venv/bin/python analysis/04_robustness_structured.py   # reaction-norm G×E (estimated)
.venv/bin/python analysis/05_robustness_unbalanced.py   # 30% missing-data check
.venv/bin/python analysis/07_robustness_factor_analytic.py  # FA rank analysis
.venv/bin/python analysis/06_summary_figure.py          # combined summary figure
```

## Headline results
| Setting | Best model | Key number |
|---|---|---|
| Simulation, new env (CV1) | MT-GBLUP | +0.03–0.12 over naive mean |
| Simulation, new lines (CV2) | MT-GBLUP | +0.005–0.10 over single-env GBLUP |
| Real wheat, new env (CV1) | MT-GBLUP | r = 0.379 vs 0.313 (naive) |
| Real wheat, new lines (CV2) | MT-GBLUP ≈ single | r = 0.443 vs 0.440 |
| Structured (reaction-norm) G×E | RN ≈ MT (both >> naive) | RN competitive only when G×E is low-rank |
| Factor-analytic rank analysis | rank predicts model choice | FA(2)=MT under RN; full MT best under CS (low rG); FA(1)=MT at high rG |
| 30% missing data | ranking unchanged | MT ≥ single-env > naive |

Real data show E1 genetically nearly independent of E2/E4/E5 (type-B rG ≈ 0.01–0.13),
exactly where cross-environment borrowing is hardest.

## Layout
- `analysis/` — `utils.py` (VanRaden G, GBLUP CV, variance components,
  Haseman–Elston, covariance bending), `simulator.py` (forward population,
  unstructured + reaction-norm G×E), `models.py` (single/pooled/MT/reaction-norm),
  and the six pipeline scripts.
- `output/` — JSON/CSV results + figures.
- `manuscript/manuscript.md` — full draft manuscript (abstract, intro, methods,
  results tables, discussion, references).
- `manuscript/verification_report.md` — data + model verification log
  (reproducibility proof for reviewers).
- `manuscript/cover_letter.md` — submission cover-letter template.
- `manuscript/journal_formatting_guide.md` — journal-specific submission guide.
- `output/figures_pub/` — publication-quality figures (300 DPI).

See `PROJECT_PLAN.md` for the full plan and candidate journals.
