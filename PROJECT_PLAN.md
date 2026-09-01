# Project Plan — Plant Breeding / Genomic Prediction Research Paper

**Goal:** ek publishable research paper, **simulation + real public data** dono se,
target = impact-factor journal (TAG / The Plant Genome / Frontiers in Plant Science tier).

**Working title (draft):**
> *Benchmarking multi-environment genomic prediction under genotype × environment
> interaction: a simulation study validated on a CIMMYT wheat panel*

---

## 1. Dataset (real, public, citable)

- **Source:** BGLR R-package `wheat` dataset → originally from **Crossa et al. (2010)**
  *Genetics* 186:713–724 (CIMMYT global wheat program).
- **Content:** 599 historical CIMMYT wheat lines, 1279 DArT markers (0/1), grain yield
  in **4 environments**, 0 missing values, plus a pedigree (A) matrix.
- Location: `data/wheat.RData` (copied from the public `gdlc/BGLR` GitHub repo).

## 2. Analyses

### Phase 1 — Real-data baseline ✅ (done, `analysis/01_real_data_baseline.py`)
- G×E structure: environment correlation matrix, two-way (line × env) variance
  components, broad-sense heritability (H²).
- GBLUP per-environment 5-fold cross-validation prediction accuracy.

**Results so far:**
- H² ≈ 0.48 on entry means.
- Env 1 is uncorrelated/negatively correlated with envs 2/4/5 → clear G×E signal.
- Per-env CV accuracy r ≈ 0.37–0.52.

### Phase 2 — Simulation study (next)
- Build a numpy quantitative-genetics G×E simulator (QTL + polygenic background +
  environment main effects + G×E), **calibrated to the wheat data**
  (same n, marker count/MAF structure, heritability, G×E correlation pattern).
- Compare prediction models across scenarios:
  1. **Single-env GBLUP** (baseline)
  2. **Across-env / pooled GBLUP**
  3. **Multi-trait GBLUP** (environments as correlated traits; MTM/MTCV)
  4. **Reaction-norm (random regression) model** on an environmental covariate
  5. **Machine learning** (random forest / gradient boosting; optional NN)
- Scenario grid: G×E magnitude (env correlation ρ = 0.2 / 0.5 / 0.8),
  heritability (0.2 / 0.4 / 0.6), training size, marker density.
- Metrics: prediction accuracy (r), mean-squared error, ranking ability.

### Phase 3 — Real-data validation
- Apply the same model comparison to the wheat dataset; check whether simulation
  conclusions generalise. (Simulation = mechanism; real data = evidence.)

## 3. Deliverables
- Fully reproducible Python pipeline (`analysis/`).
- All figures in `output/figures/`.
- `manuscript/manuscript.md` — draft skeleton: abstract, methods, results tables, discussion.

## 6. Status (all phases complete ✅)
- [x] Phase 1 — real-data baseline: H² ≈ 0.48; E1 uncorrelated with E2/E4/E5; GBLUP CV r = 0.37–0.52.
- [x] Phase 2 — G×E simulation (forward population, heritable interaction): MT-GBLUP best across
  all rG × h² scenarios; gain grows with rG.
- [x] Phase 3 — real-data validation: type-B rG estimated (E1 ≈ 0.01–0.13 vs E2/E4/E5 = 0.84–1.0);
  MT-GBLUP 0.381 vs naive 0.313 (new env); 0.443 vs 0.440 (new lines).
- [x] Full manuscript with prose (`manuscript/manuscript.md`).
- [x] Verification report (`manuscript/verification_report.md`) — data md5,
  model-vs-brute-force, RR-BLUP≡GBLUP, simulator calibration, HE unbiasedness.
- [x] Cover letter template (`manuscript/cover_letter.md`).
- [x] Robustness — structured (reaction-norm) G×E (`04_robustness_structured.py`):
  RN ≈ MT under low-rank G×E (with estimated components).
- [x] Robustness — unbalanced data, 30% MCAR (`05_robustness_unbalanced.py`):
  model ranking unchanged.
- [x] Summary figure (`06_summary_figure.py` → `output/figures/06_summary.png`).
- [x] Factor-analytic rank analysis (`07_robustness_factor_analytic.py`):
  rank diagnostic predicts model choice; FA(2)=MT under reaction-norm G×E,
  full MT best under compound symmetry (low rG), FA(1)=MT at high rG.

### Remaining for submission (author-side, not computational)
- Fill in author names, affiliations, funding, and author contributions.
- Final journal formatting pass (target: TAG / The Plant Genome / Frontiers).
- Optional further extensions: factor-analytic rank selection, explicitly
  unbalanced designs (beyond MCAR), multi-trait × multi-environment models,
  more replicates.
- Submit + handle review.

## 4. Candidate journals (with impact factor)
| Journal | ~IF |
|---|---|
| Theoretical and Applied Genetics (TAG) | ~4.4 |
| The Plant Genome | ~4.1 |
| Frontiers in Plant Science | ~4.1 |
| BMC Plant Biology | ~4.2 |
| Heredity | ~3.5 |
| Euphytica / Crop Science | ~1.6–2.0 |

## 5. Environment notes (this sandbox)
- Only **PyPI** and **GitHub** are reachable (CRAN/R, dataverse, figshare, CyVerse blocked),
  so everything is **pure Python** (no R needed) and data comes from GitHub.
- Dependencies: `requirements.txt` (numpy, pandas, scipy, scikit-learn, matplotlib,
  seaborn, pyreadr, statsmodels).
- Use `python -m venv .venv && pip install -r requirements.txt` to recreate.
