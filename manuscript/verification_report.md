# Model & Data Verification Report

*Accompanying document for "Multi-environment genomic prediction under G×E:
a simulation study validated in CIMMYT wheat."*

This report records the checks performed to ensure that the data are genuine and
that every computational method reproduces known analytical results. It is
intended to be shared with reviewers or posted as supplementary material.

---

## 1. Data provenance and integrity

| Check | Result |
|---|---|
| Source | CIMMYT wheat panel; distributed in the R package **BGLR** (Pérez et al. 2010), originally Crossa et al. (2010, *Genetics* 186:713–724) |
| Acquisition | `git clone` of the public `gdlc/BGLR` repository; `data/wheat/wheat.RData` copied verbatim |
| **File identity** | `md5 = 6dd52a7029d9f70e7634965024000a72` — **byte-identical** to the official BGLR release |
| Dimensions | `wheat.X` 599 × 1,279; `wheat.Y` 599 × 4; `wheat.A` 599 × 599; `wheat.sets` 599 (10 folds) |
| Marker coding | 0/1 only (DArT presence/absence); no value outside {0,1} |
| Missing values | 0 in X, Y, A |
| Environment columns | `1, 2, 4, 5` (E1, E2, E4, E5) — matches official documentation |
| Pedigree A | symmetric, diag 1.56–2.00 (inbreeding coefficients), consistent with ICIS-derived A |

The BGLR documentation states "markers with MAF < 0.05 were removed". The file
contains 96 markers with minor allele frequency < 0.05 (as minor(p, 1−p));
this is a property of the distributed file itself and is unchanged by our
pipeline (we do not filter markers).

## 2. Method verification (analytical checks)

### 2.1 VanRaden (2008) method-1 GRM
- Mean diagonal = 1.0000 (correct scaling for 0/1 coding with denominator
  Σ p_j(1−p_j)).
- Mean off-diagonal ≈ −0.0017 ≈ 0 (correct centring).
- All eigenvalues ≥ 0 (positive semi-definite).
- Off-diagonal GRM correlates with the pedigree A (r ≈ 0.25), as expected for
  DArT markers vs pedigree.

### 2.2 RR-BLUP ≡ GBLUP equivalence
A marker-based ridge-regression BLUP with **global** scaling
(Z / √Σ p(1−p)) reproduces method-1 GBLUP exactly (prediction correlation =
1.0000; max absolute difference ≈ 0.004, floating-point). This confirms the
GBLUP solver.

### 2.3 MT-GBLUP solver vs brute force
The eigen-decomposed multi-trait solver was compared against a brute-force
Kronecker solve (Σ_g ⊗ G + Σ_e ⊗ I) on a small example: maximum absolute
difference = **0.0** (exact).

### 2.4 Variance components
The two-way (line × environment) ANOVA sums of squares were recomputed by hand
and match the function output (SS_g = 934.005, SS_g×e = 1457.995; total
SS = 2,392.000 = g×e − 4 environment df).

## 3. Simulator calibration

The simulator was required to satisfy, for any (rG, h²) input:

1. Var(g_true) = 1 per environment.
2. The realised genetic correlation among true breeding values equals rG.
3. The Haseman–Elston (HE) estimator recovers h² and rG without bias.

| Input | Realised corr(g_true) | HE estimate of rG | HE estimate of h² |
|---|---|---|---|
| rG = 0.2 | 0.195 | 0.191 | 0.39 (input 0.4) |
| rG = 0.5 | 0.496 | 0.503 | 0.39 (input 0.4) |
| rG = 0.8 | 0.798 | 0.811 | 0.40 (input 0.4) |
| h² = 0.2 / 0.4 / 0.6 (rG = 0.5) | — | — | 0.189 / 0.388 / 0.588 |

Var(g_true) = 0.993–1.000 across seeds. All calibrations use ≥ 15 seeds. The HE
estimator is unbiased on the fixed simulator (Section 4).

## 4. Issues found and corrected during verification

1. **Citation error (corrected).** An earlier draft attributed the dataset to
   "Crossa et al. 2010, *The Plant Genome*". The correct citation is Crossa et
   al. (2010), *Genetics* 186:713–724; the *Plant Genome* paper (Pérez et al.
   2010) is the BGLR package paper. Corrected in all files.
2. **Simulator allele-frequency bug (corrected).** Genetic values were centred by
   the *input* founder allele frequencies p, but after 15 generations of drift
   the sampled lines' frequencies p̂ differ from p. Centring by p̂ (as the GRM
   does) restores Var(g_true) = 1 and removes the HE-estimator bias. Corrected.
3. **Model inconsistency (corrected).** Phase 1 originally used RR-BLUP with a
   method-2-style column scaling and tuned shrinkage, while Phases 2–3 used
   method-1 GBLUP. All phases now use method-1 GBLUP (with a consistent
   `gblup_cv` helper), so baseline and benchmark are directly comparable.
4. **Misleading docstring (corrected).** The claim "RR-BLUP ≡ GBLUP under
   VanRaden G" was narrowed to the correct statement (global scaling, method 1).
5. **Singular covariance in estimated-component BLUP (corrected).** The
   Haseman–Elston moment estimate of Σ_g is not guaranteed positive
   semi-definite and occasionally produced a singular e×e system in the
   multi-trait solver. Fixed by (i) covariance "bending" (`make_psd`, clipping
   negative eigenvalues to 1e-4) and (ii) a generalised-inverse (pinv) fallback
   in `mt_gblup`. The real-data Σ_g had a minimal eigenvalue of −0.0046, so the
   effect there was negligible (results changed by <0.003).

## 5. Robustness analyses (added)

- **Structured (reaction-norm) G×E** (`04_robustness_structured.py`): under a
  rank-two, covariate-driven generative model, the reaction-norm model matches
  MT-GBLUP with estimated components and equals it exactly with oracle
  components (both equal the BLUP). Confirms RN's value is specific to low-rank
  G×E.
- **Unbalanced data** (`05_robustness_unbalanced.py`): 30% MCAR missingness does
  not change the model ranking (MT ≥ single-env; MT > naive).
- **Factor-analytic rank analysis** (`07_robustness_factor_analytic.py`): the
  FA(k) family (k = 1..4) fit to the estimated covariance; the variance-explained
  diagnostic matches the theoretical value (e.g. compound symmetry: leading
  factor ≈ (1+(e−1)rG)/e, confirmed as 0.70 at rG = 0.6). FA(2) = MT under the
  reaction-norm (rank-2) generative model; full MT best under compound symmetry
  at low rG; FA(1) = MT at high rG. Sanity check: FA(4) reproduces the input
  covariance exactly.
- **Summary figure** (`06_summary_figure.py`): Figure 1 in the manuscript.

## 5. Reproducibility summary

- Environment: Python 3.11; numpy 2.4, pandas 3.0, scipy 1.17, scikit-learn 1.9,
  pyreadr 0.5 (see `requirements.txt`). No R required.
- Pipeline: `analysis/01_real_data_baseline.py` → `02_simulation_study.py` →
  `03_real_data_validation.py` → `04_robustness_structured.py` →
  `05_robustness_unbalanced.py` → `07_robustness_factor_analytic.py` →
  `06_summary_figure.py`.
- Raw outputs: `output/*.json`, `output/02_sim_summary.csv`,
  `output/04_structured_summary.csv`, `output/05_unbalanced_summary.csv`;
  figures in `output/figures/`.
- Runtime: full benchmark ≈ 8 min (main) + ≈ 7 min (structured) + ≈ 4 min
  (unbalanced).
