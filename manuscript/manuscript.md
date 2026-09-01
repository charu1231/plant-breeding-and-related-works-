# Multi-environment genomic prediction under genotype × environment interaction: a simulation study validated in CIMMYT wheat

*Author(s): [Name(s), affiliations, ORCID]*

*Target journal: Theoretical and Applied Genetics (primary); also suitable for The
Plant Genome, Frontiers in Plant Science, BMC Plant Biology.*

---

## Abstract

Genomic prediction in plant breeding must contend with genotype × environment
interaction (G×E): the genetic value of a line is not constant across
environments. Multi-trait (MT) models that treat environments as correlated
traits exploit the genetic covariance between environments, but their benefit
relative to simpler approaches is not well characterised. We benchmarked
single-environment GBLUP, the multi-environment phenotypic mean, a reaction-norm
model, and MT-GBLUP in a forward simulation calibrated to a real wheat panel
(599 lines, 1,279 markers, four environments), with G×E controlled through the
type-B genetic correlation (rG = 0.2–0.8) at heritabilities h² = 0.2–0.6.
MT-GBLUP was consistently the most accurate method: for untested lines its
advantage over single-environment GBLUP grew with rG (up to 0.10) and was
largest at low heritability; for an untested environment it beat the naive mean
by 0.02–0.12. The reaction-norm model added little under compound-symmetry G×E
and could be worse than the naive mean at low rG. A factor-analytic rank
analysis showed model choice is predictable from the estimated covariance — two
factors sufficed under reaction-norm G×E (FA(2) = MT), whereas the variance was
more spread under compound symmetry, favouring the full model (or FA(1) at high
rG) — and the model ranking was robust to 30% missing data. Validation in the
CIMMYT wheat panel (Crossa et al. 2010) showed one environment genetically
nearly independent of the others (type-B rG ≈ 0.01–0.13 vs 0.84–0.97), and
MT-GBLUP improved prediction of the untested environment from r = 0.313 to
0.379. These results guide multi-environment model choice.

**Keywords:** genomic selection · GBLUP · genotype × environment interaction ·
multi-trait model · reaction norm · factor-analytic · wheat · simulation

---

## 1. Introduction

Genomic selection has transformed plant breeding by enabling the prediction of
genetic merit from genome-wide marker data (Meuwissen et al. 2001; Habier et al.
2007). The workhorse method, genomic best linear unbiased prediction (GBLUP),
regresses phenotypes on a marker-derived relationship matrix and has been
validated across crops (Crossa et al. 2010; VanRaden 2008). Most applications,
however, are single-trait and single-environment: a model is trained and
evaluated within one environment, which is rarely how breeding programmes
actually operate.

Breeding programmes evaluate candidate lines in multi-environment trials (MET)
precisely because the ranking of lines changes across environments — the
phenomenon of genotype × environment interaction (G×E). Statistically, G×E
means that the genetic correlation between the same trait measured in two
environments (the type-B genetic correlation, rG) is less than one (Falconer &
Mackay 1996). When rG is low, a line selected in one environment may perform
poorly in another, and prediction accuracy for untested environments decays.
Consequently, a central question is how best to model and exploit the genetic
covariance among environments when making predictions.

Several model families address this question. The simplest approach — averaging
phenotypes across environments, or a single line main effect — treats
environment effects as noise and borrows information only through a common mean.
Multi-trait GBLUP instead treats the trait in each environment as a correlated
trait and models the full genetic covariance matrix (Burgueño et al. 2012;
Crossa et al. 2016). Reaction-norm (random regression) models express the
genetic value as a function of an environmental covariate, while factor-analytic
models impose low-rank structure on the environment covariance matrix (Jarquín
et al. 2014; Lopez-Cruz et al. 2015; Lado et al. 2016). Although these methods
have each been shown to outperform naive approaches in specific data sets, a
systematic benchmark across a controlled G×E gradient, with methods placed on an
equal footing, remains scarce. In particular, the conditions under which the
added complexity of multi-trait or reaction-norm models actually pays off —
relative to a simple multi-environment mean or single-environment GBLUP — are
not well delineated.

The objectives of this study were therefore: (i) to benchmark four
multi-environment prediction models under a factorial simulation in which G×E
magnitude (rG) and heritability (h²) were explicitly controlled; (ii) to
validate the conclusions on a public CIMMYT wheat panel; and (iii) to translate
the results into practical guidance for breeders. We used a forward simulation
calibrated to the real wheat panel (marker count, allele-frequency spectrum,
number of lines, and number of environments) so that the simulation results are
directly interpretable alongside the real-data results. All models used the same
genomic relationship matrix and, in the simulation, the same (oracle) variance
components, so that differences in accuracy reflect model structure rather than
estimation artefacts.

## 2. Materials and methods

### 2.1 Real data

We used the CIMMYT wheat panel distributed with the BGLR package (Pérez et al.
2010), originally from Crossa et al. (2010): 599 historical CIMMYT wheat lines
genotyped with 1,279 DArT markers (coded 0/1 for presence/absence) and
phenotyped for grain yield in four target mega-environments, here denoted E1,
E2, E4, E5 following the original files. The data contain no missing values, and
a pedigree-derived additive relationship matrix (A) is also provided. All
analyses used the VanRaden (2008) method-1 genomic relationship matrix
**G = ZZ′/Σ p_j(1−p_j)**, where **Z** is the centred marker matrix and p_j is the
sample allele frequency at marker j.

### 2.2 Simulation

We simulated a population that mimics the wheat panel (all details in the
accompanying code). In brief, marker allele frequencies were resampled from the
real panel's frequency spectrum; 100 diploid founders were generated; the
population was random-mated for 15 generations with recombination (Poisson
crossovers on 21 chromosomes); and 599 doubled-haploid (DH) inbred lines were
sampled from the final generation. This produced related lines with realistic
linkage disequilibrium and 0/1 (homozygous) marker coding. Genetic values were
generated under an infinitesimal model in which all markers are causal.

G×E was parameterised by the type-B genetic correlation rG through a
compound-symmetry factor model:

```
g_ij = √rG · c_i + √(1−rG) · t_ij,
```

where c_i is the main (marker-based) genetic value of line i and t_ij is an
environment-specific, marker-based interaction term; both have covariance equal
to **G**. This gives an oracle genetic covariance
**Σ_g = (1−rG) I + rG J** across environments (unit genetic variance per
environment) and a genetic correlation rG between any pair of environments.
Residuals were drawn as e_ij ~ N(0, σ²_e) with σ²_e = (1−h²)/h², so that the
per-environment heritability equals h². The scenario grid was
rG ∈ {0.2, 0.4, 0.6, 0.8} × h² ∈ {0.2, 0.4, 0.6}, with 15 independent
replicates per scenario.

**Structured (reaction-norm) G×E.** To delineate the conditions under which
reduced-rank reaction-norm models are appropriate, we also simulated a
covariate-driven model g_ij = c_i + b_i z_j, where c_i (intercept) and b_i
(slope) are independent marker-based breeding values (both with covariance ∝ G)
and z is an observable environmental covariate that also drives the environment
mean (so it is readable from the environment-mean phenotype). The genetic
covariance is therefore rank-two: Σ_g = σ²_c J + σ²_b z z′. G×E strength is set
by κ = σ²_b/σ²_c, with the average genetic variance scaled to 1. We used
κ ∈ {0.3, 0.8} × h² ∈ {0.3, 0.6}, 20 replicates; here all models used variance
components **estimated** from data (Section 2.5), matching practice, with the
oracle-component versions computed as references.

**Unbalanced (missing) data.** As a screening check we masked 30% of phenotype
cells completely at random (MCAR) in the compound-symmetry simulation
(rG ∈ {0.3, 0.7}, h² = 0.4, 20 replicates). Single-environment GBLUP and the
naive mean used the observed cells natively, while MT-GBLUP was fitted on
environment-mean-imputed data.

### 2.3 Models benchmarked

| Model | Description |
|---|---|
| Single-env GBLUP | per-environment GBLUP with λ = σ²_e / σ²_g |
| Naive mean | unweighted phenotypic mean across training environments (baseline) |
| RN-GBLUP | rank-2 linear reaction norm on the environment mean (a constrained MT model) |
| FA(k)-GBLUP | factor-analytic rank-k MT model (ΛΛ′ + Ψ; k = 1..e, with k = e ≡ MT) |
| MT-GBLUP | environments as correlated traits; Cov(vec **U**) = **Σ_g** ⊗ **G** |

MT-GBLUP was solved exactly by eigen-decomposing **G** and solving e×e systems
per eigenvector. In the simulation all models used the oracle **Σ_g** and σ²_e;
on real data **Σ_g** and σ²_e were estimated (Section 2.5). The reaction-norm
model constrains the environment covariance to the span of {1, z}, where z is an
environmental covariate (the environment mean, in both simulation and real
data), implemented as **Σ_rr = P_z Σ_g P_z**. The factor-analytic FA(k) models
re-express the (estimated) covariance as **Σ_fa = ΛΛ′ + Ψ**, where Λ is the
e×k matrix of the leading k loadings and Ψ is diagonal (environment-specific
variances); FA(e) is the full model, so the FA family is nested within MT-GBLUP
(Burgueño et al. 2012). For FA models the rank k was chosen by the proportion of
genetic variance explained by the leading k factors of the estimated covariance.

### 2.4 Validation schemes

Two prediction tasks were considered. **Task A (CV1, new environment):** predict
the genetic value of all lines in a held-out environment from the other three.
**Task B (CV2, new lines):** predict the genetic value of untested lines, using
5-fold cross-validation in the simulation and the predefined 10-fold
`wheat.sets` partition on real data. Accuracy was measured as the Pearson
correlation between predicted and true breeding values (simulation) or between
predicted values and the observed phenotype (real data, where the phenotype is
the observable proxy for the breeding value and the accuracy is therefore scaled
by √h²).

### 2.5 Real-data variance components

We estimated the environment genetic covariance **Σ_g** and residual variances
by the multivariate Haseman–Elston moment estimator (Haseman & Elston 1972; Lee
et al. 2012): the genetic covariance between environments j and k is the
regression slope of the empirical cross-product of centred phenotypes on the
off-diagonal of **G**. Moment estimates need not be positive semi-definite in
finite samples, so we applied covariance "bending" (clipping negative
eigenvalues to a small positive value) before use. This estimator was validated
on simulated data before use (Section 3.3 and the verification report).

### 2.6 Software and reproducibility

All analyses were implemented in Python 3 (numpy/scipy/pandas; no R required).
Code, data-processing steps, and all figure/table outputs are provided
(`analysis/` and `output/`). The dataset file is byte-identical to the official
BGLR release. A model-verification report accompanies this manuscript. Figure 1
(summary) is produced by `analysis/06_summary_figure.py`; the robustness
analyses are `analysis/04_robustness_structured.py` and
`analysis/05_robustness_unbalanced.py`.

## 3. Results

### 3.1 Real-data structure (Phase 1)

The two-way (line × environment) analysis of the wheat panel (Comstock & Moll
1963) gave a broad-sense heritability on entry means of H² ≈ 0.48. The four environments showed markedly
heterogeneous genetic relationships: E1 was phenotypically uncorrelated or
negatively correlated with the others (r = −0.19 to −0.02), whereas E2, E4 and
E5 were moderately to strongly correlated (r = 0.39–0.66). Single-environment
GBLUP cross-validation accuracies ranged from r = 0.39 (E4) to r = 0.53 (E1),
consistent with this G×E structure.

### 3.2 Simulation — Task A (predict new environment)

Mean prediction accuracy versus true breeding value (15 replicates):

| rG | h² | Naive mean | MT-GBLUP | RN-GBLUP |
|----|----|------------|----------|----------|
| 0.2 | 0.2 | 0.147 | **0.191** | 0.145 |
| 0.2 | 0.4 | 0.202 | **0.231** | 0.170 |
| 0.2 | 0.6 | 0.241 | **0.256** | 0.184 |
| 0.4 | 0.2 | 0.287 | **0.367** | 0.331 |
| 0.4 | 0.4 | 0.383 | **0.430** | 0.372 |
| 0.4 | 0.6 | 0.444 | **0.467** | 0.389 |
| 0.6 | 0.2 | 0.417 | **0.520** | 0.502 |
| 0.6 | 0.4 | 0.542 | **0.598** | 0.565 |
| 0.6 | 0.6 | 0.616 | **0.642** | 0.590 |
| 0.8 | 0.2 | 0.538 | **0.656** | 0.651 |
| 0.8 | 0.4 | 0.685 | **0.745** | 0.735 |
| 0.8 | 0.6 | 0.768 | **0.794** | 0.777 |

MT-GBLUP was the most accurate method in every scenario. The MT−naive gap grew
with rG (from ≈0.02–0.05 at rG = 0.2 to ≈0.03–0.12 at rG = 0.8). The
reaction-norm model matched the naive mean at low rG (and could be slightly
worse, e.g. 0.170 vs 0.202 at rG = 0.2, h² = 0.4) and approached MT-GBLUP at
high rG, but never exceeded it.

### 3.3 Simulation — Task B (predict new lines)

| rG | h² | Single-env GBLUP | MT-GBLUP |
|----|----|------------------|----------|
| 0.2 | 0.2 | 0.514 | **0.519** |
| 0.2 | 0.4 | 0.628 | **0.632** |
| 0.2 | 0.6 | 0.704 | **0.707** |
| 0.4 | 0.2 | 0.515 | **0.540** |
| 0.4 | 0.4 | 0.630 | **0.647** |
| 0.4 | 0.6 | 0.706 | **0.717** |
| 0.6 | 0.2 | 0.517 | **0.572** |
| 0.6 | 0.4 | 0.632 | **0.669** |
| 0.6 | 0.6 | 0.707 | **0.731** |
| 0.8 | 0.2 | 0.517 | **0.613** |
| 0.8 | 0.4 | 0.633 | **0.699** |
| 0.8 | 0.6 | 0.708 | **0.753** |

MT-GBLUP was at least as accurate as single-environment GBLUP everywhere; the
advantage was negligible at rG = 0.2 (≈0.003–0.005) and grew to 0.04–0.10 at
rG = 0.8, with the relative gain largest at low heritability.

### 3.4 Real-data validation (Phase 3)

The Haseman–Elston estimator was first validated on simulated data, where it
recovered the true h² and rG without bias (e.g. estimated h² = 0.19/0.39/0.59 at
true 0.2/0.4/0.6; estimated rG = 0.19/0.50/0.81 at true 0.2/0.5/0.8). Applied to
the wheat panel, it gave the estimated type-B genetic correlation matrix:

| | E1 | E2 | E4 | E5 |
|--|----|----|----|----|
| E1 | 1.00 | 0.13 | 0.08 | 0.01 |
| E2 | | 1.00 | 0.97 | 0.84 |
| E4 | | | 1.00 | 0.94 |
| E5 | | | | 1.00 |

E1 is genetically nearly independent of E2/E4/E5 — matching the published
description of this panel (Crossa et al. 2010) — while E2, E4 and E5 are
strongly correlated. On Task A (predict an untested environment), MT-GBLUP
achieved r = 0.379 versus r = 0.313 for the naive mean and r = 0.262 for the
reaction-norm model. On Task B (predict untested lines, 10-fold), MT-GBLUP
(r = 0.443) matched single-environment GBLUP (r = 0.440); the MT advantage was
concentrated in the most-correlated environment E4 (0.410 vs 0.379). Thus the
real data reproduce the simulation ranking: MT-GBLUP ≥ naive ≥ reaction-norm on
Task A, and MT-GBLUP ≈ single-env with a small edge on Task B.

### 3.5 Structured (reaction-norm) G×E

Under the reaction-norm generative model (Section 2.2), with variance
components estimated from data, the reduced-rank reaction-norm model matched
MT-GBLUP in every scenario (Table, Task A): at κ = 0.3, RN 0.719/0.846 vs MT
0.717/0.846 (h² = 0.3/0.6); at κ = 0.8, RN 0.690/0.834 vs MT 0.688/0.833. Both
models far outperformed the naive mean (0.346–0.670), and with oracle
components the two were identical to three decimals (e.g. 0.725 = 0.725), as
expected because the true covariance is already in the reaction-norm subspace.
On Task B the same pattern held (RN ≈ MT > single-env). Notably, at κ = 0.8 the
realised genetic correlations spanned negative to strongly positive values
(−0.28 to 0.89), yet the reaction-norm model remained as accurate as the full
multi-trait model, because the rank-two structure captured the genuine pattern
of G×E. These results contrast with the compound-symmetry setting (Section 3.2),
where the reaction-norm model was clearly inferior, and identify low-rank,
covariate-driven G×E as the regime in which reaction-norm models are
competitive.

| κ | h² | Naive | MT-GBLUP | RN-GBLUP | MT (oracle) | RN (oracle) |
|---|---|-------|----------|----------|-------------|-------------|
| 0.3 | 0.3 | 0.538 | 0.717 | **0.719** | 0.725 | 0.725 |
| 0.3 | 0.6 | 0.670 | 0.846 | **0.846** | 0.855 | 0.855 |
| 0.8 | 0.3 | 0.346 | 0.688 | **0.690** | 0.698 | 0.698 |
| 0.8 | 0.6 | 0.452 | 0.834 | **0.834** | 0.844 | 0.844 |

### 3.6 Unbalanced data (30% missing)

With 30% MCAR missingness, the qualitative conclusions were unchanged
(Table). On Task A, MT-GBLUP (imputed) outperformed the naive mean at both
rG = 0.3 (0.310 vs 0.238) and rG = 0.7 (0.629 vs 0.515); on Task B it matched
single-environment GBLUP at rG = 0.3 (0.574 vs 0.575) and exceeded it at
rG = 0.7 (0.628 vs 0.575). Thus moderate missing data, handled by simple
environment-mean imputation, does not upset the model ranking.

| rG | Task A: naive | Task A: MT | Task B: single | Task B: MT |
|----|---------------|------------|----------------|------------|
| 0.3 | 0.238 | **0.310** | 0.575 | 0.574 |
| 0.7 | 0.515 | **0.629** | 0.575 | **0.628** |

### 3.7 Factor-analytic rank analysis

The rank of the estimated environment covariance predicted model performance
across both G×E architectures (Figure 1d; the proportion of genetic variance
explained by the leading k factors, FA1/FA2/FA3, is reported). Under the
reaction-norm generative model, two factors explained ≈100% of the genetic
variance (0.996 at κ = 0.3 and κ = 0.8), and FA(2), FA(3), RN and the full MT
model all attained identical accuracy on Task A (0.843 at κ = 0.3, 0.832 at
κ = 0.8), while FA(1) was markedly worse (0.712 / 0.485) and the naive mean
worse still (0.665 / 0.446). Under compound-symmetry G×E the variance was
spread over more factors: at rG = 0.3 the leading factor explained only ≈0.50
of the genetic variance (three factors were needed for ≈0.90), and no
reduced-rank model beat the full MT-GBLUP (0.318 vs FA1 0.304, FA2 0.296,
FA3 0.295), while the reaction-norm model was worst (0.271). At rG = 0.6 the
leading factor explained ≈0.70, and FA(1) matched MT-GBLUP (0.590 vs 0.589),
outperforming FA(2)/FA(3) (0.573/0.573) — consistent with compound symmetry
being an FA(1) structure (one common factor plus environment-specific
variances), so that adding spurious factors is not helpful. These results
support the practical rule that the factor-analytic rank diagnostic — the
proportion of variance carried by the leading factors of the estimated
covariance — can guide model choice: a small number of dominant factors
favours the corresponding FA(k) model, whereas a spread spectrum favours the
full multi-trait model, and the reaction-norm model is only appropriate when a
measured environmental covariate drives G×E.

| Scenario | varExp(1/2/3) | naive | FA1 | FA2 | FA3 | RN | MT |
|---|---|---|---|---|---|---|---|
| CS, rG = 0.3 | 0.50/0.75/0.90 | 0.297 | 0.304 | 0.296 | 0.295 | 0.271 | **0.318** |
| CS, rG = 0.6 | 0.70/0.86/0.95 | 0.540 | **0.590** | 0.573 | 0.573 | 0.556 | 0.589 |
| RN, κ = 0.3 | 0.75/1.00/1.00 | 0.665 | 0.712 | **0.843** | 0.843 | 0.843 | 0.843 |
| RN, κ = 0.8 | 0.57/1.00/1.00 | 0.446 | 0.485 | **0.832** | 0.832 | 0.832 | 0.832 |

(Task A accuracy, new environment, 10 replicates; CS = compound symmetry,
RN = reaction norm. The rank diagnostic is the fraction of the estimated
genetic covariance explained by the leading 1/2/3 factors.)

## 4. Discussion

**Why MT-GBLUP wins.** Multi-trait GBLUP is the natural generative model for
correlated environments: it weights each training environment by its genetic
covariance with the target environment and shrinks the estimate through **G**.
When environments are strongly correlated (high rG), the same genetic signal
reappears in several environments, and MT-GBLUP pools it efficiently — which is
why its advantage over single-environment GBLUP grows with rG. When rG is low,
there is little shared signal to pool, and MT-GBLUP correctly reduces to
near-single-environment behaviour. The naive mean cannot make either adjustment:
it gives every environment equal weight regardless of its genetic relevance, so
it is systematically worse than MT-GBLUP for predicting a new environment.

**The reaction-norm model: when it helps and when it does not.** The two G×E
architectures we simulated frame the value of reduced-rank models precisely.
Compound-symmetry G×E (Σ_g = (1−rG)I + rG J) is an FA(1) structure: a single
common factor plus environment-specific variances (1−rG) along the diagonal. The
reaction-norm model discards those environment-specific variances and
constrains the remaining covariance to the span of {1, z}; because z (the
environment mean) carries no genetic information in this setting and the
specific variances are important — especially at low rG — the reaction-norm
model can be worse than even the naive mean. The factor-analytic family makes
this transparent: FA(1), which retains the specific variances, beats the naive
mean and, at rG = 0.6, matches the full MT-GBLUP; adding factors beyond the true
structure (FA(2), FA(3)) does not help. Under covariate-driven (reaction-norm)
G×E the true covariance is genuinely rank-two, and the reaction-norm model
matched MT-GBLUP exactly with oracle components and with estimated components
matched the corresponding FA(2) model (Jarquín et al. 2014). The practical
corollary is that the factor-analytic rank diagnostic — the proportion of the
estimated genetic covariance carried by its leading factors — is a useful and
cheap guide to model choice before fitting many candidate models.

**Agreement between simulation and real data.** The real wheat panel behaved
exactly as the simulation predicts. E1 is genetically nearly independent of the
other environments, and it is precisely there that cross-environment borrowing
fails; the correlated block {E2, E4, E5} is where MT-GBLUP earns its advantage.
The qualitative ranking of methods was identical in simulation and real data,
supporting the generality of the conclusions beyond the specific panel. The
robustness analyses reinforce this: the ranking was unchanged under 30% missing
data, and the structured-G×E results show that the reaction-norm model's
inferiority is specific to compound-symmetry G×E rather than to reduced-rank
models per se.

**Limitations.** (i) With one replicate per line × environment cell, σ²_g×e and
σ²_e are confounded in the ANOVA-based heritability; the marker-based
Haseman–Elston estimator separates them but relies on the GRM being a faithful
relationship measure. (ii) The simulation used oracle variance components so
that differences reflect model structure only; in practice these must be
estimated, and estimation error will reduce (but not reorder) accuracies.
(iii) Real-data Task A accuracy is reported against the observed phenotype and
is therefore scaled by √h² relative to the simulation's true-breeding-value
metric. (iv) We used a single crop panel, DArT (0/1) markers, and additive-only
models; the DArT coding and the finite marker panel influence the GRM. (v) The
reaction-norm model used only the environment mean as covariate; richer
environmental covariates could improve it.

**Practical recommendations.** When a multi-environment training set exists,
MT-GBLUP should be the default model: it is never worse than single-environment
GBLUP and its advantage grows with the genetic correlation among environments.
Breeders should estimate the type-B genetic correlation to decide whether MET
data can help a new target environment: if rG with the target is low, additional
phenotyping in that target environment is more valuable than borrowing from
dissimilar environments. Reaction-norm models should be reserved for settings
with low-rank, covariate-driven G×E. Future work should extend the benchmark to
larger environmental panels, factor-analytic rank selection, explicitly
unbalanced designs beyond MCAR missingness, and higher-order (multi-trait ×
multi-environment) models.

## 5. Conclusions

Across a controlled simulation gradient and a real CIMMYT wheat panel,
multi-trait GBLUP was the most accurate multi-environment genomic-prediction
method. Its advantage over single-environment GBLUP and the naive
multi-environment mean grows with the genetic correlation among environments and
is largest at low heritability. A reduced-rank reaction-norm model is inferior
under compound-symmetry G×E but matches MT-GBLUP when G×E is genuinely low-rank
and covariate-driven. A factor-analytic rank diagnostic predicts which regime
applies and hence which model to use. The model ranking is robust to moderate
missing data. These results give breeders a simple, evidence-based rule for
model choice in multi-environment genomic selection.

---

## 6. Declarations

**Data availability.** The wheat dataset is public (BGLR R package; Crossa et
al. 2010) and is provided in `data/wheat.RData` (byte-identical to the official
release). **Code availability.** Full reproducible pipeline in `analysis/`
(Python 3; dependencies in `requirements.txt`). **Conflicts of interest.** None
declared. **Funding.** [to be completed]. **Author contributions.** [to be
completed].

## 7. References

1. Meuwissen THE, Hayes BJ, Goddard ME (2001) Prediction of total genetic value
   using genome-wide dense marker maps. Genetics 157:1819–1829.
2. Habier D, Fernando RL, Dekkers JCM (2007) The impact of genetic relationship
   information on genome-assisted breeding values. Genetics 177:2389–2397.
   doi:10.1534/genetics.107.081190
3. VanRaden PM (2008) Efficient methods to compute genomic predictions. J Dairy
   Sci 91:4414–4423. doi:10.3168/jds.2007-0980
4. Crossa J, de los Campos G, Pérez P, Gianola D, Burgueño J, Araus JL, et al.
   (2010) Prediction of genetic values of quantitative traits in plant breeding
   using pedigree and molecular markers. Genetics 186:713–724.
   doi:10.1534/genetics.110.118521
5. Pérez P, de los Campos G, Crossa J, Gianola D (2010) Genomic-enabled
   prediction based on molecular markers and pedigree using the Bayesian Linear
   Regression Package in R. Plant Genome 3:106–116.
   doi:10.3835/plantgenome2010.04.0005
6. Burgueño J, de los Campos G, Weigel K, Crossa J (2012) Genomic prediction of
   breeding values when modeling genotype × environment interaction using
   pedigree and dense molecular markers. Crop Sci 52:707–719.
   doi:10.2135/cropsci2011.06.0299
7. Jarquín D, Crossa J, Lacaze X, du Cheyron P, Daucourt J, Lorgeou J, et al.
   (2014) A reaction norm model for genomic selection using high-dimensional
   genomic and environmental data. Theor Appl Genet 127:595–607.
   doi:10.1007/s00122-013-2243-1
8. Lopez-Cruz M, Crossa J, Bonnett D, Dreisigacker S, Poland J, Jannink JL, et
   al. (2015) Increased prediction accuracy in wheat breeding trials using a
   marker × environment interaction genomic selection model. G3 5:569–582.
   doi:10.1534/g3.114.016097
9. Crossa J, Pérez-Rodríguez P, Cuevas J, Montesinos-López OAF, Jarquín D, de
   los Campos G, et al. (2016) A genomic Bayesian multi-trait and
   multi-environment model. G3 6:2725–2744. doi:10.1534/g3.116.032359
10. Lado B, Barrios PG, Quincke M, Silva P, Gutiérrez L (2016) Modeling genotype
    × environment interaction for genomic selection with unbalanced data from a
    wheat breeding program. Crop Sci 56:2165–2179. doi:10.2135/cropsci2015.04.0207
11. Haseman JK, Elston RC (1972) The investigation of linkage between a
    quantitative trait and a marker locus. Behav Genet 2:3–19.
    doi:10.1007/BF01066731
12. Lee SH, Yang J, Goddard ME, Visscher PM, Wray NR (2012) Estimation of
    pleiotropy between complex diseases using single-nucleotide
    polymorphism-derived genomic relationships and restricted maximum likelihood.
    Bioinformatics 28:2540–2542. doi:10.1093/bioinformatics/bts474
13. Comstock RE, Moll RH (1963) Genotype-environment interactions. In: Statistical
    Genetics and Plant Breeding, NAS–NRC Publ. 982:164–196.
14. Falconer DS, Mackay TFC (1996) Introduction to Quantitative Genetics, 4th edn.
    Longman, Harlow.
