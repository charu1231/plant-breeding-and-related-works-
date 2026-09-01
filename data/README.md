# Data dictionary — `data/wheat.RData`

## Provenance

- **Source:** CIMMYT Global Wheat Program; distributed with the R package **BGLR**
  (Pérez et al. 2010, *The Plant Genome* 3:106–116); originally published in
  Crossa et al. (2010), *Genetics* 186:713–724.
- **Acquisition:** `git clone https://github.com/gdlc/BGLR`; the file
  `data/wheat/wheat.RData` was copied verbatim.
- **Integrity:** `md5 = 6dd52a7029d9f70e7634965024000a72` — byte-identical to the
  official BGLR release (verified against the public repository).

## Objects

| Object | Shape | Content |
|---|---|---|
| `wheat.X` | 599 × 1,279 | DArT marker genotypes, coded **0/1** (marker presence/absence) |
| `wheat.Y` | 599 × 4 | Grain yield (standardised, mean 0, SD 1) in four mega-environments |
| `wheat.A` | 599 × 599 | Pedigree-derived additive relationship matrix (ICIS Browse) |
| `wheat.sets` | 599 × 1 | Predefined 10-fold cross-validation partition (values 1–10) |

## Environment columns (`wheat.Y`)

Columns are named `1`, `2`, `4`, `5` in the file — these correspond to four
CIMMYT target mega-environments, referred to as **E1, E2, E4, E5** in the
manuscript.

## Marker properties

- 0/1 coding (dominant DArT); no value outside {0, 1}.
- No missing genotypes (imputed by the data providers).
- Presence frequency `p_j = mean(X_j)`: median ≈ 0.57, range ≈ 0.008–0.987.
- Minor allele frequency `min(p, 1−p) < 0.05` for 96 of 1,279 markers (a
  property of the distributed file; the pipeline does not filter markers).

## Usage

```python
import pyreadr
res = pyreadr.read_r("data/wheat.RData")
X = res["wheat.X"].values.astype(float)   # 599 x 1279 markers (0/1)
Y = res["wheat.Y"].values.astype(float)   # 599 x 4 environments
A = res["wheat.A"].values.astype(float)   # 599 x 599 pedigree A
sets = res["wheat.sets"].values.ravel().astype(int)  # 10-fold partition
```
