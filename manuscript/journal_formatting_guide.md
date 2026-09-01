# Journal Formatting & Submission Guide

How to turn the current manuscript into a submission-ready package for the
target journal(s). All computational content is finished; what remains is
author-specific information and journal-specific cosmetics.

---

## 1. Target journals (with key facts)

| Journal | ~IF | Figure format | Max abstract | Style |
|---|---|---|---|---|
| **Theoretical and Applied Genetics (TAG)** | ~4.4 | TIFF/EPS ≥300 DPI, ~8.6 cm (1 col) / ~17.8 cm (2 col) | 250 words (structured) | Springer; numbered refs |
| **The Plant Genome** | ~4.1 | PNG/TIFF ≥300 DPI | 250 words | Wiley (ACSESS); author–date refs |
| **Frontiers in Plant Science** | ~4.1 | TIFF/JPG ≥300 DPI | ~350 words | Frontiers; author–date |
| **BMC Plant Biology** | ~4.2 | PNG/TIFF ≥300 DPI | 350 words (structured) | BMC; numbered refs |

*Recommendation:* **TAG** — the benchmark/methods-contribution framing fits its
scope best, and our structured abstract already matches its format.

## 2. What to fill in (author-side placeholders)

In `manuscript/manuscript.md` (or the editable `manuscript/manuscript.docx`):
- Title-page block: author names, affiliations, ORCID, corresponding-author email.
- `## 6. Declarations`: funding statement, author contributions, conflicts.
- Date in `manuscript/cover_letter.md`; suggested reviewers (replace examples).
- `CITATION.cff`: replace the author placeholder with the real author names.

## 3. Figures

Publication versions (300 DPI, journal font sizes) are pre-generated in
`output/figures_pub/`:

| File | Caption suggestion (for legend) |
|---|---|
| `fig1_simulation.png` | Fig. 1 — Accuracy vs rG (untested environment; 3 h² panels) |
| `fig2_simulation_lines.png` | Fig. 2 — Accuracy vs rG (untested lines; 3 h² panels) |
| `fig3_real_data.png` | Fig. 3 — Wheat type-B rG and validation accuracies |
| `fig4_robustness.png` | Fig. 4 — Structured G×E, FA rank, unbalanced data |
| `fig5_overview.png` | Fig. 5 (optional) — combined overview |

Working versions (160 DPI, for review/PDF) are in `output/figures/` (01–07).

## 4. Tables (already formatted, verified)

- Table 1 — Model definitions (Methods 2.3).
- Table 2 — Main simulation, Task A (Results 3.2).
- Table 3 — Main simulation, Task B (Results 3.3).
- Table 4 — Real-data type-B rG + validation (Results 3.4).
- Table 5 — Structured G×E (Results 3.5).
- Table 6 — Unbalanced data (Results 3.6).
- Table 7 — Factor-analytic rank analysis (Results 3.7).

All table cells are machine-verified against `output/*.csv/json` by
`analysis/check_all_tables.py`.

## 5. Supplementary material to upload

- `manuscript/verification_report.md` — data/method verification log.
- `analysis/` — full reproducible Python pipeline (or a Zenodo/figshare DOI).
- `data/wheat.RData` — dataset (public; cite Crossa et al. 2010).
- `output/02_sim_results.json`, `04_structured_results.json`,
  `05_unbalanced_results.json`, `07_fa_results.json` — per-replicate raw results.

## 6. Submission checklist

- [ ] Authors/affiliations/ORCID added.
- [ ] Funding + author contributions + conflict-of-interest statements.
- [ ] Abstract ≤ journal word limit (ours = 250 words).
- [ ] References in target journal style (currently numbered, TAG-compatible;
      convert to author–date for Wiley/Frontiers).
- [ ] Figures at ≥300 DPI in required format (done, `output/figures_pub/`).
- [ ] Cover letter personalised (`manuscript/cover_letter.md`).
- [ ] Suggested reviewers chosen.
- [ ] Data/code availability statement (already drafted).
- [ ] Run `analysis/check_all_tables.py` once more after any text edit.
