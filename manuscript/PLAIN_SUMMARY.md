# Ek Page Mein — Aapka Research Paper

*(Simple language summary — paper submit karne se pehle ye padh lo)*

## Kya kiya hai?

Humne plant breeding (genomic selection) ka ek **research paper** banaya hai.
Simple baat: **ek hi variety alag-alag khet/mausam mein alag perform karti hai** —
ise bolte hain *genotype × environment interaction (G×E)*. Breeding mein asli
challenge ye hai ki ek environment mein achhi variety dusre mein fail ho sakti hai.

Paper ka sawaal: **kaunsa statistical model sabse achha prediction deta hai jab
multiple environments ka data ho?**

## Kya kiya (kaam)

1. **Asli data** — CIMMYT ka wheat dataset (599 lines, 1279 markers, 4 environments).
   Ye duniya ka famous public dataset hai (Crossa et al. 2010), bilkul genuine.
2. **Simulation** — computer par asli wheat jaisa data banaya, jisme G×E ki taakat
   (rG) aur heritability (h²) ko control kiya. Isse pata chalta hai kaunsa model
   kab behtar hai.
3. **4 models compare** kiye: single-environment GBLUP, simple average (naive),
   reaction-norm, aur multi-trait GBLUP (MT-GBLUP).

## Kya mila (results)

| Sawaal | Jawaab |
|---|---|
| Sabse achha model? | **MT-GBLUP** (environments ko correlated maanta hai) — har scenario mein sabse accurate |
| Ye kab zyada fayda deta hai? | Jab environments ke beech **genetic correlation (rG) zyada** ho, aur heritability **kam** ho |
| Reaction-norm model kab use karein? | Sirf jab G×E ka structure **low-rank / covariate-driven** ho — warna simple average se bhi kharab |
| Missing data (30%) se? | Ranking nahi badalti — result robust hai |
| Real wheat data par? | Simulation wahi ranking dikhata hai (MT-GBLUP 0.379 vs naive 0.313) |

**Bonus finding:** humne dikhaya ki **factor-analytic rank diagnostic** se pehle hi
pata lag jata hai kaunsa model best hoga — ye breeders ke liye practical rule hai.

## Kya ready hai (deliverables)

- **Manuscript** — `manuscript/manuscript.docx` (Word, editable) aur `.md`
- **7 figures** — `output/figures_pub/` (300 DPI, journal-ready)
- **Saara code** — `analysis/` (Python, reproducible) + `run_all.sh` (one-command)
- **Verification report** — proof ki data aur math sab sahi hai
- **GitHub PR** — https://github.com/charu1231/plant-breeding-and-related-works-/pull/1

## Aapko kya karna hai (sirf 3 kaam)

1. **Apna naam + affiliation + funding** daal do manuscript mein (abhi placeholder hai).
2. **Journal choose karo** — suggest: *Theoretical and Applied Genetics* (IF ~4.4).
3. **Submit** karo (cover letter ready hai).

Baaki sab — data, analysis, results, tables, figures, references — **ho chuka hai aur
verified hai**.
