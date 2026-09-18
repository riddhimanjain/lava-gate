# Denominator selection bias in local genetic correlation

Code, results and manuscript for:

> **Significance filtering induces denominator selection bias in local genetic
> correlation: closed-form characterisation, a gate-aware correction, and an
> applicability diagnostic.**
> Dana Paquin (Stanford University), Riddhiman Jain (Independent Researcher). Submitted to *Mathematical Biosciences and Engineering* (AIMS Press).
>
> Published deposit: DOI [10.5281/zenodo.22830134](https://doi.org/10.5281/zenodo.22830134)

LAVA reports a local genetic correlation only where both traits clear a univariate
heritability filter at *P* < 0.05. That filter acts on the denominator of a ratio
estimator, and it biases every estimate that survives it. This repository derives the
bias in closed form, supplies a correction and a diagnostic that says when the
correction may be trusted, and validates both against LAVA's own code path on six
psychiatric trait pairs.

The headline result: the closed form predicts each pair's measured attenuation slope
to within **0.007**, across a range from 0.579 to 0.923, with no parameter fitted to
the panel.

---

## The tool

If you only want to correct your own LAVA output, you need one file:

```r
source("code/R/lavagate.R")

d <- lava_gate_diagnostic(K, p1, p2, r_e, p1_all = ..., p2_all = ...)
if (d$usable) out <- lava_gate_correct(rho_hat, K, p1, p2, r_e)
```

`lavagate.R` is dependency-free base R. It consumes only what LAVA already returns —
`process.locus()$K`, `run.univ()$p`, `run.bivar()$rho` and its interval — plus the
sample-overlap correlation from cross-trait LD score regression. No re-run, no LD
reference, no individual-level data, so it applies to **already-published** results.

Run the diagnostic first. Where it fails, no corrected magnitude is defensible and
dividing by a genome-wide attenuation factor makes things worse, not better.

---

## Verifying before you trust any of it

All four run in a couple of minutes and need nothing but this repository:

```bash
python code/python/16_audit_numbers.py
```
```bash
python code/python/14_pair_pipeline.py --selftest
```
```bash
Rscript code/R/test_lavagate.R
```
```bash
Rscript code/replication/V0_exactness_check.R
```

- **16_audit_numbers.py** checks all 300 numeric claims in the manuscript against the
  result file each came from, *in both directions* — that the file gives the value, and
  that the manuscript still contains it. The second direction is what catches a re-run
  analysis whose prose was not updated. Exits non-zero on any failure.
- **14_pair_pipeline.py --selftest** checks the generic per-pair pipeline reproduces the
  worked pair to floating point, including its diagnostic verdict.
- **test_lavagate.R** checks the R implementation against the Python used for the paper's
  numbers, over 500 parameter combinations, and re-derives the worked pair's published
  slope, intercept and interval widths.
- **V0_exactness_check.R** proves the fast direct formulas used by the attribution
  control *are* LAVA's, by pushing identical effect matrices through both paths. Exits
  non-zero if agreement is not exact. It is a precondition for V1.

If any of them fails, stop.

---

## Layout

```
.
├── code/
│   ├── R/
│   │   ├── lavagate.R              the tool: correction + diagnostic
│   │   ├── test_lavagate.R         R-vs-Python cross-validation
│   │   ├── 00b, 01b                LDSC overlap and harmonisation for the panel
│   │   └── 00–04                   the worked pair's original scripts
│   ├── python/
│   │   ├── 01–13                   derivation, Tier A simulation, correction, figures
│   │   ├── 14                      generic per-pair pipeline (+ --selftest)
│   │   ├── 15                      panel table builder
│   │   ├── 16                      numeric audit (300 claims, both directions)
│   │   ├── 17                      Figure 6
│   │   └── V2–V5                   verification: nulls, calibration, theory vs measurement
│   └── replication/
│       ├── R1_multipair.R          the six-pair LAVA panel run
│       ├── R2_tierB_extended*.R    Tier B validation, worked pair and panel
│       ├── V0_exactness_check.R    direct formulas vs LAVA's own
│       ├── V1_gate_attribution.R   the 2x2 attribution control
│       └── run_panel.sh, resume_panel.sh, overnight.sh
├── figures/                        Figures 1–6 + S1, PNG and vector PDF
├── results/                        every file underlying a number in the paper
├── data/                           small inputs, and download locations for the bulk data
├── docs/                           LAVA source audit and the frozen pre-registrations
└── CORRECTIONS_2026-08-26.md       the verification audit: what was wrong, what changed
```

**Read `CORRECTIONS_2026-08-26.md`.** The paper was re-audited from source to method on
the assumption its findings were false until demonstrated otherwise. The central claim
survived; six numbers did not, including one fabricated sentence and one abstract claim
that quoted a prediction the paper itself reports as falsified. That file records all of
it, and lists which files elsewhere still carry superseded numbers.

---

## Running the analyses

Everything above runs from this repository alone. Reproducing the panel end-to-end also
needs the bulk GWAS and LD reference (~16 GB), which are not redistributed here —
`data/DATA_SOURCES_MULTIPAIR.md` has the figshare DOIs and download notes.

Scripts that need those files read the project root from an environment variable:

```bash
export LAVAGATE_BASE=/path/to/project        # the folder that contains SUBMISSION/
export LAVAGATE_LDSC_REF=/path/to/ldsc_ref   # only for the LDSC scripts
```

Order:

1. `code/R/01b_harmonize_panel.R` — converts the releases to LAVA input.
2. `code/R/00b_sample_overlap_ldsc_panel.R` — cross-trait LDSC for all eight phenotypes.
3. `code/replication/R1_multipair.R <GROUP>` — one genome pass per phenotype group,
   roughly 12–15 h each. `run_panel.sh` and `resume_panel.sh` orchestrate this.
4. `code/replication/R2_tierB_extended.R` and `R2_tierB_extended_multipair.R <GROUP>` —
   Tier B validation.
5. `code/python/14_pair_pipeline.py` per pair, then `15_build_table8.py`.
6. `V0` → `V1` → `V2`–`V5`, then `17_figure6_validation.py`.
7. `code/python/16_audit_numbers.py` last.

Three legacy scripts (`code/R/01_harmonize_sumstats.R`, `03_genomewide_bivar.R`,
`04_sample_overlap.R`) were originally run against an earlier project root. They now read
`LAVAGATE_BASE` like the rest; point it at wherever their inputs live.

`overnight.sh` rebuilds the panel inputs, verifies them, and launches all four passes
unattended. It will not launch the panel unless the rebuilt inputs pass their checks.

---

## Environment

- R 4.6.1 with **LAVA v0.1.5** (not redistributed — install from
  https://github.com/josefin-werme/LAVA), `data.table`, `jsonlite`
- Python 3.13 with `numpy`, `scipy`, `pandas`, `matplotlib`, `statsmodels`
- LDSC (Python 2.7) for the cross-trait sample-overlap step only

The verification and audit scripts need only Python and base R plus LAVA for `V0`.

---

## Licence

MIT, for the code and derived result files in this repository only. LAVA and the GWAS
summary statistics are not covered — see `LICENSE`.
