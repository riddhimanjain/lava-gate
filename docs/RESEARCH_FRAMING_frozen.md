# Research Framing — ASD × SCZ, Gate-Induced Bias in Local Genetic Correlation

**Written:** 14 Aug 2026, before any result in the paper existed.

> ## ⚠️ FROZEN RECORD — deliberately not updated
>
> This is the framing, deliverable list and kill conditions as committed **before** the analysis
> ran. Its entire value is that it was not edited afterwards, so it can be read against what
> actually happened. **Do not "fix" the claims below to match the results.** Where it was
> superseded, that is recorded here and nowhere in the body:
>
> | Stated below | What the data gave |
> |---|---|
> | Three attrition channels (§3) | **Four.** A fourth was found during execution: `process.locus()` silently drops a *phenotype*, not the locus, on negative variance. It turned out to be the dominant channel (ASD lost at 59.5% of loci). |
> | Channel 3 (`param.lim`) "most likely to be novel to readers" (§3) | Real, undocumented, and it removed **zero** loci here. Demoted from a finding to a mechanism that will bind for other trait pairs. |
> | "The negative result is evidence of non-detection" (§1) | **Contradicted.** Within the analysable 3.1% the result is not null: 57 of 78 loci positive, *P* = 5.6 × 10⁻⁵, all 10 sign-determinate loci positive. The paper follows the data and favours the convergent model *where testable*. |
> | The bias is "non-multiplicative" (§2, §4-A) | Falsified as stated. It is near-multiplicative; the leading-order denominator-only account is nonetheless wrong by 8–33%, which is why the second-order form exists. |
> | Kill condition "if the bias turns out to be negligible" (§6) | Not triggered. Attenuation slope 0.551 — roughly 45% of signal destroyed. |
> | Kill condition "if prior art exists" (§6) | Not triggered. Novelty check completed against primary sources; see `TASKS_EXECUTION_LOG.md` Task 5. |
> | `[PRIOR]` markers below | All resolved. Every number in the paper was independently re-derived; nothing from the discarded round survives in it. |
>
> Results live in `analysis/out/`, `analysis/out_R/` and the paper itself. Execution detail is in
> `TASKS_EXECUTION_LOG.md`.

**Status:** framing only. No results in this file. Numbers marked `[PRIOR]` come from the discarded
round and were **not** usable — each had to be independently re-derived before appearing anywhere.

---

## 1. The claim, in one sentence

The diametric-vs-convergent question for autism and schizophrenia has been tested with an
estimator that, at current sample sizes, **cannot return the answer** — because local genetic
correlation is only reported at loci that survive a significance gate, and that gate both selects
which loci report and attenuates the values at the survivors.

The negative result is therefore not evidence for convergence. It is evidence of non-detection,
and those are different claims.

## 2. Why this is worth a paper

Two things have to both be true, and the second is the one that can sink us:

1. **The bias is real and material.** Testable. Sections 4–6 below.
2. **Nobody has characterised it for local r_g.** As far as located so far, yes — the LAVA
   documentation describes the univariate filter as a way to "filter out non-associated loci that
   may yield unstable correlation estimates," i.e. it is presented as a *stabiliser*, with no
   accompanying characterisation of the bias it induces in the surviving estimates. That is the
   gap. **This must be re-verified against the published LAVA supplement and the 2023 benchmarking
   paper before we commit** — it is the single highest-risk assumption in the project.

The nearest prior art is the GWAS winner's-curse literature. The distinction we must defend:
winner's curse selects on the estimand itself; here selection acts on the **denominator**, which
is correlated with but not identical to the numerator. That asymmetry is what makes the bias
non-obvious and, we predict, non-multiplicative.

## 3. The mechanism — three attrition channels, not one

Located in the LAVA source, and all three must be in the paper:

| # | Channel | Where | Effect |
|---|---|---|---|
| 1 | **Upstream coverage** | `process.locus()` | Locus dropped before any test — insufficient SNPs after harmonisation with the LD panel |
| 2 | **The univariate gate** | `analysis_functions.R:42–44`, `univ.thresh=.05` | Bivariate test runs only if *both* traits clear P<0.05 univariate |
| 3 | **Parameter truncation** | `filter.params()`, `analysis_functions.R:384–385`, `param.lim=1.25` | ρ̂ set to `NA` post-hoc when out of bounds — censoring the tail of the very distribution we are characterising |

Channel 3 is the one most likely to be novel to readers. It is not a filter on the data; it is a
filter on the answer.

## 4. Deliverables

**A — Analytic.** Closed-form expression for `E[ρ̂ | gate]` under LAVA's own sampling model.
Must express the attenuation in terms of: true ρ, per-locus SNP count K, effective N per trait,
local h² for each trait, and sample overlap. Prediction to be stated *before* simulation:
the bias is **not** a simple multiplicative attenuation, because numerator and denominator are
built from shared Z-scores.

**B — Tier A simulation.** Synthetic, generated from the assumed sampling model. Validates A
against its own assumptions. Cheap, fast, and proves only internal consistency.

**C — Tier B simulation.** Real loci, real per-locus K, real Σ, run through the actual LAVA code
path. This is the one that counts — it validates A against reality rather than against itself.
`[PRIOR]` found real K median ≈ 287.5 against an assumed 100, which if it replicates means Tier A
alone would have been misleading.

**D — Coverage.** LAVA's own `ci.bivariate` intervals, pre- and post-gate. Report both coverage
*and* interval width — an interval can achieve nominal coverage by being vacuous, and
`[PRIOR]` suggests that is exactly what happens here. Coverage without width is not a result.

**E — Applied.** Re-run ASD × SCZ. Report the detectability ceiling: how many of 2,495 blocks are
jointly estimable, and which trait sets the ceiling. Then the diametric/convergent test, with the
power to detect it stated explicitly rather than assumed.

## 5. Test plan, in order

1. Rebuild the input object; verify against `Section 7/results/` baselines as a regression check.
2. Locate and quantify all three attrition channels on the real ASD × SCZ data.
3. Derive A. State predictions in writing before running B or C.
4. Run B. Then C. Compare against the *pre-registered* predictions from step 3.
5. Run D — coverage and width, both gates.
6. Run E — the applied re-analysis and the power statement.
7. Adversarial pass: every number re-derived independently before it enters the manuscript.

Step 3-before-4 ordering is not bureaucratic. A derivation checked after seeing the simulation is
not a derivation, it is a curve fit.

## 6. Kill conditions

Stated now, while it is cheap to be honest about them:

- **If the bias turns out to be negligible** (attenuation within a few percent of 1), Deliverable A
  is still publishable as a negative methodological result, but the paper's framing collapses —
  the ASD×SCZ null would then stand as a real null. We report that. We do not rescue the framing.
- **If prior art exists** characterising this for local r_g, the contribution reduces to
  replication plus the ASD×SCZ application. Survivable, much smaller. Check this early.
- **If the derivation and Tier B disagree** and the disagreement cannot be traced to a stated
  assumption, we report the disagreement. `[PRIOR]` records the ρ-dependence prediction *failing*
  once already — the coupling term came out near-negligible. That failure is data, and if it
  reproduces it goes in the paper.

## 7. Constraints carried over

- `ASD_SCZ_Introduction_Methodology.pdf` is already with the professor. Title, Introduction and the
  nine Methods subsections are fixed. New work fits inside that structure or the change gets
  communicated explicitly.
- European-ancestry only, forced by the UKB LD reference and the Grove ASD GWAS. State it as a
  limitation; do not let a reviewer find it first.
- `ASSISTANCE_LOG.md` stays current from here. Venue AI-disclosure policy is decided at submission,
  but the record it is written from has to be accurate as it happens, not reconstructed afterwards.
