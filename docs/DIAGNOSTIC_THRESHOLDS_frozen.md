# Applicability-diagnostic thresholds — frozen before the panel was run

**Written:** 23 August 2026, after `R2_tierB_extended.R` was launched and after the
panel GWAS were harmonised, but **before any of `R1_multipair.R`'s P1, P2a, P2b, P3 or
P5 output existed**. Nothing in this file may be edited to match an outcome. It shares
the status of `PREDICTIONS_PREREGISTERED_frozen.md` and `RESEARCH_FRAMING_frozen.md`:
its entire value is that it was not revised afterwards.

---

## Why this file exists

Section 3.5 of the manuscript specifies the applicability diagnostic **qualitatively** —
a stability check and a floor-distance check, "both must pass before a corrected
magnitude is reported" — but states no numeric threshold. That was adequate while one
pair had been run and the answer was a clear failure on both checks. It is not adequate
for scoring five further pairs: without a threshold fixed in advance, the "Usable?" row
of Table 8 would be assigned after seeing the numbers, and the diagnostic would be
unfalsifiable decoration rather than a decision rule.

## The calibration problem, stated honestly

These thresholds were chosen with knowledge of the completed pair (P4, ASD × SCZ). That
is unavoidable — P4 is the only pair whose diagnostic outcome was known, so it is the
only available calibration point, and the manuscript already declares it a failure. What
is *not* contaminated is the five pairs that had not been run when this file was written.

The honest description, which the manuscript must use, is therefore: **the thresholds
were set on the one completed pair and frozen before the remaining five were computed.**
They are not derived from first principles and are not claimed to be. A reader who
disagrees with a cut-point can read the underlying continuous quantities, which are
reported per pair in `pair_summary.json` and never collapsed only to a verdict.

---

## The three criteria

A pair is scored **usable** — meaning a corrected magnitude may be reported for it — only
if **all three** pass. Any failure means the uncorrected estimate is reported with the
attenuation stated as a bound, as was done for P4.

### C1 — Stability of the aggregate slope

Recompute the aggregate correction slope under *K* stratifications into 1, 2, 3, 4, 6 and
8 strata (`11_aggregate_correction.py`'s existing grid).

> **Pass if** max(slope) / min(slope) **< 1.25**

*Interpretation.* The corrected magnitude is a division by this slope, so a fold change of
1.25 means the reported magnitude is pinned only to about ±25%. Beyond that the choice of
an arbitrary analysis nuisance — how finely to bin *K* — moves the headline number more
than the correction itself is worth.

*P4 value:* 2.66. **Fails.**

### C2 — Distance from the selection floor

For each trait, compare observed mean θ̂ at the reporting loci against the conditional
distribution of θ̂ under θ = 0 (what a null locus reports having passed the gate), within
*K* bins, then pool across bins weighting by the number of loci.

> **Pass if, for BOTH traits,** pooled *z* ≥ **5** **and** the loci-weighted fraction
> falling below the null conditional mean ≤ **0.40**

*Interpretation.* Stage 1 inverts *g*~K~, which flattens as θ → 0; near the floor the
inversion is not identified. Pooling by loci rather than taking the minimum across bins is
deliberate: a single sparse *K* bin should not veto a pair, and a single sparse bin should
not rescue one either. Both traits must clear it, because the ratio's denominator contains
both.

*P4 values:* weaker trait (ASD) *z* spans −1.67 to +3.05 with 50–100% of loci below the
null mean in most bins; stronger trait (SCZ) spans +1.95 to +48.3. **Fails on ASD.**

### C3 — The gate output is not mostly null

From the univariate pass rates alone: under a complete local null the gate passes at
exactly α = 0.05, so the fraction of observed passes attributable to chance is
α / (observed pass rate).

> **Pass if, for BOTH traits,** the fraction of gate passes expected null **≤ 0.40**
> (equivalently, enrichment over α ≥ 2.5×)

*Interpretation.* This is independent of the correction machinery entirely — it uses only
*P*-values and counting. It is included because C1 and C2 both depend on the same
inversion, so a pair could conceivably pass them for the wrong reason. C3 cannot.

*P4 value:* ASD enrichment 2.05×, ≈49% of passes expected null. **Fails.**

---

## Pre-registered expectations for the five unrun pairs

Recorded so they can be scored against the outcome, not so they can be quietly dropped if
wrong. Two of the previous five pre-registered predictions were falsified and both stayed
in the paper; the same rule applies here.

| Pair | Expectation | Reasoning |
|---|---|---|
| **P1** SCZ × BIP | **Passes all three.** | Both traits well powered; this is the pair included specifically to demonstrate the correction succeeding. If P1 fails, the paper has no worked success case and must say so. |
| **P2a** MDD × BIP (full) | Passes C2/C3; **C1 uncertain.** | Both well powered, but the manufactured intercept is large here, and the aggregate slope's stability under a large intercept has not been tested. |
| **P2b** MDD × BIP (no UKBB) | Same as P2a. | Power is only slightly lower than P2a. |
| **P3** ASD × ADHD | **Fails**, on ASD, as in P4. | ASD is the weaker trait again and its Channel 2 loss rate is a property of the ASD GWAS, not of its partner. |
| **P5** AN × ASD | **Fails**, on both traits. | Both weakly powered; this pair is included as the floor case and is expected to yield very few testable loci. |

**The P2 prediction that matters most is not about usability at all.** It is that the
fitted intercept moves from ≈ 0.68 × *r*~e~(P2a) to ≈ 0.68 × *r*~e~(P2b) when UK Biobank
participants are removed from both studies, with *r*~e~ measured by cross-trait LDSC in
`00b_sample_overlap_ldsc_panel.R` and the disorders, LD reference and blocks held fixed.
That prediction is scored regardless of whether either arm is usable.

---

## What is deliberately NOT a criterion

- **Number of reported loci.** A pair with few reporting loci is imprecise, not
  unidentified. Precision is reported separately; it is not part of applicability.
- **Whether the corrected magnitude looks plausible.** Scoring a diagnostic on whether it
  produces an agreeable answer is the failure mode the diagnostic exists to prevent.
- **Agreement between the aggregate and per-locus corrections.** The per-locus variant is
  known to be worse (bias 0.106 against 0.026) and is not recommended; disagreement
  between them carries no information about applicability.
