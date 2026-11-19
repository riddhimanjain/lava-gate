# Verification audit, 26 August 2026 — what was checked, what was wrong, what changed

Before submission the paper was re-audited from source to method on the working
assumption that its findings were **false until demonstrated otherwise**. The
specific worry driving it: that we had attributed to LAVA an effect that was
really an artefact of our own code, our own post-processing, or our own choice of
null.

The central claim survived. Six numbers did not. This file records all of it,
including the two errors that would have gone to press.

---

## 1. The central claim was tested directly, and holds

The paper says LAVA's univariate filter causes the attenuation of ρ̂. Every slope
in the original draft was measured on **gated replicates only**, which cannot
distinguish that claim from two rivals: that LAVA's `cap.estimates` / `param.lim`
truncation causes it, or that ρ̂ is simply a biased ratio estimator.

`code/replication/V1_gate_attribution.R` runs a 2 × 2 factorial on one replicate
stream — {estimable, gated} × {no truncation, LAVA's truncation} — for all six
pairs. `V0_exactness_check.R` first proves the fast direct formulas are LAVA's,
comparing identical **δ** matrices through both paths: **0.000e+00 difference on
all 10,000 univariate *P*-values and all 1,909 defined ρ̂ values**, after matching
the `signif(., 6)` rounding LAVA applies to its own output.

Result: with truncation disabled entirely the gated slope is still 0.647 (P4) and
0.601–0.957 across the panel. Truncation moves the slope by only 0.020–0.034.
**The attribution to the filter is correct.** Added to the paper as Tier C
(Section 2.6), Table 5 and Figure 6b.

Two bonuses. Arm D reproduced every published Tier B slope to within 0.003 through
an independent implementation, a different seed and twice the replicates — an
unplanned reproducibility check on the headline numbers. And the exactness check
turned up a fifth, undocumented loss mechanism: **`run.bivar()` aborts** rather
than returning `NA` on a locus with negative local variance
(`analysis_functions.R:375`, then `387`), so the filter is load-bearing for the
code to run at all. Now Table 1, row 5.

---

## 2. Errors found and corrected

| # | Where | Was | Is | How it was found |
|---|---|---|---|---|
| 1 | Abstract | "closed-form slope ranges from 0.33 to 0.51" | **0.455 to 0.690** | 0.33 and 0.51 are `attenuation_1_over_sqrt_a1a2` and the v1 `pred_slope` — values from prediction **P2, which the paper itself reports as falsified**. The abstract was quoting the withdrawn theory. |
| 2 | §3.3, Table 11 | "Predicted intercept = 0.68 × *r*~e~" for all six pairs | closed-form *I*~locus~ per pair | The rule over-predicts P1 by 0.117. Mean absolute error 0.0432 against 0.0050 for the closed form. 0.68 was read off Tier A at *K* = 300 in a low-power cell; *I*~locus~ ∝ Δ~i~ → 0 as power grows, so it cannot be constant. |
| 3 | §3.3 | "UK Biobank contributes 59,851 of 113,154 cases+controls to the full MDD release, and 41,917 of 371,549 to the full BIP release" | real per-release counts read from `multipair.info.txt`; differences 23,004 and 59,567 | **Fabricated.** 41,917/371,549 are simply BIP's cases and controls, not a UKBB contribution. `DATA_SOURCES_MULTIPAIR.md` explicitly says these cells are "read from file… Do not fill these in by hand." |
| 4 | §3.7, C3 | weaker trait "enrichment 2.05×, ≈49% of passes expected null, 31.8 joint passes" | **1.00×, ≈100% null, 65.2 joint passes** | A *conditional* pass rate (10.26%, among Channel-2 survivors) was compared against an *unconditional* null (0.05). Channel 2 keeps exactly the loci with *U* ≥ *K*, so the right null is 0.05/*P*(χ²~*K*~ ≥ *K*) ≈ 0.103. The paper warns against this exact conflation two paragraphs earlier. |
| 5 | §2.6 | "the prediction that the bias would be non-multiplicative" was falsified | **P5** was falsified: attenuation was predicted to *worsen* as *K* grows; it eases (0.4555 at *K* = 100 → 0.6897 at *K* = 500) | The frozen file lists five claims P1–P5. P3 predicted near-multiplicativity and was **confirmed**. The draft named the wrong prediction. |
| 6 | §4.1 item 3 | "slope varies from 0.35 to 0.81 across real loci" | **0.371 to 0.971** | Stale: the 18-locus run's range, superseded by the extended run in §3.2 but not updated in the checklist. |

Two further softenings, not errors but overclaims:

- **§2.3 "exactly linear"** now separates what is exact (linearity of the closed
  form in ρ — an algebraic identity, no ρ² term is generated) from what is
  approximate (the linearisation of *S*₁₂ on the diagonals, hence the *values* of
  *S*~locus~ and *I*~locus~). The empirical support is unchanged and strong
  (quadratic term *P* = 0.936, 0.68, 0.96 on the three pairs tested).
- **Abstract and Conclusion** no longer generalise P4's 38% attenuation to the
  method. The panel range is 0.579–0.923.

---

## 3. A new result, from testing the theory properly

The closed form had only ever been checked against Tier A — simulation drawn from
the model the derivation assumes, which as §2.6 says "cannot detect an error in
them". `V4_closedform_vs_tierB.py` checks it against Tier B on the same loci, per
pair. It predicts the measured slope of **all six pairs to within 0.007** (mean
0.0048), across a range from 0.579 to 0.923, with no parameter fitted to the
panel. This is now the paper's strongest evidence (Table 11, Figure 6a) and it was
absent from the draft.

---

## 4. One anomaly we could not fully explain, now a stated limitation

Under a complete local null a phenotype survives Channel 2 at
*P*(χ²~*K*~ ≥ *K*) ≈ 48.7%, and real signal can only raise that. ASD sits *below*
it in all three pairs it appears in: 40.5%, 42.0%, 40.0%. Parsing every LAVA
warning in the run logs excludes a second drop mechanism (1,407/1,407 for P3,
1,474/1,474 for P5 carry the documented negative-variance warning), so the χ²
reference is itself ~2% mis-calibrated for that trait. A one-parameter variance
rescaling does not identify against "the trait has signal" (every well-powered
trait fits the opposite side of 1). Reported in §4.4 with its bounded consequence:
the enrichment moves between 1.02 and 1.36 depending on the null, against 3.4–8.1
for every trait with real signal, so no verdict changes.

Separately, and cleanly: ASD's local heritability *P*-values are
**indistinguishable from a complete local null** (Kolmogorov–Smirnov *P* = 0.47,
0.078, 0.86 in P4/P3/P5) while all nine other trait-analyses reject at
*P* < 10⁻¹²⁰. The three pairs the diagnostic rejects are exactly the three
containing ASD.

---

## 5. Pre-registration, scored honestly

- Of five frozen claims (`docs/PREDICTIONS_PREREGISTERED_frozen.md`): **P2 and P5
  falsified**, P1/P3/P4 confirmed.
- Of five frozen per-pair usability expectations: four correct; **P5 partly wrong**
  (predicted to fail "on both traits"; fails on ASD only — AN clears C2 and C3).
- The separately frozen 0.68 × *r*~e~ intercept prediction: **falsified**.

Frozen files were **not** edited. Neither were `docs/EXECUTION_LOG.md` or
`docs/SESSION_LOG_2026-08-23.md`, which are historical records.

---

## 6. Everything else checked and found correct

LAVA source citations (`input_processing.R:38, 125–128, 134–166, 169–171`;
`analysis_functions.R:24, 42–44, 64, 172–176, 375, 387`), v0.1.5, `nref.scale`
= 1, the noncentral Wishart covariance and the algebra of *S*~locus~ / *I*~locus~
(re-derived independently), Table 11's every cell against each pair's raw
`attrition_audit.csv`, all Channel-2 losses on the processed denominator, the
diagnostic C1/C2/C3 values, Tables 3, 6, 8, 9, coverage and width, the Fisher
odds ratio (2.37 is the conditional MLE, which is what its exact interval
1.32–4.53 belongs to — the sample OR is 2.376), case/control counts, SNP
intersections, and Tier B's SEs, *R*², cell count and *K* range.

---

## 7. Verification is now machine-checked

`code/python/16_audit_numbers.py` grew from **94 to 269 claims**. Table 11 had
*zero* coverage before; it and the diagnostic verdicts, Tier C, the closed-form
comparison, the null baselines and the data provenance are all covered now. The
manuscript-side match was also tightened: it previously used a plain substring
test, which passes spuriously ("0.62" sits inside "0.6236"), and now requires the
match not be flanked by a digit. All 269 pass under the strict matcher.

`V0_exactness_check.R` exits non-zero if it fails and is a precondition for V1.

---

## 8. Files carrying superseded numbers — read the paper, not these

These were **not** rewritten. They are talk notes and study aids, not submission
material, and rewriting nine of them would risk introducing fresh errors:

`presentation/TALK_10min_LAVA_gate.md`, `STUDY_MODULES/00_START_HERE.md`,
`05_winners_curse_comparison.md`, `10_sample_overlap.md`,
`16_project_status_and_panel.md`, `18_cheat_sheet.md`,
`FIGURES_WALKTHROUGH.md`, `FINAL_10MIN_REVIEW.md`, `FULL_SIMPLE_WALKTHROUGH.md`.

Where any of them states 0.68 × *r*~e~, a 2.05× enrichment, 49% of passes null,
31.8 joint passes, a 0.35–0.81 per-locus range, or a 0.33–0.51 closed-form range,
**`manuscript/manuscript.md` supersedes it.** Section 2 above gives the
replacements.
