# Tier A Results — gate-induced bias in LAVA local genetic correlation

**Generated:** 2026-08-17 14:20 UTC from `out/*.csv`. Predictions were frozen in `out/PREDICTIONS_PREREGISTERED.md` before any simulation ran.

**Scope.** Pure sampling-model simulation: validates the derivation against its own assumptions. It does **not** validate the assumptions. That is Tier B, which needs the real LD panel and the real LAVA code path.

> ### ⚠️ Superseded as the headline — read this first
> Every slope below is **Tier A's** answer, and Tier A has since been superseded on the headline number. Tier B (real loci, real per-locus *K*, real Σ, LAVA's own code path) gives an attenuation slope of **0.5512** (95% CI 0.514–0.589), and Tier A's **0.6064 falls outside that interval**. Per the rule committed to before either ran, Tier B wins.
>
> **The paper's headline attenuation is 0.551, not 0.606.** The numbers here remain valid as Tier A's own output and are what the Tier A vs Tier B comparison is built from — but do not quote them as the result. Tier B: `out_R/tierB.txt`.
>
> One further caveat on the `param.lim` section below: on real data `param.lim` discarded **zero** loci. Its effect here is a statement about the simulated configurations, not about ASD × SCZ.


## Scorecard against the pre-registered predictions

| # | Claim | Verdict |
|---|---|---|
| P1 | Attenuation toward zero | **Confirmed** |
| P2 | Attenuation = 1/sqrt(a1 a2), denominator only | **Falsified** |
| P3 | Near-multiplicative in rho | **Confirmed (5/6)** |
| P4 | Sample overlap creates an intercept | **Direction confirmed, magnitude badly wrong** |
| P5 | Attenuation worsens as K grows | **Falsified — and the prose contradicted my own formula** |

## P1 — Attenuation is real and large

| Scenario | K | fitted slope | interpretation |
|---|---|---|---|
| S1_asdlike_x_sczlike | 300 | 0.6064 | true rho=0.5 reports as 0.303 |
| S2_matched_moderate | 300 | 0.5721 | true rho=0.5 reports as 0.286 |
| S3_overlap_re05 | 300 | 0.6062 | true rho=0.5 reports as 0.303 |
| S4_overlap_re10 | 300 | 0.6010 | true rho=0.5 reports as 0.300 |
| S5_smallK | 100 | 0.4507 | true rho=0.5 reports as 0.225 |
| S6_largeK | 500 | 0.6880 | true rho=0.5 reports as 0.344 |

At the ASD-like x SCZ-like configuration the slope is **0.6064** — roughly 39% of the true signal is lost before it is ever reported.


## P2 — Falsified: the denominator account is incomplete

| Scenario | sim slope | 1/sqrt(a1 a2) | rel. error |
|---|---|---|---|
| S1_asdlike_x_sczlike | 0.6064 | 0.5392 | 12.5% |
| S2_matched_moderate | 0.5721 | 0.4944 | 15.7% |
| S3_overlap_re05 | 0.6062 | 0.5392 | 12.4% |
| S4_overlap_re10 | 0.6010 | 0.5392 | 11.5% |
| S5_smallK | 0.4507 | 0.3390 | 33.0% |
| S6_largeK | 0.6880 | 0.6351 | 8.3% |

Errors of 8–33% mean the gate inflates the **numerator** too, partially offsetting the denominator inflation. The true attenuation is *less* severe than the leading-order account predicts. `04_refined_theory.py` derives the second-order form from the exact noncentral Wishart covariance.


## Refined theory, validated on held-out parameters

Six parameter settings never used in building the refinement, 7 values of rho each, 300,000 simulated loci per cell.

- leading order (v1): mean abs error **0.0492**, max 0.1709
- second order (v2): mean abs error **0.0066**, max 0.0348
- 88% of held-out cells within 0.01

v2 remains statistically distinguishable from the simulation (max |z| = 47 at these Monte Carlo sizes). It is a second-order expansion, not an identity, and should be reported as such.


## P4 — Sample overlap manufactures correlation from nothing

| r_e | true rho | simulated E[rho_hat \| G] | v1 predicted | v2 predicted |
|---|---|---|---|---|
| 0.2 | 0 | 0.1562 | ~0.0002 | 0.1597 |
| 0.08 | 0 | 0.0346 | ~0.0002 | 0.0373 |
| 0.05 | 0 | 0.0340 | ~0.0002 | — |
| 0.1 | 0 | 0.0666 | ~0.0002 | — |

**This is the most consequential single result.** Two traits with *zero* true local genetic correlation, measured in overlapping samples, report a positive local correlation after the gate. At r_e = 0.20 the spurious value is ~0.156. The leading-order formula under-predicted this by a factor of ~280; the second-order form captures it.


## P5 — Falsified, and the error was mine

The pre-registration prose asserted attenuation *worsens* as K grows. The closed form it was supposedly describing predicted the opposite, and the simulation agrees with the formula:

| K | sim slope |
|---|---|
| 100 | 0.4507 |
| 300 | 0.6064 |
| 500 | 0.6880 |

At fixed theta, larger K means larger noncentrality (lambda = K*theta), so the gate selects less aggressively. Recorded as a falsified prediction rather than silently corrected, because the prose and the mathematics disagreed and only the simulation revealed it.


## Channel 3: param.lim, tested separately

**Hypothesis tested and rejected:** that v2's residual error was caused by `param.lim` truncation. Re-running with the truncation disabled barely moved the error (mean |err| 0.00678 without vs 0.00670 with). The residual is delta-method curvature, not channel 3.

`param.lim` is nonetheless a real and separate distortion:

| \|rho\| | fraction dropped | shift in reported mean |
|---|---|---|
| 0.00 | 0.0008 | -0.0001 |
| 0.20 | 0.0015 | -0.0008 |
| 0.30 | 0.0018 | +0.0013 |
| 0.45 | 0.0041 | -0.0026 |
| 0.60 | 0.0070 | +0.0044 |
| 0.70 | 0.0131 | -0.0068 |
| 0.90 | 0.0300 | -0.0116 |

Note the reporting rule is discontinuous: rho_hat in (1, 1.25] is **capped to 1**, while rho_hat > 1.25 becomes **NA**. Two nearly identical estimates receive opposite treatment.


## What Tier B must test

Everything above assumes the sampling model is correct. Tier B checks that assumption against reality:

1. Real per-locus K from the LD panel — not the assumed values swept here.
2. Real per-locus Sigma and the actual ASD/SCZ sample overlap.
3. The real LAVA code path, including `process.locus` upstream attrition (channel 1), which this simulation does not model at all.
4. Whether the empirical distribution of theta across the 2,495 blocks puts ASD x SCZ in the severe-attenuation regime.

