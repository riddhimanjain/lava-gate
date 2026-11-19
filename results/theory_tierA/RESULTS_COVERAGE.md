# Deliverable D — Coverage, width, and sign determinacy

Built from `out/coverage.csv`, `out/coverage_summary.json` and
`out/sign_determinacy.json`, all produced by `06_coverage.py`.
Replicates LAVA's `ci.bivariate()` parametric bootstrap exactly
(`confidence_intervals.R`), including the correlation clipping at ±0.99999 and
the final clip of the interval to [−1, 1].

> **Note on provenance.** Unlike `RESULTS_TIERA.md`, this file has **no generator
> script** — `06_coverage.py` writes the CSV and JSON but not this write-up, so this
> prose was assembled by hand and cannot be regenerated. Every figure in §§1–3 was
> re-verified against the JSON on 17 Aug 2026 and all match. If you change any number
> here, change it from the JSON.

---

## 1. Coverage is fine. That is not the story.

| | coverage | mean width | fraction with width > 1 |
|---|---|---|---|
| pre-gate | 0.974 | 1.651 | 0.936 |
| post-gate | 0.949 | 1.386 | 0.869 |

The pre-registered expectation was that gating would break calibration. **It does
not.** Post-gate coverage is 94.9% against a nominal 95%. Pre-gate intervals are
mildly conservative at 97.4%.

Reporting this plainly matters: the honest finding is that LAVA's intervals are
*correctly calibrated*, and any claim that the gate breaks their coverage would be
wrong.

## 2. The intervals are nearly vacuous

rho is bounded in [−1, 1], so the widest possible interval has width 2. The mean
post-gate width is **1.386** — roughly 69% of the entire parameter space. Nearly
nine in ten reported intervals are wider than 1.

An interval achieves nominal coverage trivially if it is wide enough. Coverage
and width have to be read together, which is why they are never separated here.

## 3. The consequence: the sign is undetermined

The diametric model predicts rho < 0; the convergent model predicts rho > 0. So
the question the paper asks is, at minimum, a question about **sign**. Fraction of
reported intervals that span zero — i.e. that cannot distinguish the two models:

| true rho | spans zero (r_e = 0) | sign determinate | spans zero (r_e = 0.10) |
|---|---|---|---|
| −0.4 | **0.861** | 0.139 | 0.927 |
| −0.2 | 0.915 | 0.085 | 0.937 |
| 0.0 | 0.937 | 0.063 | 0.925 |
| +0.2 | 0.920 | 0.080 | 0.913 |
| +0.4 | 0.877 | 0.123 | 0.838 |
| +0.6 | 0.779 | 0.221 | 0.701 |

**Even when the diametric model is true and strongly so (rho = −0.4), 86% of loci
that clear the gate return an interval that cannot tell you the sign.** Only about
1 locus in 7 is informative about the direction at all.

Sample overlap makes it worse, not better: at r_e = 0.10 and true rho = −0.4, 93%
of intervals span zero. Overlap both shifts the point estimate upward (see
`RESULTS_TIERA.md` §P4) and widens the effective uncertainty about sign.

## 4. What this does and does not license

> **⚠️ Revised 17 Aug 2026.** This section previously concluded that the applied result
> is an uninformative null, and called that "the paper's central claim." **The data
> contradicted it** and the paper no longer says this. Sign concordance came out at
> **57 of 78 loci positive (73.1%, P = 5.6 × 10⁻⁵)**, not near 50%, with all 10
> sign-determinate loci positive. The conditional argument below is still valid; the
> conclusion drawn from it was not. See `../out_R/applied_summary.txt`.

**The valid part.** If a real diametric signal of rho = −0.4 yields sign-determinate
intervals at only ~14% of loci, then a sign-concordance statistic near 50% *would
have been* uninformative — consistent with the convergent model, the diametric model,
and no local correlation at all. That was the anticipated outcome, and this section
was written to explain why it would not discriminate.

**What actually happened.** Concordance was not near 50%. It was 73.1% positive at
P = 5.6 × 10⁻⁵, and the direction of the determinate loci was unanimous (10 positive,
0 negative). So the applied result is **not** a non-detection within the window the
method can see; it favours the convergent model there.

**Two things this section's logic still does establish**, and both are load-bearing:

1. **The count of sign-determinate loci cannot discriminate.** Read the table above:
   determinacy is 0.139 at rho = −0.4 and 0.123 at +0.4 — nearly symmetric. A
   pre-registered test built on that count was retired for exactly this reason.
   **Direction** discriminates; the count does not.
2. **Non-detection elsewhere is still the correct reading for the other 96.9% of the
   genome.** Wide intervals mean a diametric relationship outside the analysable
   window would be invisible. "Unsupported where testable, untested elsewhere" is the
   claim — not "we can't tell anything," and not "convergence is established."

## 5. Caveats

- Sampling-model simulation only. Real per-locus K, Sigma and upstream attrition
  are Tier B's job; if the real K distribution differs materially from the values
  swept here, these fractions move.
- `n_boot` reduced from LAVA's 10,000 to 4,000 for the sign-determinacy table
  (compute budget). Coverage table uses the full 10,000.
- Intervals here are evaluated against the *true* rho. LAVA users compare them
  against zero, which is a different and easier question.
