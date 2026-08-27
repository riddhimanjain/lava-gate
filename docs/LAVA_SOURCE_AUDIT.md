# Further undocumented behaviours in LAVA v0.1.5

A systematic read of `analysis_functions.R`, `input_processing.R` and
`confidence_intervals.R` beyond the four attrition channels already in the paper.
Each item cites file and line, states the trigger condition, and says whether it
affects the ASD × SCZ analysis.

Ordered by how much they matter.

---

## 1. The reference-panel correction is dead code — `nref.scale` is always 1

`input_processing.R:38`

```r
process.locus = function(locus, input, phenos=NULL, min.K=2, ... , cap.estimates=T) {nref.correction=F
```

`nref.correction=F` is assigned **inside the function body, on the same line as the
opening brace**, where it reads as part of the signature. It is not a parameter and
cannot be set by the caller. Consequently, at line 134:

```r
loc$nref.scale = ifelse(nref.correction, (input$reference$sample.size - loc$K - 1)/(input$reference$sample.size - 1), 1)
```

**always evaluates to 1.** The degrees-of-freedom correction for finite LD reference
panel size — written, named, and applied at three separate points (`h2.obs` at :143,
`omega` at :166, and the univariate statistic at `analysis_functions.R:63`) — is
silently inert in this release.

**Trigger:** always.
**Magnitude:** at *K* = 270 and *N*<sub>ref</sub> = 99,339 the intended correction is
(99339 − 271)/99338 = 0.99727, i.e. 0.27%.
**Effect on our results:** none material. Tier B ran through LAVA's own code path, so it
already reflects `nref.scale = 1`. Tier A used the nominal correction, a 0.27% shift in
the gate threshold, well inside the 1.1% mean error of the closed form. **But the paper's
Methods must not present *c*<sub>ref</sub> as active** — corrected in the manuscript.

**Why it matters beyond us:** anyone reading the source to reconstruct the estimator —
which is exactly what a methods paper or a reimplementation does — will assume the
correction is applied, because the code says so at every point of use. The single line
that disables it is in the one place a reader will not look.

---

## 2. A fifth loss mechanism: the K/N ratio rule nulls heritability but not the correlation

`input_processing.R:161–163`

```r
thresh.ratio = 0.1
if (any(loc$K / loc$N > thresh.ratio)) {for (var in grep("^h2", names(loc), value=T)) loc[[var]][loc$K/loc$N > thresh.ratio] = NA}
```

When *K*/*N* exceeds 0.1 for a phenotype, its `h2.obs` and `h2.latent` are set to `NA`.
But:

- `loc$omega` is computed **after** this, at line 166, from `delta` and `sigma`. It is
  **not** nulled.
- `neg.var` (line 170) tests `diag(loc$omega)`, not `h2.obs`, so the phenotype is **not**
  flagged as failed.
- `univariate.test()` (`analysis_functions.R:63`) recomputes the statistic from `delta`
  and `sigma` directly, so the **P-value is still produced**.

The result is a phenotype that reports `h2.obs = NA` while still yielding a univariate
*P*-value and entering the bivariate analysis. Two notions of "estimable" — non-NA
heritability, and non-NA *P* — that disagree at exactly these loci.

**Trigger:** *K*/*N* > 0.1 for any phenotype in the locus.
**Effect on ASD × SCZ:** none. Maximum *K* is 607 against ASD's *N* = 46,350, giving
*K*/*N* = 0.013. The rule never fires here.
**Who it hits:** small GWAS. At *N* = 5,000 it fires for every block with *K* > 500;
at *N* = 2,000, for every block with *K* > 200 — i.e. most of the genome.

**Practical consequence:** an attrition audit keyed on `h2.obs` and one keyed on the
univariate *P* will return different counts, and neither is wrong. Ours is keyed on *P*,
which is why it is unaffected; a study keyed on `h2.obs` would report additional loss
with no mechanism named.

---

## 3. Negative variance drops the *strongest* loci as well as the weakest

`input_processing.R:142–143, 170`

```r
dtd = sum(loc$delta[,i]^2)
loc$sigma[i] = (1 - dtd) / (loc$N[i] - loc$K - 1)
...
neg.var = diag(loc$sigma) < 0 | diag(loc$omega) < 0
```

The residual variance is estimated as (1 − **δ**′**δ**)/(*N* − *K* − 1). When
**δ**′**δ** > 1 — which happens when a locus carries *strong* signal, not weak — `sigma`
goes **negative**, `neg.var` fires, and the phenotype is silently dropped by the same
Channel 2 mechanism the paper already describes.

So the mechanism we characterised as removing the underpowered trait has a second branch
that removes the **overpowered** one. Both set the same flag, and the warning message
(`input_processing.R:176`) says only "Negative variance estimate for phenotype(s) ...",
without distinguishing which side of the estimator failed.

There is a further consequence when sample overlap is supplied. Line 155:

```r
loc$sigma = sqrt(loc$sigma) %*% as.matrix(input$sample.overlap[...]) %*% sqrt(loc$sigma)
```

`sqrt()` of a negative diagonal yields `NaN`, so the whole `sigma` matrix — including the
entries for the *other*, healthy phenotype — becomes `NaN` before `neg.var` is evaluated
at line 170. One trait's strong locus can therefore contaminate its partner's variance.

**Trigger:** **δ**′**δ** > 1 at a locus.
**Effect on ASD × SCZ:** not observed — our attrition audit recorded only whether the
*P*-value was NA, which cannot separate the two branches. A one-line check is given in
§7 below and should be run before the next submission draft.
**Who it hits:** loci explaining a large share of trait variance — MHC in immune traits,
*APOE* in Alzheimer's, strong-effect eQTL blocks. Precisely the loci an applied study
most wants to report.

---

## 4. Results depend on which *other* phenotypes are in the input object

`input_processing.R:125–128`

```r
max.K = floor(max.prop.K * min(loc$N))
if (max.K < min.K) max.K = min.K
if (ncol(R) > max.K) R = R[,1:max.K]
```

`min(loc$N)` is the minimum sample size across **all phenotypes in the input object**,
not across the pair being analysed. Adding a third, smaller-*N* phenotype to
`process.input()` therefore reduces *K* at every locus for **every** pair — including
pairs that do not involve it.

Since the attenuation slope depends strongly on *K* (0.4555 at *K* = 100 against 0.6897
at *K* = 500, paper Table 3), a two-trait local genetic correlation is **not reproducible
from the trait pair alone**. It depends on the composition of the input object.

**Trigger:** any input object containing a phenotype with smaller *N* than the analysed pair.
**Effect on ASD × SCZ:** none — the input object holds exactly these two traits.
**Reporting implication:** studies must state the full phenotype set used to build the
input object, not just the pair. To our knowledge none do.

---

## 5. `param.lim` voids an entire multiple regression, not the offending coefficient

`analysis_functions.R:395`

```r
if (!multreg) { data[out.of.bounds, params] = NA } else { data[,params] = NA }
```

In the bivariate branch only the offending row is nulled. In the **multireg** branch
`data[,params] = NA` nulls **every row** — all predictors' gammas, confidence intervals
and *P*-values — when any single coefficient exceeds `param.lim`. The warning message
(line 392) names only the offending pair, so the user is not told that the entire model
was discarded.

Arguably defensible (a compromised model is compromised throughout), but it is not what
the message or the documentation says.

**Effect on ASD × SCZ:** none; no multiple regression was run.

---

## 6. Smaller items

**6a. `h2.obs` is censored at zero while `omega` is not.**
`input_processing.R:158–159` caps `h2.obs` and `h2.latent` at 0 under `cap.estimates=T`
(the default), but `loc$omega` is computed afterwards at line 166 and keeps its negative
diagonal. `run.univ()` (`analysis_functions.R:97`) caps again on output. The reported
local heritability is therefore a **censored** variable with a point mass at exactly 0,
while the *P*-value beside it is computed from the uncensored statistic.

This is why the correction in our paper inverts from the reported **P-value** rather than
from `h2.obs`: the *P*-value is uncensored and the mapping is invertible. Any analysis
that regresses `h2.obs` on annotations, or treats it as continuous, is working with a
censored outcome and should say so.

**6b. Vector condition in `ci.bivariate`.**
`confidence_intervals.R:31`

```r
if (sign(out$rho.lower)!=sign(out$rho.upper)) out$r2.lower = 0
```

`out` has one row per phenotype pair. For *P* > 2 this condition has length > 1, which is
an **error** in R ≥ 4.2 and a silent first-element-only bug before that. It is not reached
from `run.bivar()`, which passes 2 × 2 submatrices (line 165), but it is reached by any
direct call to `ci.bivariate()` with more than two phenotypes. `sign(0)` is also a
distinct third value, so an interval bound of exactly 0 takes the wrong branch.

**6c. `if (any(out.of.bounds))` is NA-unsafe.**
`analysis_functions.R:386`. If the estimate is `NaN`, `any()` returns `NA` and the `if`
aborts the locus with "missing value where TRUE/FALSE needed" rather than returning NA.

**6d. `cap()` uses `1:length(values)`.**
`analysis_functions.R:401`. On a zero-length input this becomes `c(1, 0)` and indexes out
of bounds. Not reachable from the current call sites, but the idiom is unsafe.

**6e. The gate has a different shape for continuous traits.**
`analysis_functions.R:64` uses `pchisq` for binary traits and `pf(stat/K, K, N-K-1)` for
continuous ones. The selection geometry, and therefore the attenuation, differs between
the two. **Our closed form is derived for the χ² case and applies to binary traits.** The
continuous case needs the F-distribution analogue of Δ<sub>i</sub>. This is now stated as
a scope limit in the manuscript. Note also that `N - K - 1 ≤ 0` yields `NaN` silently.

**6f. The interval is computed at a different parameter value than the point estimate.**
`confidence_intervals.R:8` clips the correlation to ±0.99999 before drawing, while
`cap.estimates` reports the point estimate as exactly ±1. For every locus whose estimate
is capped — all ten sign-determinate loci in our ASD × SCZ analysis — the reported ρ̂ and
its interval come from different Θ.

---

## 7. Diagnostics worth running

Two of the above are checkable on data already on disk, and both should be resolved
before submission.

**Which branch of Channel 2 fired.** Re-run the attrition audit recording `sum(delta^2)`
per phenotype per locus, then:

```r
# strong-signal branch: sigma < 0 because delta'delta > 1
strong <- audit$dtd_ASD > 1 | audit$dtd_SCZ > 1
# weak-signal branch: omega diagonal < 0
weak   <- !strong & is.na(audit$p_ASD)
table(strong, weak)
```

If `strong` is non-empty, the paper's Channel 2 description needs a second clause, and
the finding is stronger: the mechanism removes both tails.

**Whether the K/N rule ever fired.** `max(audit$K) / min(N_ASD, N_SCZ)` — if below 0.1 it
never did, which we expect here but have not recorded.

---

## 8. What changed in the manuscript as a result

- **Methods §2.1** no longer states that *c*<sub>ref</sub> is applied. It now gives the
  general form, notes that `nref.correction` is hardcoded false in v0.1.5 so the factor
  is 1, and states the 0.27% magnitude and that no reported number depends on it.
- **Methods §2.1** now states the derivation is for binary traits (χ² gate) and flags the
  continuous F-gate as out of scope.
- **Limitations** gains the input-object-composition reproducibility hazard (§4 above),
  since it bears directly on the reporting guidance the paper gives.

Items 2, 3, 4 and 5 are candidates for a short follow-up note. They are real, they are
undocumented, and none of them is needed to make the present paper's argument.
