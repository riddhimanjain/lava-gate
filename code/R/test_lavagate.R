#!/usr/bin/env Rscript
# test_lavagate.R -- check the R implementation against the Python one.

HERE <- tryCatch(dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
         error = function(e) ".")
source(file.path(HERE, "lavagate.R"))

SUB   <- normalizePath(file.path(HERE, "..", ".."))
PYDIR <- file.path(SUB, "code", "python")
TOL   <- 1e-8

fail <- 0L
DIFFS <- list()
chk <- function(name, got, want, tol = TOL) {
  d <- max(abs(as.numeric(got) - as.numeric(want)), na.rm = TRUE)
  ok <- is.finite(d) && d <= tol
  if (!ok) fail <<- fail + 1L
  DIFFS[[name]] <<- d
  cat(sprintf("  [%s] %-46s max|diff| = %.3g\n", if (ok) "ok  " else "FAIL", name, d))
}

KS     <- c(50L, 100L, 250L, 300L, 500L)
THETAS <- c(0.005, 0.02, 0.065, 0.15, 0.30)
RES    <- c(0.0, 0.02756, 0.10, 0.20)

cat("=== R vs Python: Delta, g_K, slope, intercept ===\n")

tmp <- tempfile(fileext = ".csv")
py <- sprintf('
import os, itertools, csv
os.chdir(r"""%s""")
from importlib import machinery
R = machinery.SourceFileLoader("refined", "04_refined_theory.py").load_module()
C = machinery.SourceFileLoader("corr", "08_correction.py").load_module()
import numpy as np
rows = []
for K in %s:
    for t1 in %s:
        for t2 in %s:
            D1, _ = R.delta_i(K, t1); D2, _ = R.delta_i(K, t2)
            for re in %s:
                sl, ic = C.locus_slope_intercept(K, t1, t2, D1, D2, re)
                rows.append(dict(K=K, t1=t1, t2=t2, r_e=re,
                                 D1=float(D1), D2=float(D2),
                                 slope=float(sl), intercept=float(ic)))
with open(r"""%s""", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
',
  PYDIR,
  paste0("[", paste(KS, collapse = ","), "]"),
  paste0("[", paste(THETAS, collapse = ","), "]"),
  paste0("[", paste(THETAS, collapse = ","), "]"),
  paste0("[", paste(RES, collapse = ","), "]"),
  tmp)

pyfile <- tempfile(fileext = ".py"); writeLines(py, pyfile)
st <- system2("python", pyfile, stdout = TRUE, stderr = TRUE)
if (!file.exists(tmp)) { cat(paste(st, collapse = "\n"), "\n"); stop("python side failed") }
P <- read.csv(tmp)

P$D1_R <- mapply(function(K, t) lg_delta(K, t), P$K, P$t1)
P$D2_R <- mapply(function(K, t) lg_delta(K, t), P$K, P$t2)
si <- lava_slope_intercept(P$K, P$t1, P$t2, P$D1_R, P$D2_R, P$r_e)

chk("Delta_1", P$D1_R, P$D1)
chk("Delta_2", P$D2_R, P$D2)
chk("slope",     si$slope,     P$slope)
chk("intercept", si$intercept, P$intercept)
cat(sprintf("  (%d parameter combinations)\n", nrow(P)))

cat("\n=== Stage 1 inverts stage 1 ===\n")
rt <- do.call(rbind, lapply(KS, function(K) {
  th <- THETAS
  g <- th + vapply(th, function(t) lg_delta(K, t), numeric(1)) / K
  data.frame(K = K, theta = th, back = lava_invert_theta(g, K)$theta)
}))
chk("g_K inverted returns theta", rt$back, rt$theta, tol = 1e-5)

cat("\n=== The selection floor is strictly positive ===\n")
fl <- vapply(KS, function(K) lg_grid(K)$g[1], numeric(1))
cat(sprintf("  g_K(0) at K = %s : %s\n",
            paste(KS, collapse = ", "),
            paste(sprintf("%.5f", fl), collapse = ", ")))
if (any(fl <= 0)) { cat("  [FAIL] g_K(0) must be > 0\n"); fail <- fail + 1L
} else cat("  [ok  ] a null locus reports positive heritability once gated\n")

cat("\n=== End to end on the completed pair's 78 reported loci ===\n")
aud_f <- file.path(SUB, "results", "applied_tierB", "attrition_audit.csv")
ref_f <- file.path(SUB, "results", "correction", "applied_correction_summary.json")
if (file.exists(aud_f) && file.exists(ref_f)) {
  a <- read.csv(aud_f)
  rep_ <- a[a$gate_pass %in% TRUE & is.finite(a$rho_hat), ]
  cat(sprintf("  reported loci: %d (expected 78)\n", nrow(rep_)))
  if (nrow(rep_) != 78) { cat("  [FAIL] wrong number of reported loci\n"); fail <- fail + 1L }

  out <- lava_gate_correct(rep_$rho_hat, rep_$K, rep_$p_ASD, rep_$p_SCZ,
                           r_e = 0.02756,
                           rho_lower = rep_$rho_lower, rho_upper = rep_$rho_upper)
  ag <- attr(out, "aggregate")
  ref <- tryCatch(jsonlite::fromJSON(ref_f), error = function(e) NULL)
  if (!is.null(ref)) {
    chk("aggregate slope vs published",     ag$slope,     ref$aggregate_slope, 1e-7)
    chk("aggregate intercept vs published", ag$intercept, ref$aggregate_intercept, 1e-7)
    chk("mean rho_hat vs published",        ag$mean_rho_hat, ref$mean_rho_hat, 1e-9)
    chk("mean CI width vs published",
        mean(out$ci_width, na.rm = TRUE), ref$mean_ci_width, 1e-9)
    chk("n floored vs published",
        sum(out$slope <= 0), ref$n_floored_no_recoverable_signal, 0)
  }
  cat(sprintf("  mean width %.3f -> corrected %.3f (inflation %.2fx)\n",
              mean(out$ci_width, na.rm = TRUE),
              mean(out$ci_width_corrected, na.rm = TRUE),
              mean(out$width_inflation[is.finite(out$width_inflation)])))

  d <- lava_gate_diagnostic(rep_$K, rep_$p_ASD, rep_$p_SCZ, r_e = 0.02756,
                            p1_all = a$p_ASD[a$processed %in% TRUE &
                                             is.finite(a$p_ASD) & is.finite(a$p_SCZ)],
                            p2_all = a$p_SCZ[a$processed %in% TRUE &
                                             is.finite(a$p_ASD) & is.finite(a$p_SCZ)])
  cat(sprintf("  diagnostic: usable = %s, failed = {%s}\n",
              d$usable, paste(d$failed_criteria, collapse = ", ")))
  cat(sprintf("    C1 fold-change %.3f | C2 pooled z %.2f / %.2f | C3 null frac %.3f / %.3f\n",
              d$C1$slope_fold_change, d$C2$trait1$pooled_z, d$C2$trait2$pooled_z,
              d$C3$frac_null_trait1, d$C3$frac_null_trait2))
  if (d$usable) { cat("  [FAIL] R diagnostic calls the completed pair usable; ",
                      "the published verdict is NOT APPLICABLE\n"); fail <- fail + 1L
  } else cat("  [ok  ] R diagnostic independently reproduces the published verdict\n")
} else cat("  (skipped: results files not present)\n")

outdir <- file.path(SUB, "results", "crossvalidation")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
writeLines(jsonlite::toJSON(c(list(
    n_param_combinations = nrow(P),
    K_min = min(KS), K_max = max(KS),
    theta_min = min(THETAS), theta_max = max(THETAS),
    r_e_min = min(RES), r_e_max = max(RES),
    n_checks = length(DIFFS), n_failed = fail,
    max_diff_delta = max(DIFFS[["Delta_1"]], DIFFS[["Delta_2"]]),
    max_diff_slope_intercept = max(DIFFS[["slope"]], DIFFS[["intercept"]])),
  DIFFS), auto_unbox = TRUE, pretty = TRUE, digits = 10),
  file.path(outdir, "r_vs_python.json"))

cat(sprintf("\n%s\n", if (fail == 0L) "ALL CHECKS PASSED" else
            sprintf("%d CHECK(S) FAILED", fail)))
quit(status = if (fail == 0L) 0L else 1L)
