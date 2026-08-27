# 03_applied_summary.R  --  Task 4: the applied ASD x SCZ result

base <- Sys.getenv("LAVAGATE_BASE")
if (!nzchar(base)) stop("Set LAVAGATE_BASE to the folder that contains SUBMISSION/")
outdir <- file.path(base, "analysis/out_R")
df <- read.csv(file.path(outdir, "attrition_audit.csv"))

proc <- df[df$processed, ]
n_proc <- nrow(proc)
rep_ <- df[!is.na(df$rho_hat), ]
n_rep <- nrow(rep_)

sink(file.path(outdir, "applied_summary.txt"))
cat("APPLIED SUMMARY  ASD x SCZ\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

cat("=== A. DENOMINATOR AUDIT (settles the 'negative dependence' question) ===\n")
asd_est <- sum(!is.na(proc$p_ASD)); scz_est <- sum(!is.na(proc$p_SCZ))
both_est <- sum(!is.na(proc$p_ASD) & !is.na(proc$p_SCZ))
cat("Processed loci                          :", n_proc, "\n")
cat("  ASD estimable (p_ASD not NA)          :", asd_est,
    sprintf("(%.1f%% of processed)\n", 100*asd_est/n_proc))
cat("  SCZ estimable (p_SCZ not NA)          :", scz_est,
    sprintf("(%.1f%% of processed)\n", 100*scz_est/n_proc))
cat("  BOTH estimable                        :", both_est,
    sprintf("(%.1f%% of processed)\n", 100*both_est/n_proc))
cat("\nPass rates on the CONDITIONAL denominator (as reported in attrition_audit.txt):\n")
cat("  ASD  P<0.05 | estimable               :", sprintf("%.3f\n", mean(proc$p_ASD < 0.05, na.rm=TRUE)))
cat("  SCZ  P<0.05 | estimable               :", sprintf("%.3f\n", mean(proc$p_SCZ < 0.05, na.rm=TRUE)))
cat("\nPass rates on the COMMON denominator (all processed loci):\n")
a_all <- sum(proc$p_ASD < 0.05, na.rm=TRUE)/n_proc
s_all <- sum(proc$p_SCZ < 0.05, na.rm=TRUE)/n_proc
cat("  ASD  P<0.05 / all processed           :", sprintf("%.4f\n", a_all))
cat("  SCZ  P<0.05 / all processed           :", sprintf("%.4f\n", s_all))
cat("  product (independence prediction)     :", sprintf("%.4f  -> %.1f loci\n",
    a_all*s_all, a_all*s_all*n_proc))
cat("  OBSERVED joint                        :", sprintf("%.4f  -> %d loci\n",
    n_rep/n_proc, n_rep))
cat("\nOn loci where BOTH traits are estimable:\n")
be <- proc[!is.na(proc$p_ASD) & !is.na(proc$p_SCZ), ]
if (nrow(be) > 0) {
  a_be <- mean(be$p_ASD < 0.05); s_be <- mean(be$p_SCZ < 0.05)
  cat("  ASD pass:", sprintf("%.4f", a_be), " SCZ pass:", sprintf("%.4f", s_be),
      " product:", sprintf("%.4f", a_be*s_be),
      " observed:", sprintf("%.4f", mean(be$p_ASD<0.05 & be$p_SCZ<0.05)), "\n")
  tab <- table(ASD = be$p_ASD < 0.05, SCZ = be$p_SCZ < 0.05)
  cat("\n  2x2 contingency (both-estimable loci):\n"); print(tab)
  cat("\n  Fisher exact test for association between the two gates:\n")
  print(fisher.test(tab))
}

cat("\n=== B. PHENOTYPE-DROPOUT (the fourth attrition mode) ===\n")
cat("Processed but ASD dropped               :", sum(is.na(proc$p_ASD)),
    sprintf("(%.1f%%)\n", 100*sum(is.na(proc$p_ASD))/n_proc))
cat("Processed but SCZ dropped               :", sum(is.na(proc$p_SCZ)),
    sprintf("(%.1f%%)\n", 100*sum(is.na(proc$p_SCZ))/n_proc))
cat("Which trait is dropped more often?      :",
    ifelse(sum(is.na(proc$p_ASD)) > sum(is.na(proc$p_SCZ)), "ASD", "SCZ"), "\n")

cat("\n=== C. THE APPLIED RESULT (", n_rep, "reported loci) ===\n")
cat("Mean   rho_hat                          :", sprintf("%.4f\n", mean(rep_$rho_hat)))
cat("Median rho_hat                          :", sprintf("%.4f\n", median(rep_$rho_hat)))
cat("SD     rho_hat                          :", sprintf("%.4f\n", sd(rep_$rho_hat)))
cat("Range                                   :", sprintf("%.4f to %.4f\n",
    min(rep_$rho_hat), max(rep_$rho_hat)))
n_neg <- sum(rep_$rho_hat < 0); n_pos <- sum(rep_$rho_hat > 0)
cat("\nSIGN CONCORDANCE\n")
cat("  negative (diametric direction)        :", n_neg,
    sprintf("(%.1f%%)\n", 100*n_neg/n_rep))
cat("  positive (convergent direction)       :", n_pos,
    sprintf("(%.1f%%)\n", 100*n_pos/n_rep))
cat("  binomial test vs 50/50                :\n")
print(binom.test(n_pos, n_rep, 0.5))

cat("\nSIGN DETERMINACY (does the 95% CI exclude zero?)\n")
det <- !is.na(rep_$rho_lower) & !is.na(rep_$rho_upper) &
       (sign(rep_$rho_lower) == sign(rep_$rho_upper))
cat("  CI excludes zero                      :", sum(det),
    sprintf("(%.1f%% of %d)\n", 100*sum(det)/n_rep, n_rep))
cat("  CI spans zero (sign undetermined)     :", sum(!det),
    sprintf("(%.1f%%)\n", 100*sum(!det)/n_rep))
cat("  mean CI width                         :",
    sprintf("%.4f\n", mean(rep_$rho_upper - rep_$rho_lower, na.rm=TRUE)))
cat("  median CI width                       :",
    sprintf("%.4f\n", median(rep_$rho_upper - rep_$rho_lower, na.rm=TRUE)))
cat("  fraction wider than 1.0               :",
    sprintf("%.3f\n", mean((rep_$rho_upper - rep_$rho_lower) > 1, na.rm=TRUE)))
if (sum(det) > 0) {
  cat("\n  Of the sign-determinate loci:\n")
  cat("    negative:", sum(det & rep_$rho_hat < 0),
      "  positive:", sum(det & rep_$rho_hat > 0), "\n")
}

cat("\nBIVARIATE SIGNIFICANCE\n")
cat("  bivar P<0.05                          :", sum(rep_$bivar_p < 0.05, na.rm=TRUE), "\n")
cat("  bivar P<0.05/", n_rep, " (Bonferroni)  :",
    sum(rep_$bivar_p < 0.05/n_rep, na.rm=TRUE), "\n", sep="")

cat("\n=== D. THE KEY TEST ===\n")
cat("REVISED 17 Aug 2026. This section originally reported ONLY the test in D1 and\n")
cat("concluded from it that the data cannot separate diametric from null. That was\n")
cat("wrong on two counts -- it never tested against the null, and it never tested\n")
cat("DIRECTION. D2 and D3 below are the tests that were missing. The count of\n")
cat("sign-determinate loci is RETIRED as non-discriminating; direction is the\n")
cat("statistic that discriminates. See TASKS_EXECUTION_LOG.md Task 4 section D.\n")

cat("\n--- D1. Determinacy RATE vs a true rho = -0.4 (retired: non-discriminating) ---\n")
cat("Tier A (K=300) predicts 13.9% determinacy at a true rho of -0.4. But it also\n")
cat("predicts 12.3% at +0.4 -- determinacy is nearly SYMMETRIC in the sign of rho,\n")
cat("so this rate cannot tell the two hypotheses apart whatever it comes out at.\n\n")
cat("  Expected sign-determinate loci if diametric were true:",
    sprintf("%.1f of %d\n", 0.139*n_rep, n_rep))
cat("  Observed sign-determinate loci                      :",
    sprintf("%d of %d\n", sum(det), n_rep))
cat("\n  Binomial test, observed vs 13.9% expected under a true rho = -0.4:\n")
print(binom.test(sum(det), n_rep, 0.139))
cat("\n  Read this as uninformative, NOT as support for a null. Agreement on the rate\n")
cat("  while the SIGN is inverted (see D3) is evidence against diametric.\n")

cat("\n--- D2. Determinacy rate vs the NULL (true rho = 0) ---\n")
cat("Under a true rho of 0 a 95% interval excludes zero about 5% of the time. The\n")
cat("original section never ran this comparison, so it reported 'indistinguishable\n")
cat("from either' having tested only one side.\n\n")
print(binom.test(sum(det), n_rep, 0.05))

cat("\n--- D3. DIRECTION of the sign-determinate loci (the discriminating test) ---\n")
n_det_pos <- sum(det & rep_$rho_hat > 0)
n_det_neg <- sum(det & rep_$rho_hat < 0)
cat("A diametric architecture predicts the determinate loci are predominantly\n")
cat("NEGATIVE. A convergent one predicts positive. This is the test that separates\n")
cat("them, and unlike the count it is not symmetric in rho.\n\n")
cat("  sign-determinate and POSITIVE :", n_det_pos, "\n")
cat("  sign-determinate and NEGATIVE :", n_det_neg, "\n")
if (n_det_pos + n_det_neg > 0) {
  cat("\n  Binomial test on direction, vs a neutral 50/50 null:\n")
  print(binom.test(n_det_pos, n_det_pos + n_det_neg, 0.5))
}
cat("\nCAVEAT, and it matters. Every determinate interval has rho_upper pinned at\n")
cat("exactly 1.000 (boundary-truncated). The SIGN is determined; the MAGNITUDE is\n")
cat("entirely unconstrained. Calling these loci determinate is technically correct\n")
cat("and practically hollow, and the paper says so.\n")

cat("\n=== E. ATTENUATION-CORRECTED MAGNITUDES ===\n")
cat("Tier B pooled slope 0.5512 (95% CI 0.5138-0.5886) is the empirical\n")
cat("attenuation factor. Dividing observed estimates by it gives:\n")
cat("  mean rho_hat observed                 :", sprintf("%.4f\n", mean(rep_$rho_hat)))
cat("  mean rho_hat corrected (/0.5512)      :", sprintf("%.4f\n", mean(rep_$rho_hat)/0.5512))
cat("  corrected, CI bounds of slope         :", sprintf("%.4f to %.4f\n",
    mean(rep_$rho_hat)/0.5886, mean(rep_$rho_hat)/0.5138))
cat("\nNOTE: correcting the point estimate does NOT narrow the interval. Interval\n")
cat("width remains the binding constraint on any per-locus claim, not this.\n")

cat("\n=== F. REQUIRED ASD SAMPLE SIZE ===\n")
cat("theta = Omega/Sigma is the local signal-to-noise ratio, and Sigma scales as\n")
cat("1/N, so theta is proportional to N. Current ASD N = 46,350.\n\n")
K_med <- median(proc$K, na.rm=TRUE)
t_crit <- qchisq(0.95, K_med)
cat("Median K =", K_med, "-> chi-square 95% critical value =", sprintf("%.1f\n", t_crit))
cat("\nSolving for the multiple m of current N at which ASD's gate pass rate\n")
cat("matches SCZ's current rate:\n")
asd_rate <- mean(proc$p_ASD < 0.05, na.rm=TRUE)
scz_rate <- mean(proc$p_SCZ < 0.05, na.rm=TRUE)
lam_for <- function(rate, K) {
  f <- function(l) pchisq(qchisq(0.95, K), K, ncp=l, lower.tail=FALSE) - rate
  if (f(0) > 0) return(0)
  uniroot(f, c(0, 50*K))$root
}
lam_asd <- lam_for(asd_rate, K_med); lam_scz <- lam_for(scz_rate, K_med)
cat("  implied ASD noncentrality lambda      :", sprintf("%.2f\n", lam_asd))
cat("  implied SCZ noncentrality lambda      :", sprintf("%.2f\n", lam_scz))
cat("  ratio (required N multiple)           :", sprintf("%.2f\n", lam_scz/max(lam_asd,1e-9)))
cat("  => required ASD N                     :",
    sprintf("%.0f\n", 46350 * lam_scz/max(lam_asd,1e-9)))
cat("\nCAVEAT: this assumes local h2 per block is comparable between the two traits,\n")
cat("which it is not exactly. Treat as an order-of-magnitude figure and say so.\n")
sink()

cat("Done. See", file.path(outdir, "applied_summary.txt"), "\n")
