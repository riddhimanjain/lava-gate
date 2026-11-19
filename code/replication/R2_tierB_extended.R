#!/usr/bin/env Rscript
# R2_tierB_extended.R

suppressMessages({ library(LAVA) })

BASE <- Sys.getenv("LAVAGATE_BASE")
if (!nzchar(BASE)) stop("Set LAVAGATE_BASE to the folder that contains SUBMISSION/")
REF_PREFIX <- file.path(BASE, "Section 2_LD Reference/ukb_ref_consolidated/lava-ukb-v1.1")
LOCFILE    <- file.path(BASE, "LAVA/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
SUMDIR     <- file.path(BASE, "Section 7_LAVA Analysis/sumstats")
INFO       <- file.path(BASE, "Section 7_LAVA Analysis/input.info.txt")
OVERLAP    <- file.path(BASE, "Section 7_LAVA Analysis/sample.overlap.txt")
AUDIT      <- file.path(BASE, "analysis/out_R/attrition_audit.csv")
PHENOS     <- c("ASD", "SCZ")
OUTDIR     <- file.path(BASE, "analysis/out_R/tierB_extended")

N_LOCI      <- 120
N_REP       <- 2000
RHO_ORIG    <- c(0, 0.2, 0.4, 0.6, 0.8)
RHO_GRID    <- sort(unique(c(RHO_ORIG, -0.8, -0.6, -0.4, -0.2)))
SEED        <- 20260823
UNIV_THRESH <- 0.05

PUBLISHED_SLOPE     <- 0.5512
PUBLISHED_INTERCEPT <- 0.0283

stopifnot(file.exists(INFO), file.exists(LOCFILE), file.exists(AUDIT))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

audit_in <- utils::read.csv(AUDIT)
frame    <- audit_in[audit_in$processed & !is.na(audit_in$K), ]
stopifnot(nrow(frame) > N_LOCI)

qs <- quantile(frame$K, probs = seq(0, 1, length.out = N_LOCI + 1))
sel <- unique(unlist(lapply(seq_len(N_LOCI), function(j) {
  cand <- frame$locus[frame$K >= qs[j] & frame$K <= qs[j + 1]]
  if (!length(cand)) NA_integer_ else cand[sample.int(length(cand), 1)]
})))
sel <- sort(sel[!is.na(sel)])
cat(sprintf("Selected %d loci; K range %d-%d (median %.0f) from a frame of %d processed loci\n",
            length(sel), min(frame$K[frame$locus %in% sel]),
            max(frame$K[frame$locus %in% sel]),
            median(frame$K[frame$locus %in% sel]), nrow(frame)))

loci  <- read.loci(LOCFILE)
input <- process.input(input.info.file    = INFO,
                       sample.overlap.file = if (file.exists(OVERLAP)) OVERLAP else NULL,
                       ref.prefix          = REF_PREFIX,
                       phenos              = PHENOS,
                       input.dir           = SUMDIR)

clone.locus <- function(lo) list2env(as.list(lo, all.names = TRUE),
                                     envir = new.env(parent = globalenv()))

CKPT <- file.path(OUTDIR, "R2_checkpoint.rds")
audit <- data.frame(locus = sel, K = NA_integer_, outcome = NA_character_,
                    dropped_pheno = NA_character_, theta1 = NA_real_,
                    theta2 = NA_real_, lambda_min = NA_real_,
                    stringsAsFactors = FALSE)
rows <- list(); start_ii <- 1
if (file.exists(CKPT)) {
  s <- readRDS(CKPT)
  if (identical(s$sel, sel)) {
    audit <- s$audit; rows <- s$rows; start_ii <- s$next_ii
    cat("Resuming from checkpoint at locus", start_ii, "of", length(sel), "\n")
  }
}

t0 <- Sys.time()
for (ii in seq(start_ii, length(sel))) {
  lid  <- sel[ii]
  lrow <- loci[loci$LOC == lid, ]
  lo   <- tryCatch(process.locus(lrow, input), error = function(e) NULL)

  if (is.null(lo)) {
    audit$outcome[ii] <- "channel1_not_processed"
  } else if (is.null(lo$phenos) || length(lo$phenos) != 2) {
    audit$K[ii] <- lo$K
    audit$outcome[ii] <- "channel2_phenotype_dropped"
    audit$dropped_pheno[ii] <- setdiff(PHENOS, lo$phenos)[1]
  } else {
    K     <- lo$K
    Sigma <- as.matrix(lo$sigma)[PHENOS, PHENOS]
    h2    <- pmax(as.numeric(diag(as.matrix(lo$omega)[PHENOS, PHENOS])), 1e-6)
    theta <- h2 / diag(Sigma)
    audit$K[ii] <- K; audit$outcome[ii] <- "usable"
    audit$theta1[ii] <- theta[1]; audit$theta2[ii] <- theta[2]
    audit$lambda_min[ii] <- K * min(theta)

    Lsig <- chol(Sigma)

    for (rho in RHO_GRID) {
      Om <- matrix(c(h2[1], rho * sqrt(h2[1] * h2[2]),
                     rho * sqrt(h2[1] * h2[2]), h2[2]), 2, 2)
      if (min(eigen(Om, symmetric = TRUE, only.values = TRUE)$values) <= 0) next
      Lom <- chol(Om * K)

      rec <- numeric(0); npass <- 0L; nplim <- 0L
      for (rep in seq_len(N_REP)) {
        D <- matrix(rnorm(K * 2), K, 2) %*% Lsig
        D[1:2, ] <- D[1:2, ] + Lom

        sim <- clone.locus(lo)
        sim$delta  <- `dimnames<-`(D, list(NULL, PHENOS))
        sim$omega  <- `dimnames<-`((crossprod(D) / K - Sigma) * lo$nref.scale,
                                   list(PHENOS, PHENOS))
        sim$h2.obs <- setNames(diag(sim$omega), PHENOS)

        u <- tryCatch(run.univ(sim), error = function(e) NULL)
        if (is.null(u) || nrow(u) != 2 || any(is.na(u$p))) next
        if (!all(u$p < UNIV_THRESH)) next
        npass <- npass + 1L

        b <- tryCatch(run.bivar(sim, CIs = FALSE, p.values = FALSE),
                      error = function(e) NULL)
        if (is.null(b) || nrow(b) == 0) next
        if (is.na(b$rho[1])) { nplim <- nplim + 1L; next }
        rec <- c(rec, b$rho[1])
      }

      rows[[length(rows) + 1L]] <- data.frame(
        locus = lid, K = K, theta1 = theta[1], theta2 = theta[2],
        lambda_min = K * min(theta),
        r_e = Sigma[1, 2] / sqrt(Sigma[1, 1] * Sigma[2, 2]),
        rho_true = rho, n_rep = N_REP, n_gate_pass = npass,
        n_paramlim_drop = nplim, n_reported = length(rec),
        gate_rate = npass / N_REP,
        mean_rho_hat = if (length(rec)) mean(rec) else NA_real_,
        se_rho_hat   = if (length(rec) > 1) sd(rec) / sqrt(length(rec)) else NA_real_,
        stringsAsFactors = FALSE)
    }
  }

  saveRDS(list(sel = sel, audit = audit, rows = rows, next_ii = ii + 1L), CKPT)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  done <- ii - start_ii + 1
  cat(sprintf("  %3d/%d  locus %-6d K=%-5s %-26s  %.1f min elapsed, ~%.0f min left\n",
              ii, length(sel), lid, format(audit$K[ii]), audit$outcome[ii],
              el, el / done * (length(sel) - ii)))
  flush.console()
  if (length(rows))
    utils::write.csv(do.call(rbind, rows),
                     file.path(OUTDIR, "tierB_extended.csv"), row.names = FALSE)
}

tb <- do.call(rbind, rows)
utils::write.csv(tb,    file.path(OUTDIR, "tierB_extended.csv"),   row.names = FALSE)
utils::write.csv(audit, file.path(OUTDIR, "tierB_locus_audit.csv"), row.names = FALSE)

ok <- tb[!is.na(tb$mean_rho_hat), ]
fit_full <- lm(mean_rho_hat ~ rho_true, data = ok)
fit_orig <- lm(mean_rho_hat ~ rho_true, data = ok[ok$rho_true %in% RHO_ORIG, ])

fit_quad <- lm(mean_rho_hat ~ rho_true + I(rho_true^2), data = ok)
lin_test <- anova(fit_full, fit_quad)

per_locus <- do.call(rbind, lapply(split(ok, ok$locus), function(d) {
  if (nrow(d) < 3) return(NULL)
  f <- lm(mean_rho_hat ~ rho_true, data = d)
  data.frame(locus = d$locus[1], K = d$K[1], lambda_min = d$lambda_min[1],
             n_cells = nrow(d), slope = coef(f)[2], intercept = coef(f)[1],
             row.names = NULL)
}))
fit_mech <- lm(slope ~ log(lambda_min), data = per_locus)

ci  <- confint(fit_full)["rho_true", ]
cio <- confint(fit_orig)["rho_true", ]

summ <- list(
  n_selected = nrow(audit),
  n_usable   = sum(audit$outcome == "usable", na.rm = TRUE),
  n_channel1 = sum(audit$outcome == "channel1_not_processed", na.rm = TRUE),
  n_channel2 = sum(audit$outcome == "channel2_phenotype_dropped", na.rm = TRUE),
  survivorship_loss_pct = 100 * mean(audit$outcome != "usable", na.rm = TRUE),
  slope_full = unname(coef(fit_full)[2]), slope_full_lo = unname(ci[1]),
  slope_full_hi = unname(ci[2]), intercept_full = unname(coef(fit_full)[1]),
  slope_origgrid = unname(coef(fit_orig)[2]), slope_origgrid_lo = unname(cio[1]),
  slope_origgrid_hi = unname(cio[2]), intercept_origgrid = unname(coef(fit_orig)[1]),
  r2_full = summary(fit_full)$r.squared,
  quad_term_p = lin_test$`Pr(>F)`[2],
  per_locus_slope_min = min(per_locus$slope), per_locus_slope_max = max(per_locus$slope),
  mech_slope = unname(coef(fit_mech)[2]),
  mech_p = summary(fit_mech)$coefficients[2, 4],
  published_slope = PUBLISHED_SLOPE,
  published_in_extended_CI = PUBLISHED_SLOPE >= ci[1] && PUBLISHED_SLOPE <= ci[2])
writeLines(jsonlite::toJSON(summ, auto_unbox = TRUE, pretty = TRUE, digits = 8),
           file.path(OUTDIR, "tierB_extended_summary.json"))

sink(file.path(OUTDIR, "tierB_extended.txt"))
cat("TIER B (EXTENDED) --", format(Sys.time()), "\n")
cat("Real loci, real K, real Sigma, LAVA's own run.univ()/run.bivar() path.\n\n")
cat(sprintf("Loci selected      : %d   (stratified across measured K in attrition_audit.csv)\n",
            nrow(audit)))
print(table(audit$outcome, useNA = "ifany"))
cat(sprintf("\nSURVIVORSHIP: %.1f%% of selected loci were lost before Tier B could use\n",
            summ$survivorship_loss_pct))
cat("them. The published 18-locus slope was fitted on survivors only.\n")
cat("Channel 2 drops, by trait:\n"); print(table(audit$dropped_pheno, useNA = "no"))
cat(sprintf("\nUsable loci        : %d   K range %d - %d (median %.0f)\n",
            summ$n_usable, min(ok$K), max(ok$K), median(ok$K)))
cat(sprintf("Cells fitted       : %d over rho grid {%s}\n",
            nrow(ok), paste(RHO_GRID, collapse = ", ")))

cat("\n--- FIT ON THE ORIGINAL GRID (directly comparable to 0.5512) ---\n")
print(summary(fit_orig))
cat(sprintf("slope %.5f  95%% CI %.4f - %.4f   intercept %.5f\n",
            summ$slope_origgrid, cio[1], cio[2], summ$intercept_origgrid))

cat("\n--- FIT ON THE FULL SYMMETRIC GRID ---\n")
print(summary(fit_full))
cat(sprintf("slope %.5f  95%% CI %.4f - %.4f   intercept %.5f   R2 %.4f\n",
            summ$slope_full, ci[1], ci[2], summ$intercept_full, summ$r2_full))

cat("\n--- EXACT LINEARITY IN RHO ---\n")
cat("The closed form claims E[rho.hat | gate] is exactly linear in rho. Adding a\n")
cat("quadratic term to the pooled fit should not improve it.\n")
print(lin_test)
cat(sprintf("quadratic term P = %.4g  (large P supports exact linearity)\n", summ$quad_term_p))

cat("\n--- SURVIVORSHIP MECHANISM ---\n")
cat("Per-locus slope regressed on log of the weaker trait's noncentrality K*theta.\n")
cat("Channel 2 removes LOW-noncentrality loci. If this slope is POSITIVE, then\n")
cat("restricting to survivors raises the estimated slope, i.e. UNDERSTATES\n")
cat("attenuation -- which is precisely the caveat on 0.5512.\n")
print(summary(fit_mech))
cat(sprintf("mechanism slope %+.5f  P = %.4g\n", summ$mech_slope, summ$mech_p))
cat(sprintf("Per-locus slope range: %.4f - %.4f\n",
            summ$per_locus_slope_min, summ$per_locus_slope_max))

cat("\n--- VERDICT AGAINST THE PUBLISHED VALUE ---\n")
cat(sprintf("Published (18 loci) : slope %.4f, intercept %.4f\n",
            PUBLISHED_SLOPE, PUBLISHED_INTERCEPT))
cat(sprintf("Extended (%d loci)  : slope %.4f, 95%% CI %.4f - %.4f\n",
            summ$n_usable, summ$slope_full, ci[1], ci[2]))
cat(sprintf("Published value lies %s the extended 95%% CI.\n",
            if (summ$published_in_extended_CI) "INSIDE" else "OUTSIDE"))
if (!summ$published_in_extended_CI)
  cat("=> The manuscript's headline slope MUST be updated to the extended value,\n   in section 3.2, Table 8 and the Conclusion, before anything else is done.\n")
sink()

utils::write.csv(per_locus, file.path(OUTDIR, "tierB_per_locus_slopes.csv"), row.names = FALSE)
if (file.exists(CKPT)) file.remove(CKPT)
cat("\nDONE ->", OUTDIR, "\n")
writeLines(readLines(file.path(OUTDIR, "tierB_extended.txt")))
