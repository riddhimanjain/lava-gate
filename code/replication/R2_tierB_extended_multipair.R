#!/usr/bin/env Rscript
# R2_tierB_extended_multipair.R

suppressMessages({ library(LAVA) })

BASE <- Sys.getenv("LAVAGATE_BASE")
if (!nzchar(BASE)) stop("Set LAVAGATE_BASE to the folder that contains SUBMISSION/")
REF_PREFIX <- file.path(BASE, "Section 2_LD Reference/ukb_ref_consolidated/lava-ukb-v1.1")
LOCFILE    <- file.path(BASE, "LAVA/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
SUMDIR     <- file.path(BASE, "Section 7_LAVA Analysis/sumstats")
INFO       <- file.path(BASE, "Section 7_LAVA Analysis/multipair.info.txt")
OVERLAP    <- file.path(BASE, "Section 7_LAVA Analysis/multipair.sample.overlap.txt")
MPROOT     <- file.path(BASE, "analysis/out_R/multipair")
OUTROOT    <- file.path(BASE, "analysis/out_R/tierB_extended_multipair")

N_LOCI      <- 120
N_REP       <- 2000
RHO_ORIG    <- c(0, 0.2, 0.4, 0.6, 0.8)
RHO_GRID    <- sort(unique(c(RHO_ORIG, -0.8, -0.6, -0.4, -0.2)))
SEED        <- 20260823
UNIV_THRESH <- 0.05

GROUPS <- list(
  P1 = list(phenos = c("SCZ", "BIP"),
            pairs  = list(P1 = c("SCZ", "BIP")), audit_pair = "P1"),
  P2 = list(phenos = c("MDD", "BIP", "MDD_noUKB", "BIP_noUKB"),
            pairs  = list(P2a = c("MDD", "BIP"),
                          P2b = c("MDD_noUKB", "BIP_noUKB")), audit_pair = "P2a"),
  P3 = list(phenos = c("ASD", "ADHD"),
            pairs  = list(P3 = c("ASD", "ADHD")), audit_pair = "P3"),
  P5 = list(phenos = c("AN", "ASD"),
            pairs  = list(P5 = c("AN", "ASD")), audit_pair = "P5")
)

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !(args[1] %in% names(GROUPS)))
  stop(sprintf("usage: Rscript R2_tierB_extended_multipair.R <%s>",
               paste(names(GROUPS), collapse = "|")))
GRP    <- args[1]
PHENOS <- GROUPS[[GRP]]$phenos
PAIRS  <- GROUPS[[GRP]]$pairs
AUDIT  <- file.path(MPROOT, GROUPS[[GRP]]$audit_pair, "attrition_audit.csv")

stopifnot(file.exists(INFO), file.exists(OVERLAP), file.exists(LOCFILE), file.exists(AUDIT))
dir.create(OUTROOT, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

cat("=========================================================\n")
cat("GROUP     :", GRP, "\n")
cat("PHENOTYPES:", paste(PHENOS, collapse = ", "), "\n")
cat("PAIRS     :", paste(sapply(PAIRS, paste, collapse = " x "), collapse = " | "), "\n")
cat("=========================================================\n")

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
                       sample.overlap.file = OVERLAP,
                       ref.prefix          = REF_PREFIX,
                       phenos              = PHENOS,
                       input.dir           = SUMDIR)

clone.locus <- function(lo) list2env(as.list(lo, all.names = TRUE),
                                     envir = new.env(parent = globalenv()))

CKPT <- file.path(OUTROOT, sprintf("R2_checkpoint_%s.rds", GRP))
audit <- data.frame(locus = sel, K = NA_integer_, outcome = NA_character_,
                    stringsAsFactors = FALSE)
rows <- setNames(vector("list", length(PAIRS)), names(PAIRS))
start_ii <- 1
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
  } else {
    K <- lo$K
    audit$K[ii] <- K
    present <- lo$phenos
    audit$outcome[ii] <- "processed"

    for (pn in names(PAIRS)) {
      tr <- PAIRS[[pn]]
      if (is.null(present) || !all(tr %in% present)) next

      Sigma <- as.matrix(lo$sigma)[tr, tr]
      h2    <- pmax(as.numeric(diag(as.matrix(lo$omega)[tr, tr])), 1e-6)
      theta <- h2 / diag(Sigma)
      Lsig  <- chol(Sigma)

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
          sim$phenos <- tr
          sim$delta  <- `dimnames<-`(D, list(NULL, tr))
          sim$omega  <- `dimnames<-`((crossprod(D) / K - Sigma) * lo$nref.scale,
                                     list(tr, tr))
          sim$h2.obs <- setNames(diag(sim$omega), tr)
          sim$sigma  <- `dimnames<-`(Sigma, list(tr, tr))

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

        rows[[pn]][[length(rows[[pn]]) + 1L]] <- data.frame(
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
  }

  saveRDS(list(sel = sel, audit = audit, rows = rows, next_ii = ii + 1L), CKPT)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  done <- ii - start_ii + 1
  cat(sprintf("[%s] %3d/%d  locus %-6d K=%-5s %-24s  %.1f min elapsed, ~%.0f min left\n",
              GRP, ii, length(sel), lid, format(audit$K[ii]), audit$outcome[ii],
              el, el / done * (length(sel) - ii)))
  flush.console()
}

for (pn in names(PAIRS)) {
  tb <- do.call(rbind, rows[[pn]])
  d  <- file.path(OUTROOT, pn); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  if (is.null(tb) || !nrow(tb)) {
    cat(sprintf("%s: no usable loci -- nothing to fit.\n", pn)); next
  }
  utils::write.csv(tb, file.path(d, "tierB_extended.csv"), row.names = FALSE)

  ok <- tb[!is.na(tb$mean_rho_hat), ]
  fit_full <- lm(mean_rho_hat ~ rho_true, data = ok)
  fit_orig <- lm(mean_rho_hat ~ rho_true, data = ok[ok$rho_true %in% RHO_ORIG, ])
  fit_quad <- lm(mean_rho_hat ~ rho_true + I(rho_true^2), data = ok)
  lin_test <- anova(fit_full, fit_quad)

  per_locus <- do.call(rbind, lapply(split(ok, ok$locus), function(g) {
    if (nrow(g) < 3) return(NULL)
    f <- lm(mean_rho_hat ~ rho_true, data = g)
    data.frame(locus = g$locus[1], K = g$K[1], lambda_min = g$lambda_min[1],
               n_cells = nrow(g), slope = coef(f)[2], intercept = coef(f)[1],
               row.names = NULL)
  }))
  fit_mech <- if (nrow(per_locus) >= 3 && length(unique(per_locus$lambda_min)) > 1)
    lm(slope ~ log(lambda_min), data = per_locus) else NULL

  ci  <- confint(fit_full)["rho_true", ]
  cio <- confint(fit_orig)["rho_true", ]

  n_usable_loci <- length(unique(ok$locus))
  summ <- list(
    pair = pn, group = GRP, traits = PAIRS[[pn]],
    n_selected = nrow(audit), n_usable_loci = n_usable_loci,
    slope_full = unname(coef(fit_full)[2]), slope_full_lo = unname(ci[1]),
    slope_full_hi = unname(ci[2]), intercept_full = unname(coef(fit_full)[1]),
    slope_origgrid = unname(coef(fit_orig)[2]), slope_origgrid_lo = unname(cio[1]),
    slope_origgrid_hi = unname(cio[2]), intercept_origgrid = unname(coef(fit_orig)[1]),
    r2_full = summary(fit_full)$r.squared,
    quad_term_p = lin_test$`Pr(>F)`[2],
    per_locus_slope_min = min(per_locus$slope), per_locus_slope_max = max(per_locus$slope),
    mech_slope = if (!is.null(fit_mech)) unname(coef(fit_mech)[2]) else NA_real_,
    mech_p = if (!is.null(fit_mech)) summary(fit_mech)$coefficients[2, 4] else NA_real_)
  writeLines(jsonlite::toJSON(summ, auto_unbox = TRUE, pretty = TRUE, digits = 8),
             file.path(d, "tierB_extended_summary.json"))
  utils::write.csv(per_locus, file.path(d, "tierB_per_locus_slopes.csv"), row.names = FALSE)

  cat(sprintf("\n%s (%s x %s): %d usable loci, slope %.4f (95%% CI %.4f-%.4f), mech_slope %s\n",
              pn, PAIRS[[pn]][1], PAIRS[[pn]][2], n_usable_loci,
              summ$slope_full, ci[1], ci[2],
              if (!is.null(fit_mech)) sprintf("%.4f (P=%.3g)", summ$mech_slope, summ$mech_p) else "NA"))
}

if (file.exists(CKPT)) file.remove(CKPT)
cat("\nDONE ->", OUTROOT, "\n")
