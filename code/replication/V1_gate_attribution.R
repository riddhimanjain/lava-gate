#!/usr/bin/env Rscript
# V1_gate_attribution.R   --   ADVERSARIAL VERIFICATION OF THE PAPER'S CLAIM

.args  <- commandArgs(trailingOnly = FALSE)
.self  <- sub("^--file=", "", grep("^--file=", .args, value = TRUE)[1])
BASE   <- normalizePath(file.path(dirname(.self), "..", ".."))
OUTDIR <- file.path(BASE, "results/verification")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

N_REP    <- 4000
RHO_GRID <- c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8)
ALPHA    <- 0.05
PARAM_LIM <- 1.25
CHUNK    <- 500
SEED     <- 20260826

SRC <- list(
  P4  = file.path(BASE, "results/tierB_extended/tierB_extended.csv"),
  P1  = file.path(BASE, "results/tierB_extended_multipair/P1/tierB_extended.csv"),
  P2a = file.path(BASE, "results/tierB_extended_multipair/P2a/tierB_extended.csv"),
  P2b = file.path(BASE, "results/tierB_extended_multipair/P2b/tierB_extended.csv"),
  P3  = file.path(BASE, "results/tierB_extended_multipair/P3/tierB_extended.csv"),
  P5  = file.path(BASE, "results/tierB_extended_multipair/P5/tierB_extended.csv")
)

set.seed(SEED)

sim_cell <- function(K, th1, th2, r_e, rho, n_rep = N_REP) {
  Sigma <- matrix(c(1, r_e, r_e, 1), 2, 2)
  Om    <- matrix(c(th1, rho * sqrt(th1 * th2),
                    rho * sqrt(th1 * th2), th2), 2, 2)
  if (min(eigen(Om, symmetric = TRUE, only.values = TRUE)$values) <= 0) return(NULL)

  U   <- chol(Sigma)
  Lom <- chol(Om * K)
  thr <- qchisq(1 - ALPHA, K)

  acc <- list(A = numeric(0), B = numeric(0), C = numeric(0), D = numeric(0))
  n_est <- 0L; n_gate <- 0L; n_cap <- 0L; n_plim <- 0L

  done <- 0L
  while (done < n_rep) {
    m  <- min(CHUNK, n_rep - done); done <- done + m
    Z1 <- matrix(rnorm(m * K), m, K); Z2 <- matrix(rnorm(m * K), m, K)
    D1 <- Z1 * U[1, 1]
    D2 <- Z1 * U[1, 2] + Z2 * U[2, 2]
    D1[, 1] <- D1[, 1] + Lom[1, 1]; D1[, 2] <- D1[, 2] + Lom[2, 1]
    D2[, 1] <- D2[, 1] + Lom[1, 2]; D2[, 2] <- D2[, 2] + Lom[2, 2]

    S11 <- rowSums(D1 * D1); S22 <- rowSums(D2 * D2); S12 <- rowSums(D1 * D2)
    rm(Z1, Z2, D1, D2)

    o11 <- S11 / K - Sigma[1, 1]
    o22 <- S22 / K - Sigma[2, 2]
    o12 <- S12 / K - Sigma[1, 2]

    estimable <- o11 > 0 & o22 > 0
    gated     <- (S11 / Sigma[1, 1]) > thr & (S22 / Sigma[2, 2]) > thr

    rh <- rep(NA_real_, m)
    rh[estimable] <- o12[estimable] / sqrt(o11[estimable] * o22[estimable])

    rt <- rh
    plim <- !is.na(rt) & abs(rt) > PARAM_LIM
    rt[plim] <- NA_real_
    capped <- !is.na(rt) & abs(rt) > 1
    rt[capped] <- sign(rt[capped])

    n_est  <- n_est  + sum(estimable)
    n_gate <- n_gate + sum(gated)
    n_cap  <- n_cap  + sum(capped & gated)
    n_plim <- n_plim + sum(plim & gated)

    acc$A <- c(acc$A, rh[estimable])
    acc$B <- c(acc$B, rt[estimable & !is.na(rt)])
    acc$C <- c(acc$C, rh[gated])
    acc$D <- c(acc$D, rt[gated & !is.na(rt)])
  }

  data.frame(K = K, theta1 = th1, theta2 = th2, r_e = r_e, rho_true = rho,
             n_rep = n_rep,
             frac_estimable = n_est / n_rep, frac_gated = n_gate / n_rep,
             frac_capped_of_gated  = n_cap  / max(n_gate, 1L),
             frac_plim_of_gated    = n_plim / max(n_gate, 1L),
             mean_A = mean(acc$A), mean_B = mean(acc$B),
             mean_C = mean(acc$C), mean_D = mean(acc$D),
             med_A  = median(acc$A), med_C = median(acc$C),
             n_A = length(acc$A), n_B = length(acc$B),
             n_C = length(acc$C), n_D = length(acc$D),
             stringsAsFactors = FALSE)
}

fit_arm <- function(d, col) {
  ok <- d[is.finite(d[[col]]), ]
  if (nrow(ok) < 5) return(c(NA, NA, NA, NA))
  f <- lm(ok[[col]] ~ ok$rho_true)
  ci <- confint(f)[2, ]
  c(unname(coef(f)[2]), unname(ci[1]), unname(ci[2]), unname(coef(f)[1]))
}

all_rows <- list(); summ <- list()
for (pair in names(SRC)) {
  if (!file.exists(SRC[[pair]])) { cat("MISSING:", SRC[[pair]], "\n"); next }
  tb <- utils::read.csv(SRC[[pair]])
  loci <- unique(tb[, c("locus", "K", "theta1", "theta2", "r_e")])
  loci <- loci[!duplicated(loci$locus), ]
  usable <- unique(tb$locus[!is.na(tb$mean_rho_hat)])
  loci <- loci[loci$locus %in% usable, ]
  cat(sprintf("\n=== %s : %d loci, K %d-%d ===\n", pair, nrow(loci),
              min(loci$K), max(loci$K))); flush.console()

  CKPT <- file.path(OUTDIR, paste0("V1_ckpt_", pair, ".rds"))
  rows <- list(); start_i <- 1L
  if (file.exists(CKPT)) {
    s <- readRDS(CKPT)
    if (identical(s$loci_id, loci$locus)) {
      rows <- s$rows; start_i <- s$next_i
      cat(sprintf("  resuming %s at locus %d/%d\n", pair, start_i, nrow(loci)))
    }
  }
  t0 <- Sys.time()
  for (i in seq(start_i, nrow(loci))) {
    for (rho in RHO_GRID) {
      r <- sim_cell(loci$K[i], loci$theta1[i], loci$theta2[i], loci$r_e[i], rho)
      if (!is.null(r)) { r$pair <- pair; r$locus <- loci$locus[i]
                         rows[[length(rows) + 1L]] <- r }
    }
    saveRDS(list(loci_id = loci$locus, rows = rows, next_i = i + 1L), CKPT)
    if (i %% 5 == 0 || i == nrow(loci))
      cat(sprintf("  %s %d/%d  (%.1f min)\n", pair, i, nrow(loci),
                  as.numeric(difftime(Sys.time(), t0, units = "mins")))) ; flush.console()
  }
  d <- do.call(rbind, rows); all_rows[[pair]] <- d

  fA <- fit_arm(d, "mean_A"); fB <- fit_arm(d, "mean_B")
  fC <- fit_arm(d, "mean_C"); fD <- fit_arm(d, "mean_D")
  summ[[pair]] <- data.frame(
    pair = pair, n_loci = nrow(loci), n_cells = nrow(d),
    slope_A_estimable_raw   = fA[1], A_lo = fA[2], A_hi = fA[3], int_A = fA[4],
    slope_B_estimable_trunc = fB[1], B_lo = fB[2], B_hi = fB[3], int_B = fB[4],
    slope_C_gated_raw       = fC[1], C_lo = fC[2], C_hi = fC[3], int_C = fC[4],
    slope_D_gated_trunc     = fD[1], D_lo = fD[2], D_hi = fD[3], int_D = fD[4],
    gate_effect_at_trunc    = fD[1] - fB[1],
    gate_effect_no_trunc    = fC[1] - fA[1],
    trunc_effect_at_gate    = fD[1] - fC[1],
    trunc_effect_no_gate    = fB[1] - fA[1],
    mean_frac_capped_of_gated = mean(d$frac_capped_of_gated, na.rm = TRUE),
    mean_frac_plim_of_gated   = mean(d$frac_plim_of_gated,   na.rm = TRUE),
    mean_frac_gated           = mean(d$frac_gated, na.rm = TRUE),
    mean_frac_estimable       = mean(d$frac_estimable, na.rm = TRUE),
    stringsAsFactors = FALSE)
  print(summ[[pair]][, c("pair", "slope_A_estimable_raw", "slope_B_estimable_trunc",
                         "slope_C_gated_raw", "slope_D_gated_trunc",
                         "gate_effect_at_trunc", "trunc_effect_at_gate")])
  flush.console()
  utils::write.csv(do.call(rbind, all_rows),
                   file.path(OUTDIR, "V1_cells.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, summ),
                   file.path(OUTDIR, "V1_summary.csv"), row.names = FALSE)
}

cat("\n================ GATE ATTRIBUTION ================\n")
print(do.call(rbind, summ)[, c("pair", "slope_A_estimable_raw",
                               "slope_B_estimable_trunc", "slope_C_gated_raw",
                               "slope_D_gated_trunc", "gate_effect_at_trunc",
                               "trunc_effect_at_gate")])
cat("\nDONE ->", OUTDIR, "\n")
