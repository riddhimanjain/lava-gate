#!/usr/bin/env Rscript
# V0_exactness_check.R  --  prove the fast reimplementation used by V1 is the

suppressMessages(library(LAVA))
set.seed(7)

.args  <- commandArgs(trailingOnly = FALSE)
.self  <- sub("^--file=", "", grep("^--file=", .args, value = TRUE)[1])
BASE   <- normalizePath(file.path(dirname(.self), "..", ".."))
OUTDIR <- file.path(BASE, "results/verification")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
PH <- c("T1", "T2")

make_locus <- function(K, Sigma, D, N = c(50000, 130000)) {
  lo <- new.env(parent = globalenv())
  lo$id <- 1L; lo$K <- K; lo$P <- 2L; lo$phenos <- PH
  lo$binary <- setNames(c(TRUE, TRUE), PH)
  lo$N <- setNames(N, PH)
  lo$nref.scale <- 1
  lo$sigma <- `dimnames<-`(Sigma, list(PH, PH))
  lo$delta <- `dimnames<-`(D, list(NULL, PH))
  lo$omega <- `dimnames<-`(crossprod(D) / K - Sigma, list(PH, PH))
  lo$omega.cor <- suppressWarnings(cov2cor(lo$omega))
  lo$h2.obs <- setNames(diag(lo$omega), PH)
  lo$h2.latent <- setNames(c(NA_real_, NA_real_), PH)
  lo$ascertained.h2 <- setNames(c(FALSE, FALSE), PH)
  lo
}

cfgs <- list(
  list(K = 100, th1 = 0.02,  th2 = 0.15, r = 0.00),
  list(K = 270, th1 = 0.006, th2 = 0.08, r = 0.02756),
  list(K = 500, th1 = 0.05,  th2 = 0.05, r = 0.20),
  list(K =  39, th1 = 0.30,  th2 = 0.01, r = 0.10),
  list(K = 584, th1 = 0.001, th2 = 0.40, r = 0.2372)
)
RHOS <- c(-0.8, -0.3, 0, 0.45, 0.9)
N_REP <- 200

rows <- list()
worst_p <- 0; worst_rho <- 0; n_cmp_p <- 0L; n_cmp_rho <- 0L
gate_not_subset <- 0L; n_gate <- 0L
n_bivar_abort <- 0L; n_abort_but_estimable <- 0L

for (cf in cfgs) {
  K <- cf$K
  Sigma <- matrix(c(1, cf$r, cf$r, 1), 2, 2)
  U <- chol(Sigma)
  thr <- qchisq(0.95, K)
  for (rho in RHOS) {
    Om <- matrix(c(cf$th1, rho * sqrt(cf$th1 * cf$th2),
                   rho * sqrt(cf$th1 * cf$th2), cf$th2), 2, 2)
    if (min(eigen(Om, symmetric = TRUE, only.values = TRUE)$values) <= 0) next
    Lom <- chol(Om * K)
    for (rep in seq_len(N_REP)) {
      D <- matrix(rnorm(K * 2), K, 2) %*% U
      D[1:2, ] <- D[1:2, ] + Lom
      lo <- make_locus(K, Sigma, D)

      u <- run.univ(lo)
      b <- tryCatch(suppressMessages(run.bivar(lo, CIs = FALSE, p.values = FALSE)),
                    error = function(e) NULL)

      S   <- crossprod(D)
      p_d <- pchisq(diag(S) / diag(Sigma) * lo$nref.scale, K, lower.tail = FALSE)
      o   <- S / K - Sigma
      estimable_d <- o[1, 1] > 0 && o[2, 2] > 0
      if (is.null(b)) {
        n_bivar_abort <- n_bivar_abort + 1L
        if (estimable_d) n_abort_but_estimable <- n_abort_but_estimable + 1L
      }
      rho_d <- if (estimable_d) o[1, 2] / sqrt(o[1, 1] * o[2, 2]) else NA_real_
      p_d   <- signif(p_d, 6)
      rho_d <- signif(rho_d, 6)
      rho_dt <- rho_d
      if (!is.na(rho_dt) && abs(rho_dt) > 1.25) rho_dt <- NA_real_
      if (!is.na(rho_dt) && abs(rho_dt) > 1)    rho_dt <- sign(rho_dt)

      dp <- max(abs(u$p - p_d) / pmax(abs(p_d), 1e-300))
      worst_p <- max(worst_p, dp); n_cmp_p <- n_cmp_p + 2L
      dr <- NA_real_
      rho_lava <- if (is.null(b)) NA_real_ else b$rho[1]
      if (estimable_d) {
        if (is.na(rho_lava) != is.na(rho_dt)) {
          dr <- Inf
        } else if (!is.na(rho_lava)) {
          dr <- abs(rho_lava - rho_dt) / max(abs(rho_dt), 1e-12)
          worst_rho <- max(worst_rho, dr); n_cmp_rho <- n_cmp_rho + 1L
        }
      }
      g <- all(diag(S) / diag(Sigma) > thr)
      if (g) { n_gate <- n_gate + 1L
               if (!(o[1, 1] > 0 && o[2, 2] > 0)) gate_not_subset <- gate_not_subset + 1L }

      rows[[length(rows) + 1L]] <- data.frame(
        K = K, r_e = cf$r, rho_true = rho,
        p_lava_1 = u$p[1], p_direct_1 = p_d[1],
        p_lava_2 = u$p[2], p_direct_2 = p_d[2],
        rho_lava = rho_lava, rho_direct = rho_dt, estimable = estimable_d,
        bivar_aborted = is.null(b),
        rel_p = dp, rel_rho = dr, gated = g, stringsAsFactors = FALSE)
    }
  }
  cat(sprintf("  K=%-4d r_e=%-7.5f done\n", K, cf$r)); flush.console()
}

res <- do.call(rbind, rows)
utils::write.csv(res, file.path(OUTDIR, "V0_exactness.csv"), row.names = FALSE)

cat("\n================= V0 EXACTNESS CHECK =================\n")
cat(sprintf("replicates compared          : %d\n", nrow(res)))
cat(sprintf("univariate p, worst rel diff : %.3e  (over %d comparisons)\n",
            worst_p, n_cmp_p))
cat(sprintf("rho.hat,      worst rel diff : %.3e  (over %d comparisons)\n",
            worst_rho, n_cmp_rho))
cat(sprintf("NA/non-NA disagreements      : %d\n", sum(is.infinite(res$rel_rho))))
cat(sprintf("gate passes                  : %d\n", n_gate))
cat(sprintf("gate passes NOT estimable    : %d   (must be 0: qchisq(.95,K) > K)\n",
            gate_not_subset))
cat(sprintf("run.bivar() aborted          : %d of %d  (negative local variance)\n",
            n_bivar_abort, nrow(res)))
cat(sprintf("  ... of which WERE estimable: %d   (must be 0)\n", n_abort_but_estimable))

ok <- worst_p < 1e-12 && worst_rho < 1e-12 &&
      sum(is.infinite(res$rel_rho), na.rm = TRUE) == 0 &&
      gate_not_subset == 0 && n_abort_but_estimable == 0
cat(sprintf("\nVERDICT: %s\n", if (ok)
  "PASS -- the direct formulas ARE LAVA's; V1 is entitled to use them."
  else "FAIL -- do NOT trust V1; the reimplementation differs from LAVA."))
writeLines(as.character(ok), file.path(OUTDIR, "V0_PASS.txt"))
if (!ok) quit(status = 1)
