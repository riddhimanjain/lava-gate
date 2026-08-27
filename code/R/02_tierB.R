# 02_tierB.R  --  Tier B: real loci, real K, real Sigma, real LAVA code path

suppressMessages({library(LAVA); library(MASS)})

base <- Sys.getenv("LAVAGATE_BASE")
if (!nzchar(base)) stop("Set LAVAGATE_BASE to the folder that contains SUBMISSION/")
outdir  <- file.path(base, "analysis/out_R")
sumdir  <- file.path(base, "Section 7_LAVA Analysis/sumstats")
refpref <- file.path(base, "Section 2_LD Reference/ukb_ref_consolidated/lava-ukb-v1.1")
locfile <- file.path(base, "LAVA/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
infofile    <- file.path(base, "Section 7_LAVA Analysis/input.info.txt")
overlapfile <- file.path(base, "Section 7_LAVA Analysis/sample.overlap.txt")

RHO_GRID  <- c(0, 0.2, 0.4, 0.6, 0.8)
N_REP     <- 2000
N_LOCI_USE<- 40
SEED      <- 20260814
set.seed(SEED)

audit <- read.csv(file.path(outdir, "attrition_audit.csv"))
proc  <- audit[audit$processed & !is.na(audit$K), ]
stopifnot(nrow(proc) > 0)

qs  <- quantile(proc$K, probs = seq(0, 1, length.out = N_LOCI_USE + 1))
sel <- unique(sapply(seq_len(N_LOCI_USE), function(j) {
  cand <- proc[proc$K >= qs[j] & proc$K <= qs[j + 1], ]
  if (nrow(cand) == 0) return(NA)
  cand$locus[sample.int(nrow(cand), 1)]
}))
sel <- sel[!is.na(sel)]
cat("Selected", length(sel), "real loci; K range",
    min(proc$K[proc$locus %in% sel]), "-", max(proc$K[proc$locus %in% sel]), "\n")

input <- process.input(input.info.file = infofile,
                       sample.overlap.file = if (file.exists(overlapfile)) overlapfile else NULL,
                       ref.prefix = refpref, phenos = c("ASD", "SCZ"),
                       input.dir = sumdir)
loci  <- read.loci(locfile)

results <- list(); idx <- 1
ckpt <- file.path(outdir, "tierB_checkpoint.rds")

for (lid in sel) {
  lrow  <- loci[loci$LOC == lid, ]
  locus <- tryCatch(process.locus(lrow, input), error = function(e) NULL)
  if (is.null(locus) || is.null(locus$phenos) || length(locus$phenos) != 2) next

  K     <- locus$K
  Sigma <- as.matrix(locus$sigma)
  h2    <- pmax(as.numeric(diag(as.matrix(locus$omega))), 1e-6)
  theta <- h2 / diag(Sigma)

  for (rho in RHO_GRID) {
    Om <- matrix(c(h2[1], rho*sqrt(h2[1]*h2[2]),
                   rho*sqrt(h2[1]*h2[2]), h2[2]), 2, 2)
    ev <- eigen(Om, symmetric = TRUE)$values
    if (min(ev) <= 0) next

    Lom <- chol(Om * K)
    rec <- numeric(0); npass <- 0; nplim <- 0

    for (rep in seq_len(N_REP)) {
      M <- matrix(0, K, 2); M[1:2, ] <- Lom
      D <- M + mvrnorm(K, mu = c(0, 0), Sigma = Sigma)

      sim <- locus
      sim$delta <- D
      S  <- t(D) %*% D
      sim$omega  <- (S/K - Sigma) * locus$nref.scale
      sim$h2.obs <- diag(sim$omega)

      u <- tryCatch(run.univ(sim), error = function(e) NULL)
      if (is.null(u) || nrow(u) != 2 || any(is.na(u$p))) next
      if (!all(u$p < 0.05)) next
      npass <- npass + 1

      b <- tryCatch(run.bivar(sim, CIs = FALSE, p.values = FALSE),
                    error = function(e) NULL)
      if (is.null(b) || nrow(b) == 0) next
      if (is.na(b$rho[1])) { nplim <- nplim + 1; next }
      rec <- c(rec, b$rho[1])
    }

    results[[idx]] <- data.frame(
      locus = lid, K = K, theta1 = theta[1], theta2 = theta[2],
      sigma12 = Sigma[1,2], r_e = Sigma[1,2]/sqrt(Sigma[1,1]*Sigma[2,2]),
      rho_true = rho, n_rep = N_REP, n_gate_pass = npass,
      n_paramlim_drop = nplim, n_reported = length(rec),
      gate_rate = npass/N_REP,
      mean_rho_hat = if (length(rec)) mean(rec) else NA_real_,
      se_rho_hat   = if (length(rec) > 1) sd(rec)/sqrt(length(rec)) else NA_real_
    )
    idx <- idx + 1
  }
  saveRDS(results, ckpt)
  cat("locus", lid, "K =", K, "done\n")
}

tb <- do.call(rbind, results)
write.csv(tb, file.path(outdir, "tierB.csv"), row.names = FALSE)

ok <- tb[!is.na(tb$mean_rho_hat), ]
fit <- lm(mean_rho_hat ~ rho_true, data = ok)

sink(file.path(outdir, "tierB.txt"))
cat("TIER B  --  real loci, real K, real Sigma, real LAVA code path\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")
cat("Loci used:", length(unique(ok$locus)), " K range:",
    min(ok$K), "-", max(ok$K), " median:", median(ok$K), "\n")
cat("Replicates per cell:", N_REP, "\n\n")
print(summary(fit))
cat("\nPooled slope     :", signif(coef(fit)[2], 5), "\n")
cat("Pooled intercept :", signif(coef(fit)[1], 5), "\n")
cat("\nCompare against Tier A (out/RESULTS_TIERA.md):\n")
cat("  Tier A K=300 slope was 0.6064 under the assumed sampling model.\n")
cat("  A materially different Tier B slope means the sampling model, not the\n")
cat("  gate analysis, is what needs revisiting.\n")
cat("\nPer-locus slopes:\n")
for (l in unique(ok$locus)) {
  d <- ok[ok$locus == l, ]
  if (nrow(d) >= 3) {
    f <- lm(mean_rho_hat ~ rho_true, data = d)
    cat(sprintf("  locus %-8s K=%5d  slope=%.4f  intercept=%+.4f\n",
                l, d$K[1], coef(f)[2], coef(f)[1]))
  }
}
sink()

cat("\nDone. See", file.path(outdir, "tierB.txt"), "\n")
