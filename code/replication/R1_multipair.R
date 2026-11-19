#!/usr/bin/env Rscript
# R1_multipair.R  --  the trait-pair panel, one pass per PHENOTYPE GROUP.

suppressMessages({ library(LAVA) })

BASE <- Sys.getenv("LAVAGATE_BASE")
if (!nzchar(BASE)) stop("Set LAVAGATE_BASE to the folder that contains SUBMISSION/")
REF_PREFIX <- file.path(BASE, "Section 2_LD Reference/ukb_ref_consolidated/lava-ukb-v1.1")
LOCFILE    <- file.path(BASE, "LAVA/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
SUMDIR     <- file.path(BASE, "Section 7_LAVA Analysis/sumstats")
OUTROOT    <- file.path(BASE, "analysis/out_R/multipair")

INFO_FILE    <- file.path(BASE, "Section 7_LAVA Analysis/multipair.info.txt")
OVERLAP_FILE <- file.path(BASE, "Section 7_LAVA Analysis/multipair.sample.overlap.txt")

UNIV_THRESH <- 0.05
CKPT_EVERY  <- 25

GROUPS <- list(
  P1 = list(phenos = c("SCZ", "BIP"),
            pairs  = list(P1 = c("SCZ", "BIP"))),
  P2 = list(phenos = c("MDD", "BIP", "MDD_noUKB", "BIP_noUKB"),
            pairs  = list(P2a = c("MDD", "BIP"),
                          P2b = c("MDD_noUKB", "BIP_noUKB"))),
  P3 = list(phenos = c("ASD", "ADHD"),
            pairs  = list(P3 = c("ASD", "ADHD"))),
  P5 = list(phenos = c("AN", "ASD"),
            pairs  = list(P5 = c("AN", "ASD")))
)

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !(args[1] %in% names(GROUPS)))
  stop(sprintf("usage: Rscript R1_multipair.R <%s> [N_LOCI_MAX]",
               paste(names(GROUPS), collapse = "|")))
GRP        <- args[1]
N_LOCI_MAX <- if (length(args) > 1) as.numeric(args[2]) else Inf
PHENOS     <- GROUPS[[GRP]]$phenos
PAIRS      <- GROUPS[[GRP]]$pairs

cat("=========================================================\n")
cat("GROUP     :", GRP, "\n")
cat("PHENOTYPES:", paste(PHENOS, collapse = ", "), "\n")
cat("PAIRS     :", paste(sapply(PAIRS, paste, collapse = " x "), collapse = " | "), "\n")
if (is.finite(N_LOCI_MAX)) cat("PILOT     : first", N_LOCI_MAX, "loci only\n")
cat("=========================================================\n")

stopifnot(file.exists(INFO_FILE), file.exists(OVERLAP_FILE), file.exists(LOCFILE))
dir.create(OUTROOT, showWarnings = FALSE, recursive = TRUE)

loci <- read.loci(LOCFILE)
n_all <- nrow(loci)
if (is.finite(N_LOCI_MAX)) loci <- loci[seq_len(min(N_LOCI_MAX, nrow(loci))), ]

inp <- process.input(input.info.file     = INFO_FILE,
                     sample.overlap.file = OVERLAP_FILE,
                     ref.prefix          = REF_PREFIX,
                     phenos              = PHENOS,
                     input.dir           = SUMDIR)

N_SNPS_COMMON <- length(inp$sum.stats[[PHENOS[1]]]$SNP)
cat(sprintf("\nSNPs common to %s + reference: %d\n", paste(PHENOS, collapse=" + "),
            N_SNPS_COMMON))

Nmin       <- min(inp$info$N[inp$info$phenotype %in% PHENOS], na.rm = TRUE)
MAX_PROP_K <- 0.75
K_CEILING  <- 700
stopifnot("max.prop.K would bind: K is being capped by the smallest GWAS." =
            floor(MAX_PROP_K * Nmin) > K_CEILING,
          "K/N rule would fire: h2.obs will be NA'd for the smallest GWAS." =
            (K_CEILING / Nmin) < 0.1)
cat(sprintf("Combined-pass checks OK. min(N)=%d -> max.K=%d, max K/N=%.4f\n\n",
            Nmin, floor(MAX_PROP_K * Nmin), K_CEILING / Nmin))

n <- nrow(loci)
blank <- function(ph) {
  d <- data.frame(
    locus = loci$LOC, chr = loci$CHR, start = loci$START, stop = loci$STOP,
    processed = FALSE, K = NA_integer_, n_snps = NA_integer_,
    h2_1 = NA_real_, h2_2 = NA_real_, p_1 = NA_real_, p_2 = NA_real_,
    gate_pass = FALSE, rho_hat = NA_real_, rho_raw = NA_real_,
    paramlim_drop = NA, rho_lower = NA_real_, rho_upper = NA_real_,
    bivar_p = NA_real_, n_pheno_returned = NA_integer_,
    dropped_pheno = NA_character_, stringsAsFactors = FALSE)
  names(d)[names(d) == "h2_1"] <- paste0("h2_", ph[1])
  names(d)[names(d) == "h2_2"] <- paste0("h2_", ph[2])
  names(d)[names(d) == "p_1"]  <- paste0("p_",  ph[1])
  names(d)[names(d) == "p_2"]  <- paste0("p_",  ph[2])
  d
}

CKPT <- file.path(OUTROOT, sprintf("checkpoint_%s.rds", GRP))
res <- lapply(PAIRS, blank)
survivors <- data.frame(locus = loci$LOC, K = NA_integer_, n_snps = NA_integer_,
                        present = NA_character_, secs = NA_real_,
                        stringsAsFactors = FALSE)
start_i <- 1L
if (file.exists(CKPT)) {
  s <- readRDS(CKPT)
  if (identical(s$phenos, PHENOS) && identical(nrow(s$survivors), nrow(survivors))) {
    res <- s$res; survivors <- s$survivors; start_i <- s$next_i
    cat(sprintf("Resuming from checkpoint at locus %d of %d (%.1f%% done)\n",
                start_i, n, 100 * (start_i - 1) / n))
  } else cat("Checkpoint does not match this group; starting fresh.\n")
}

save_all <- function(final = FALSE) {
  for (pn in names(PAIRS)) {
    d <- file.path(OUTROOT, pn); dir.create(d, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(res[[pn]], file.path(d, "attrition_audit.csv"), row.names = FALSE)
  }
  utils::write.csv(survivors, file.path(OUTROOT, sprintf("channel2_survivors_%s.csv", GRP)),
                   row.names = FALSE)
}

t0 <- Sys.time()
for (i in seq(start_i, n)) {
  ti <- Sys.time()
  lo <- tryCatch(process.locus(loci[i, ], inp), error = function(e) NULL)

  if (!is.null(lo)) {
    survivors$K[i] <- lo$K
    survivors$n_snps[i] <- if (!is.null(lo$n.snps)) lo$n.snps else length(lo$snps)
    survivors$present[i] <- paste(lo$phenos, collapse = ",")
    u <- tryCatch(run.univ(lo), error = function(e) NULL)

    for (pn in names(PAIRS)) {
      ph <- PAIRS[[pn]]; r <- res[[pn]]
      cn_h <- paste0("h2_", ph); cn_p <- paste0("p_", ph)
      r$processed[i] <- TRUE
      r$K[i] <- lo$K
      r$n_snps[i] <- survivors$n_snps[i]

      present <- ph[ph %in% lo$phenos]
      r$n_pheno_returned[i] <- length(present)
      if (!is.null(u)) for (k in 1:2) {
        j <- which(u$phen == ph[k])
        if (length(j)) { r[[cn_h[k]]][i] <- u$h2.obs[j[1]]; r[[cn_p[k]]][i] <- u$p[j[1]] }
      }

      if (length(present) < 2) {
        r$dropped_pheno[i] <- paste(setdiff(ph, present), collapse = ",")
        res[[pn]] <- r; next
      }
      pass <- isTRUE(r[[cn_p[1]]][i] < UNIV_THRESH) && isTRUE(r[[cn_p[2]]][i] < UNIV_THRESH)
      if (!pass) { res[[pn]] <- r; next }
      r$gate_pass[i] <- TRUE

      b_raw <- tryCatch(run.bivar(lo, phenos = ph, param.lim = Inf, CIs = FALSE,
                                  p.values = FALSE, cap.estimates = FALSE),
                        error = function(e) NULL)
      if (!is.null(b_raw) && nrow(b_raw) > 0) r$rho_raw[i] <- b_raw$rho[1]

      b <- tryCatch(run.bivar(lo, phenos = ph), error = function(e) NULL)
      if (!is.null(b) && nrow(b) > 0) {
        r$rho_hat[i]   <- b$rho[1]
        r$rho_lower[i] <- b$rho.lower[1]
        r$rho_upper[i] <- b$rho.upper[1]
        r$bivar_p[i]   <- b$p[1]
        r$paramlim_drop[i] <- is.na(b$rho[1]) && !is.na(r$rho_raw[i])
      }
      res[[pn]] <- r
    }
  }
  survivors$secs[i] <- as.numeric(difftime(Sys.time(), ti, units = "secs"))

  if (i %% CKPT_EVERY == 0 || i == n) {
    saveRDS(list(phenos = PHENOS, res = res, survivors = survivors, next_i = i + 1L), CKPT)
    save_all()
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    done <- i - start_i + 1
    cat(sprintf("[%s] %d/%d  %.1f min elapsed, %.1f s/locus, ~%.1f h remaining\n",
                GRP, i, n, el, el * 60 / done, el / done * (n - i) / 60))
    flush.console()
  }
}
save_all(final = TRUE)

panel <- list()
for (pn in names(PAIRS)) {
  ph <- PAIRS[[pn]]; r <- res[[pn]]
  cn_p <- paste0("p_", ph)
  d <- file.path(OUTROOT, pn); dir.create(d, showWarnings = FALSE, recursive = TRUE)

  pro  <- sum(r$processed)
  e1   <- sum(r$processed & !is.na(r[[cn_p[1]]]))
  e2   <- sum(r$processed & !is.na(r[[cn_p[2]]]))
  both <- sum(r$processed & !is.na(r[[cn_p[1]]]) & !is.na(r[[cn_p[2]]]))
  gp   <- sum(r$gate_pass, na.rm = TRUE)
  rep_ <- sum(!is.na(r$rho_hat))
  plim <- sum(r$paramlim_drop, na.rm = TRUE)
  det  <- !is.na(r$rho_lower) & (r$rho_lower > 0 | r$rho_upper < 0)
  w    <- r$rho_upper - r$rho_lower
  bonf <- sum(r$processed & r[[cn_p[1]]] < 0.05 / pro & r[[cn_p[2]]] < 0.05 / pro,
              na.rm = TRUE)
  weaker <- ph[which.min(c(e1, e2))]

  panel[[pn]] <- data.frame(
    pair = pn, traits = paste(ph, collapse = " x "), group = GRP,
    n_snps_common = N_SNPS_COMMON,
    blocks = nrow(r), processed = pro,
    est_1 = e1, est_2 = e2, both_estimable = both,
    gate_pass = gp, reported = rep_,
    pct_reported = round(100 * rep_ / nrow(r), 2),
    paramlim_drops = plim, bonferroni_survivors = bonf,
    ch2_loss_1 = round(100 * mean(r$processed & is.na(r[[cn_p[1]]])), 1),
    ch2_loss_2 = round(100 * mean(r$processed & is.na(r[[cn_p[2]]])), 1),
    weaker_trait = weaker,
    mean_rho = round(mean(r$rho_hat, na.rm = TRUE), 4),
    pct_positive = round(100 * mean(r$rho_hat > 0, na.rm = TRUE), 1),
    n_determinate = sum(det, na.rm = TRUE),
    mean_ci_width = round(mean(w, na.rm = TRUE), 4),
    median_K = median(r$K[r$processed], na.rm = TRUE),
    stringsAsFactors = FALSE)

  sink(file.path(d, "attrition_audit.txt"))
  cat("ATTRITION AUDIT  ", pn, "  ", paste(ph, collapse = " x "), "\n",
      format(Sys.time()), "\n\n", sep = "")
  cat(sprintf("Phenotype group in the pass         : %s\n", paste(PHENOS, collapse = ", ")))
  cat(sprintf("SNPs common to group + reference    : %d\n", N_SNPS_COMMON))
  cat(sprintf("Blocks in partition                 : %d\n", nrow(r)))
  cat(sprintf("CHANNEL 1 processed                 : %d (%.1f%%)\n", pro, 100*pro/nrow(r)))
  cat(sprintf("  %-10s estimable                : %d (%.1f%% of processed)\n",
              ph[1], e1, 100*e1/max(pro,1)))
  cat(sprintf("  %-10s estimable                : %d (%.1f%% of processed)\n",
              ph[2], e2, 100*e2/max(pro,1)))
  cat(sprintf("CHANNEL 2 both estimable            : %d (%.1f%% of processed)\n",
              both, 100*both/max(pro,1)))
  cat(sprintf("CHANNEL 3 both clear P<%.2f          : %d (%.1f%% of both-estimable)\n",
              UNIV_THRESH, gp, 100*gp/max(both,1)))
  cat(sprintf("CHANNEL 4 param.lim drops           : %d\n", plim))
  cat(sprintf("FINAL     reported rho              : %d (%.2f%% of all blocks)\n",
              rep_, 100*rep_/nrow(r)))
  cat(sprintf("\nWeaker trait (sets the ceiling)     : %s\n", weaker))
  cat(sprintf("Bonferroni (0.05/%d) survivors      : %d\n", pro, bonf))
  cat("\nReal K distribution (processed loci):\n"); print(summary(r$K[r$processed]))
  cat("\nNOTE: do NOT test sign against a 50/50 binomial null for this pair.\n")
  cat("Run 14_pair_pipeline.py with this pair's r_e to obtain the filtered null.\n")
  sink()
  cat("wrote", d, "\n")
}

pan <- do.call(rbind, panel)
utils::write.csv(pan, file.path(OUTROOT, sprintf("PANEL_SUMMARY_%s.csv", GRP)),
                 row.names = FALSE)
cat("\n================ PANEL SUMMARY:", GRP, "================\n"); print(pan)
cat("\nMerge the per-group PANEL_SUMMARY_*.csv files to populate Table 8.\n")
if (file.exists(CKPT) && !is.finite(N_LOCI_MAX)) file.remove(CKPT)
