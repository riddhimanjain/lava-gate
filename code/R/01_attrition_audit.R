# 01_attrition_audit.R  --  Quantify every attrition channel on real data (Task 2)

suppressMessages(library(LAVA))

base <- Sys.getenv("LAVAGATE_BASE")
if (!nzchar(base)) stop("Set LAVAGATE_BASE to the folder that contains SUBMISSION/")
sumdir  <- file.path(base, "Section 7_LAVA Analysis/sumstats")
refpref <- file.path(base, "Section 2_LD Reference/ukb_ref_consolidated/lava-ukb-v1.1")
locfile <- file.path(base, "LAVA/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
outdir  <- file.path(base, "analysis/out_R")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

infofile    <- file.path(base, "Section 7_LAVA Analysis/input.info.txt")
overlapfile <- file.path(base, "Section 7_LAVA Analysis/sample.overlap.txt")

stopifnot(file.exists(infofile), file.exists(locfile))

cat("Reading input...\n")
input <- process.input(input.info.file  = infofile,
                       sample.overlap.file = if (file.exists(overlapfile)) overlapfile else NULL,
                       ref.prefix = refpref,
                       phenos = c("ASD", "SCZ"),
                       input.dir = sumdir)

loci <- read.loci(locfile)
n_loci <- nrow(loci)
cat("Loci in partition:", n_loci, "\n")

if (!is.null(input$sample.overlap)) {
  cat("\nSample overlap matrix in use:\n"); print(input$sample.overlap)
  write.csv(as.data.frame(input$sample.overlap),
            file.path(outdir, "sample_overlap_actual.csv"))
}

ckpt <- file.path(outdir, "attrition_audit_checkpoint.rds")
CKPT_EVERY <- 25

res <- vector("list", n_loci)
start_i <- 1
if (file.exists(ckpt)) {
  saved <- readRDS(ckpt)
  res <- saved$res
  start_i <- saved$next_i
  cat("Resuming from checkpoint: locus", start_i, "of", n_loci,
      "(", round(100*(start_i-1)/n_loci,1), "% already done)\n")
}

t_start <- Sys.time()
pb  <- txtProgressBar(min = 0, max = n_loci, style = 3, initial = start_i - 1)

for (i in seq(start_i, n_loci)) {
  setTxtProgressBar(pb, i)
  row <- list(locus = loci$LOC[i], chr = loci$CHR[i],
              start = loci$START[i], stop = loci$STOP[i],
              processed = FALSE, K = NA_integer_, n_snps = NA_integer_,
              h2_ASD = NA_real_, h2_SCZ = NA_real_,
              p_ASD = NA_real_, p_SCZ = NA_real_,
              gate_pass = NA, rho_hat = NA_real_, rho_raw = NA_real_,
              paramlim_drop = NA, rho_lower = NA_real_, rho_upper = NA_real_,
              bivar_p = NA_real_)

  locus <- tryCatch(process.locus(loci[i, ], input), error = function(e) NULL)

  if (!is.null(locus)) {
    row$processed <- TRUE
    row$K      <- locus$K
    row$n_snps <- if (!is.null(locus$n.snps)) locus$n.snps else NA_integer_

    u <- tryCatch(run.univ(locus), error = function(e) NULL)
    if (is.null(u)) {
      row$gate_pass <- FALSE
    } else {
      h2_asd <- u$h2.obs[u$phen == "ASD"]; h2_scz <- u$h2.obs[u$phen == "SCZ"]
      p_asd  <- u$p[u$phen == "ASD"];      p_scz  <- u$p[u$phen == "SCZ"]
      row$h2_ASD <- if (length(h2_asd)) h2_asd[1] else NA_real_
      row$h2_SCZ <- if (length(h2_scz)) h2_scz[1] else NA_real_
      row$p_ASD  <- if (length(p_asd))  p_asd[1]  else NA_real_
      row$p_SCZ  <- if (length(p_scz))  p_scz[1]  else NA_real_
      row$gate_pass <- length(p_asd) > 0 && length(p_scz) > 0 &&
                        p_asd[1] < 0.05 && p_scz[1] < 0.05
    }

    if (row$gate_pass) {
      b_raw <- tryCatch(run.bivar(locus, param.lim = Inf, CIs = FALSE,
                                  p.values = FALSE, cap.estimates = FALSE),
                        error = function(e) NULL)
      if (!is.null(b_raw) && nrow(b_raw) > 0) row$rho_raw <- b_raw$rho[1]

      b <- tryCatch(run.bivar(locus), error = function(e) NULL)
      if (!is.null(b) && nrow(b) > 0) {
        row$rho_hat   <- b$rho[1]
        row$rho_lower <- b$rho.lower[1]
        row$rho_upper <- b$rho.upper[1]
        row$bivar_p   <- b$p[1]
        row$paramlim_drop <- is.na(b$rho[1]) && !is.na(row$rho_raw)
      }
    }
  }
  res[[i]] <- as.data.frame(row, stringsAsFactors = FALSE)

  if (i %% CKPT_EVERY == 0 || i == n_loci) {
    saveRDS(list(res = res, next_i = i + 1), ckpt)
    elapsed <- as.numeric(Sys.time() - t_start, units = "secs")
    done_now <- i - start_i + 1
    rate <- if (done_now > 0) elapsed / done_now else NA_real_
    eta_min <- if (!is.na(rate)) round(rate * (n_loci - i) / 60, 1) else NA
    cat(sprintf("\n[checkpoint] locus %d/%d done, %.1f sec/locus, ~%s min remaining\n",
                i, n_loci, rate, format(eta_min)))
    flush.console()
  }
}
close(pb)

df <- do.call(rbind, res)
write.csv(df, file.path(outdir, "attrition_audit.csv"), row.names = FALSE)
if (file.exists(ckpt)) file.remove(ckpt)

n_proc     <- sum(df$processed)
n_gate     <- sum(df$gate_pass, na.rm = TRUE)
n_reported <- sum(!is.na(df$rho_hat))
n_plim     <- sum(df$paramlim_drop, na.rm = TRUE)
asd_pass   <- mean(df$p_ASD[df$processed] < 0.05, na.rm = TRUE)
scz_pass   <- mean(df$p_SCZ[df$processed] < 0.05, na.rm = TRUE)

sink(file.path(outdir, "attrition_audit.txt"))
cat("ATTRITION AUDIT  ASD x SCZ\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")
cat("Blocks in partition                :", n_loci, "\n")
cat("CHANNEL 1  processed by process.locus:", n_proc,
    sprintf("(%.1f%%)  -- lost upstream: %d\n", 100*n_proc/n_loci, n_loci - n_proc))
cat("CHANNEL 2  passed univariate gate    :", n_gate,
    sprintf("(%.1f%% of processed)\n", 100*n_gate/max(n_proc,1)))
cat("           ASD alone clears P<0.05   :", sprintf("%.1f%%\n", 100*asd_pass))
cat("           SCZ alone clears P<0.05   :", sprintf("%.1f%%\n", 100*scz_pass))
cat("CHANNEL 3  dropped by param.lim      :", n_plim, "\n")
cat("FINAL      loci with a reported rho  :", n_reported,
    sprintf("(%.2f%% of all blocks)\n", 100*n_reported/n_loci))
cat("\nWhich trait sets the ceiling? ",
    ifelse(asd_pass < scz_pass, "ASD", "SCZ"), "\n")
cat("\nReal K distribution (processed loci):\n")
print(summary(df$K[df$processed]))
cat("\nBonferroni check: loci surviving univariate gate at 0.05/",
    n_proc, " = ", signif(0.05/max(n_proc,1), 3), "\n", sep = "")
cat("  ", sum(df$p_ASD < 0.05/n_proc & df$p_SCZ < 0.05/n_proc, na.rm = TRUE), "\n")
sink()

cat("\nDone. See", file.path(outdir, "attrition_audit.txt"), "\n")
cat("\nIMPORTANT: the K distribution above feeds Tier B. Report it -- the Tier A\n")
cat("sweep assumed K values that may not match reality.\n")
