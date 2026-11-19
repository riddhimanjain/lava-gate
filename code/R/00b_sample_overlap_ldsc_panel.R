#!/usr/bin/env Rscript
# 00b_sample_overlap_ldsc_panel.R  --  cross-trait LDSC for the whole panel
suppressMessages({ library(GenomicSEM); library(data.table) })

BASE <- Sys.getenv("LAVAGATE_BASE")
if (!nzchar(BASE)) stop("Set LAVAGATE_BASE to the folder that contains SUBMISSION/")
ADIR     <- file.path(BASE, "Section 7_LAVA Analysis")
SUMDIR   <- file.path(ADIR, "sumstats")
MUNGEDIR <- file.path(ADIR, "munged")
LDBASE <- Sys.getenv("LAVAGATE_LDSC_REF")
if (!nzchar(LDBASE)) stop("Set LAVAGATE_LDSC_REF to the LDSC reference folder")
LD       <- file.path(LDBASE, "ld_clean")
WLD      <- file.path(LDBASE, "wld_clean")
HM3      <- file.path(LDBASE, "w_hm3.snplist.txt")

PUBLISHED_RE_ASD_SCZ <- 0.02756
RE_TOL <- 0.002

stopifnot(dir.exists(LD), dir.exists(WLD), file.exists(HM3))

POP_PREV <- c(ASD = 0.015, SCZ = 0.010, BIP = 0.010, BIP_noUKB = 0.010,
              MDD = 0.150, MDD_noUKB = 0.150, ADHD = 0.050, AN = 0.010)

info <- read.table(file.path(ADIR, "multipair.info.txt"), header = TRUE,
                   stringsAsFactors = FALSE)
stopifnot("run 01b_harmonize_panel.R first" = nrow(info) >= 8)
info$N <- info$cases + info$controls
info$sample.prev <- info$cases / info$N
stopifnot("a phenotype has no population prevalence entry" =
            all(info$phenotype %in% names(POP_PREV)))

dir.create(MUNGEDIR, showWarnings = FALSE, recursive = TRUE)
todo <- info[!file.exists(file.path(MUNGEDIR, paste0(info$phenotype, ".sumstats.gz"))), ]
if (nrow(todo)) {
  cat("Munging", nrow(todo), "file(s):", paste(todo$phenotype, collapse = ", "), "\n")
  wd <- setwd(MUNGEDIR); on.exit(setwd(wd), add = TRUE)
  for (i in seq_len(nrow(todo))) {
    cat(sprintf("  [%d/%d] %s\n", i, nrow(todo), todo$phenotype[i]))
    munge(files       = file.path(SUMDIR, todo$filename[i]),
          hm3         = HM3,
          trait.names = todo$phenotype[i],
          N           = NA_real_)
    gc()
  }
  setwd(wd)
} else cat("All phenotypes already munged.\n")

munged <- file.path(MUNGEDIR, paste0(info$phenotype, ".sumstats.gz"))
stopifnot("munging did not produce every file" = all(file.exists(munged)))

cat("\n=== cross-trait LDSC over", nrow(info), "phenotypes ===\n")
print(info[, c("phenotype", "cases", "controls", "N", "sample.prev")])

out <- ldsc(traits          = munged,
            sample.prev     = info$sample.prev,
            population.prev = unname(POP_PREV[info$phenotype]),
            ld = LD, wld = WLD,
            trait.names     = info$phenotype,
            sep_weights     = TRUE)

I <- out$I; S <- out$S
dimnames(I) <- list(info$phenotype, info$phenotype)
dimnames(S) <- list(info$phenotype, info$phenotype)

cat("\n--- Intercept matrix I ---\n"); print(round(I, 5))
cat("\n--- Genetic covariance matrix S ---\n"); print(signif(S, 4))

RE <- round(cov2cor(I), 5)
dimnames(RE) <- list(info$phenotype, info$phenotype)
cat("\n--- Sampling-correlation matrix cov2cor(I) : this is r_e ---\n")
print(RE)

naive <- outer(info$N, info$N, function(a, b) 1 / sqrt(a * b)) * I
cat("\n--- 'intercept / sqrt(N1*N2)' for comparison (NOT used; not a correlation) ---\n")
print(signif(naive, 3))

re_p4 <- RE["ASD", "SCZ"]
cat(sprintf("\nASD x SCZ r_e : panel run %.5f vs published %.5f (diff %.5f)\n",
            re_p4, PUBLISHED_RE_ASD_SCZ, re_p4 - PUBLISHED_RE_ASD_SCZ))
if (abs(re_p4 - PUBLISHED_RE_ASD_SCZ) > RE_TOL) {
  cat("\n*** WARNING ***\n")
  cat("The panel run does not reproduce the published ASD x SCZ overlap. LDSC\n")
  cat("intercepts are estimated pairwise, so this cell should be identical to\n")
  cat("the two-trait run. A discrepancy means something changed in the inputs\n")
  cat("(munging, reference panel, or the harmonised sumstats) and the panel\n")
  cat("MUST NOT be used until it is explained.\n")
} else {
  cat("Regression test PASSES: the panel reproduces the published P4 overlap.\n")
}

f_overlap <- file.path(ADIR, "multipair.sample.overlap.txt")
write.table(RE, f_overlap, quote = FALSE)

PAIRS <- list(P1 = c("SCZ","BIP"), P2a = c("MDD","BIP"),
              P2b = c("MDD_noUKB","BIP_noUKB"), P3 = c("ASD","ADHD"),
              P4 = c("ASD","SCZ"), P5 = c("AN","ASD"))
pair_re <- do.call(rbind, lapply(names(PAIRS), function(pn) {
  p <- PAIRS[[pn]]
  data.frame(pair = pn, trait1 = p[1], trait2 = p[2],
             r_e = RE[p[1], p[2]], ldsc_intercept = I[p[1], p[2]],
             predicted_gate_intercept = 0.68 * RE[p[1], p[2]],
             rg = S[p[1], p[2]] / sqrt(S[p[1], p[1]] * S[p[2], p[2]]),
             stringsAsFactors = FALSE)
}))
f_pairs <- file.path(ADIR, "multipair.pair_re.csv")
write.csv(pair_re, f_pairs, row.names = FALSE)

cat("\n================ PER-PAIR SAMPLE OVERLAP ================\n")
print(pair_re, row.names = FALSE)
cat(sprintf("\nP2 CONTRAST: r_e falls from %.5f (P2a, with UKBB) to %.5f (P2b, without).\n",
            RE["MDD","BIP"], RE["MDD_noUKB","BIP_noUKB"]))
cat(sprintf("Predicted gate intercept therefore moves %.4f -> %.4f.\n",
            0.68*RE["MDD","BIP"], 0.68*RE["MDD_noUKB","BIP_noUKB"]))
cat("Same two disorders, same LD reference, same blocks. R1 measures whether the\n")
cat("fitted intercept actually moves as predicted.\n")

saveRDS(out, file.path(ADIR, "multipair.ldsc.rds"))
cat("\nwrote", f_overlap, "\n     ", f_pairs, "\n")
cat("\nDONE.\n")
