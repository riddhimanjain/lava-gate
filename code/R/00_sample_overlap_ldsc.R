#!/usr/bin/env Rscript
# 00_sample_overlap_ldsc.R  --  cross-trait LDSC -> sample overlap r_e   (Task 1)
suppressMessages({library(GenomicSEM)})

proj <- Sys.getenv("LAVAGATE_BASE")
if (!nzchar(proj)) stop("Set LAVAGATE_BASE to the folder that contains SUBMISSION/")
adir    <- file.path(proj, "Section 7_LAVA Analysis")
mungedir<- file.path(adir, "munged")
ldbase <- Sys.getenv("LAVAGATE_LDSC_REF")
if (!nzchar(ldbase)) stop("Set LAVAGATE_LDSC_REF to the LDSC reference folder")
ld      <- file.path(ldbase, "ld_clean")
wld     <- file.path(ldbase, "wld_clean")

f_asd <- file.path(mungedir, "ASD.sumstats.gz")
f_scz <- file.path(mungedir, "SCZ.sumstats.gz")
stopifnot(file.exists(f_asd), file.exists(f_scz), dir.exists(ld), dir.exists(wld))

ASD_CASES <- 18381; ASD_CONTROLS <- 27969; N1 <- ASD_CASES + ASD_CONTROLS
SCZ_CASES <- 53386; SCZ_CONTROLS <- 77258; N2 <- SCZ_CASES + SCZ_CONTROLS
sample.prev     <- c(ASD_CASES / N1, SCZ_CASES / N2)
population.prev <- c(0.015, 0.01)

cat("=== cross-trait ldsc(): ASD x SCZ ===\n")
cat("ASD file:", f_asd, "\n")
cat("SCZ file:", f_scz, "\n")
cat("N1 (ASD) =", N1, " N2 (SCZ) =", N2, "\n")
cat("ld  =", ld, "\n")
cat("wld =", wld, "\n\n")

LDSCoutput <- ldsc(traits = c(f_asd, f_scz),
                    sample.prev = sample.prev, population.prev = population.prev,
                    ld = ld, wld = wld, trait.names = c("ASD", "SCZ"), sep_weights = TRUE)

cat("\n--- Genetic covariance / heritability matrix (S) ---\n"); print(LDSCoutput$S)
cat("\n--- Intercept matrix (I) ---\n"); print(LDSCoutput$I)

I <- LDSCoutput$I
rownames(I) <- colnames(I) <- c("ASD", "SCZ")
gcov_int <- I["ASD", "SCZ"]
overlap_term <- gcov_int / sqrt(N1 * N2)
cat("\ngcov_int (cross-trait LDSC intercept) =", gcov_int, "\n")
cat("intercept / sqrt(N1*N2) =", overlap_term, "\n")

mat <- round(cov2cor(I), 5)
rownames(mat) <- colnames(mat) <- c("ASD", "SCZ")
cat("\n--- Standardised sampling-correlation matrix (cov2cor(I), 1s on diagonal) ---\n")
print(mat)

f_out <- file.path(adir, "sample.overlap.txt")
write.table(mat, f_out, quote = FALSE)
cat("\nWrote sample overlap matrix ->", f_out, "\n")
cat("\nDONE.\n")
