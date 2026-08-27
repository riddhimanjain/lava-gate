#!/usr/bin/env Rscript
# 04_sample_overlap.R -- Estimate ASD/SCZ sample overlap via cross-trait LDSC
suppressMessages({library(GenomicSEM); library(data.table)})

base <- Sys.getenv("LAVAGATE_BASE")
if (!nzchar(base)) stop("Set LAVAGATE_BASE to the folder that contains SUBMISSION/")
adir   <- file.path(base, "Section 7_LAVA Analysis")
sumdir <- file.path(adir, "sumstats")
ldscdir<- file.path(base, "Section 3_Downstream Genetic Architecture/LDSC_reference (Zenodo 10515792)")
hm3    <- file.path(ldscdir, "w_hm3.snplist.txt")
ld     <- file.path(ldscdir, "ld_clean")
wld    <- file.path(ldscdir, "wld_clean")
mungedir <- file.path(adir, "munged")
dir.create(mungedir, showWarnings = FALSE)
setwd(mungedir)

f_asd <- file.path(sumdir, "ASD_Grove2019.lava.sumstats.gz")
f_scz <- file.path(sumdir, "SCZ_PGC3_EUR.lava.sumstats.gz")

cat("=== munge() ASD (effect=OR) ===\n")
munge(files = f_asd, hm3 = hm3, trait.names = "ASD",
      column.names = list(SNP="SNP", A1="A1", A2="A2", effect="OR", P="P", N="N"),
      info.filter = 0.9, maf.filter = 0.01)

cat("=== munge() SCZ (effect=BETA) ===\n")
munge(files = f_scz, hm3 = hm3, trait.names = "SCZ",
      column.names = list(SNP="SNP", A1="A1", A2="A2", effect="BETA", P="P", N="N"),
      info.filter = 0.9, maf.filter = 0.01)

cat("=== ldsc() ===\n")
ASD_CASES <- 18381; ASD_CONTROLS <- 27969
SCZ_CASES <- 53386;  SCZ_CONTROLS <- 77258
sample.prev     <- c(ASD_CASES/(ASD_CASES+ASD_CONTROLS), SCZ_CASES/(SCZ_CASES+SCZ_CONTROLS))
population.prev <- c(0.015, 0.01)

traits <- c(file.path(mungedir, "ASD.sumstats.gz"), file.path(mungedir, "SCZ.sumstats.gz"))
LDSCoutput <- ldsc(traits = traits, sample.prev = sample.prev, population.prev = population.prev,
                   ld = ld, wld = wld, trait.names = c("ASD","SCZ"), sep_weights = TRUE)

cat("\n--- Genetic covariance (S) ---\n"); print(LDSCoutput$S)
cat("\n--- Intercept matrix (I) ---\n");   print(LDSCoutput$I)

I <- LDSCoutput$I
rownames(I) <- colnames(I) <- c("ASD", "SCZ")
mat <- round(cov2cor(I), 5)
cat("\n--- Standardised sampling-correlation matrix ---\n"); print(mat)

f_out <- file.path(adir, "sample.overlap.txt")
write.table(mat, f_out, quote = FALSE)
cat("\nWrote sample overlap matrix ->", f_out, "\n")
cat("\nDONE.\n")
