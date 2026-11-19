#!/usr/bin/env Rscript
# 01_harmonize_sumstats.R
suppressMessages(library(data.table))

base <- Sys.getenv("LAVAGATE_BASE")
if (!nzchar(base)) stop("Set LAVAGATE_BASE to the folder that contains SUBMISSION/")
in_asd <- file.path(base, "Section 1_Primary GWAS/Autism Spectrum (Grove et al. 2019)/iPSYCH-PGC_ASD_Nov2017")
in_scz <- file.path(base, "Section 1_Primary GWAS/Schizophrenia (Trubetskoy et al. 2022)/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv")
outdir <- file.path(base, "Section 7_LAVA Analysis/sumstats")

ASD_CASES <- 18381L; ASD_CONTROLS <- 27969L; ASD_N <- ASD_CASES + ASD_CONTROLS

cat("[ASD] reading", in_asd, "\n")
asd <- fread(in_asd, data.table = FALSE)
cat("[ASD] raw rows:", nrow(asd), " cols:", paste(names(asd), collapse=","), "\n")
stopifnot(all(c("SNP","A1","A2","OR","P") %in% names(asd)))
asd_out <- data.frame(
  SNP = asd$SNP, A1 = asd$A1, A2 = asd$A2,
  OR  = asd$OR,  P  = asd$P,  N  = ASD_N
)
ok <- complete.cases(asd_out) & asd_out$OR > 0 & asd_out$P > 0 & asd_out$P <= 1
cat("[ASD] dropping", sum(!ok), "rows failing QC\n")
asd_out <- asd_out[ok, ]
f_asd <- file.path(outdir, "ASD_Grove2019.lava.sumstats.gz")
fwrite(asd_out, f_asd, sep = "\t")
cat("[ASD] wrote", nrow(asd_out), "rows ->", f_asd, "\n\n")

cat("[SCZ] reading", in_scz, "\n")
scz <- fread(in_scz, skip = "CHROM", data.table = FALSE)
cat("[SCZ] raw rows:", nrow(scz), " cols:", paste(names(scz), collapse=","), "\n")
stopifnot(all(c("ID","A1","A2","BETA","PVAL","NCAS","NCON") %in% names(scz)))
scz_out <- data.frame(
  SNP  = scz$ID, A1 = scz$A1, A2 = scz$A2,
  BETA = scz$BETA, P = scz$PVAL,
  N    = scz$NCAS + scz$NCON
)
ok <- complete.cases(scz_out) & is.finite(scz_out$BETA) & scz_out$P > 0 & scz_out$P <= 1 & scz_out$N > 0
cat("[SCZ] dropping", sum(!ok), "rows failing QC\n")
scz_out <- scz_out[ok, ]
f_scz <- file.path(outdir, "SCZ_PGC3_EUR.lava.sumstats.gz")
fwrite(scz_out, f_scz, sep = "\t")
cat("[SCZ] wrote", nrow(scz_out), "rows ->", f_scz, "\n\n")

scz_cases    <- as.integer(median(scz$NCAS, na.rm = TRUE))
scz_controls <- as.integer(median(scz$NCON, na.rm = TRUE))
cat("[SCZ] nominal cases/controls (median per-SNP):", scz_cases, "/", scz_controls, "\n")

info <- data.frame(
  phenotype = c("ASD", "SCZ"),
  cases     = c(ASD_CASES,   scz_cases),
  controls  = c(ASD_CONTROLS, scz_controls),
  filename  = c(basename(f_asd), basename(f_scz))
)
f_info <- file.path(base, "Section 7_LAVA Analysis/input.info.txt")
write.table(info, f_info, row.names = FALSE, quote = FALSE, sep = "\t")
cat("[INFO] wrote input info ->", f_info, "\n")
print(info)
cat("\nDONE.\n")
