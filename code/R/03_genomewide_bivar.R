#!/usr/bin/env Rscript
# 03_genomewide_bivar.R -- Genome-wide ASD x SCZ local genetic correlation
suppressMessages(library(LAVA))

base <- Sys.getenv("LAVAGATE_BASE")
if (!nzchar(base)) stop("Set LAVAGATE_BASE to the folder that contains SUBMISSION/")
adir    <- file.path(base, "Section 7_LAVA Analysis")
info    <- file.path(adir, "input.info.txt")
locfile <- file.path(base, "LAVA/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
refpref <- file.path(base, "Section 2_LD Reference/ukb_ref_consolidated/lava-ukb-v1.1")
sumdir  <- file.path(adir, "sumstats")
outdir  <- file.path(adir, "results")
logf    <- file.path(outdir, "genomewide_progress.log")
ovl     <- file.path(adir, "sample.overlap.txt")

UNIV_THR <- 0.05
sample.overlap.file <- if (file.exists(ovl)) ovl else NULL

logmsg <- function(...) { cat(sprintf(...), "\n", file = logf, append = TRUE); cat(sprintf(...), "\n") }
cat("", file = logf)
if (is.null(sample.overlap.file))
  logmsg("*** WARNING: no sample.overlap.txt found -> running with NULL overlap. Estimates NOT corrected for sample overlap. ***")

rds <- file.path(adir, "input_object.rds")
if (file.exists(rds)) {
  logmsg("[%s] loading cached input object ...", format(Sys.time(), "%H:%M:%S"))
  input <- readRDS(rds)
} else {
  logmsg("[%s] process.input ...", format(Sys.time(), "%H:%M:%S"))
  input <- process.input(input.info.file = info, sample.overlap.file = sample.overlap.file,
                         ref.prefix = refpref, input.dir = sumdir)
  saveRDS(input, rds)
}
logmsg("[%s] input ready: %d analysis SNPs", format(Sys.time(), "%H:%M:%S"), length(input$analysis.snps))

loci <- read.loci(locfile)
NL <- nrow(loci)
logmsg("[%s] %d loci to process", format(Sys.time(), "%H:%M:%S"), NL)

univ.res <- list(); bivar.res <- list(); n.ok <- 0L; n.biv <- 0L
for (i in seq_len(NL)) {
  locus <- tryCatch(process.locus(loci[i, ], input), error = function(e) NULL)
  if (!is.null(locus)) {
    n.ok <- n.ok + 1L
    ub <- tryCatch(run.univ.bivar(locus, univ.thresh = UNIV_THR), error = function(e) NULL)
    if (!is.null(ub)) {
      if (!is.null(ub$univ))  univ.res[[length(univ.res)+1]] <- data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, ub$univ)
      if (!is.null(ub$bivar) && nrow(ub$bivar) > 0) { bivar.res[[length(bivar.res)+1]] <- data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, ub$bivar); n.biv <- n.biv + nrow(ub$bivar) }
    }
  }
  if (i %% 100 == 0 || i == NL) {
    logmsg("[%s] %d/%d loci | analysable=%d | bivariate=%d", format(Sys.time(), "%H:%M:%S"), i, NL, n.ok, n.biv)
    if (length(univ.res))  write.table(do.call(rbind, univ.res),  file.path(outdir, "genomewide_univ.txt"),  row.names=FALSE, quote=FALSE, sep="\t")
    if (length(bivar.res)) write.table(do.call(rbind, bivar.res), file.path(outdir, "genomewide_bivar.txt"), row.names=FALSE, quote=FALSE, sep="\t")
  }
}
logmsg("[%s] DONE. analysable loci=%d, bivariate estimates=%d", format(Sys.time(), "%H:%M:%S"), n.ok, n.biv)
