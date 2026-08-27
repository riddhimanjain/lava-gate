#!/usr/bin/env Rscript
# 01b_harmonize_panel.R
suppressMessages(library(data.table))

BASE <- Sys.getenv("LAVAGATE_BASE")
if (!nzchar(BASE)) stop("Set LAVAGATE_BASE to the folder that contains SUBMISSION/")
GWAS   <- file.path(BASE, "Section 1_Primary GWAS")
ADIR   <- file.path(BASE, "Section 7_LAVA Analysis")
OUTDIR <- file.path(ADIR, "sumstats")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

PANEL <- list(
  list(id = "BIP",       ukbb = TRUE,  a1 = "A1",  a2 = "A2",
       ev = 'VCF header: BETA is "beta or ln(OR) of A1"',
       src = file.path(GWAS, "Bipolar (Mullins et al. 2021)/pgc-bip2021-all.vcf.tsv.gz")),
  list(id = "BIP_noUKB", ukbb = FALSE, a1 = "A1",  a2 = "A2",
       ev = "Ricopili daner: OR is the odds ratio for A1",
       src = file.path(GWAS, "Bipolar noUKBB (Mullins et al. 2021)/daner_bip_pgc3_nm_noukbiobank.gz")),
  list(id = "MDD",       ukbb = TRUE,  a1 = "A1",  a2 = "A2",
       ev = "Ricopili daner: OR is the odds ratio for A1",
       src = file.path(GWAS, "Major Depression (Wray et al. 2018)/MDD2018_ex23andMe.gz")),
  list(id = "MDD_noUKB", ukbb = FALSE, a1 = "A1",  a2 = "A2",
       ev = "Ricopili daner: OR is the odds ratio for A1",
       src = file.path(GWAS, "Major Depression rmUKBB (Wray et al. 2018)/daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz")),
  list(id = "ADHD",      ukbb = FALSE, a1 = "A1",  a2 = "A2",
       ev = "Ricopili daner: OR is the odds ratio for A1",
       src = file.path(GWAS, "ADHD (Demontis et al. 2023)/ADHD2022_iPSYCH_deCODE_PGC.meta.gz")),
  list(id = "AN",        ukbb = FALSE, a1 = "ALT", a2 = "REF",
       ev = 'VCF header: BETA is "beta or ln(OR) of ALT"',
       src = file.path(GWAS, "Anorexia Nervosa (Watson et al. 2019)/pgcAN2.2019-07.vcf.tsv.gz"))
)

CAND <- list(
  SNP    = c("SNP", "ID", "RSID", "MARKERNAME", "SNPID"),
  A1     = c("A1", "EFFECT_ALLELE", "ALT"),
  A2     = c("A2", "OTHER_ALLELE", "REF"),
  OR     = c("OR"),
  BETA   = c("BETA", "LOGOR", "LOG_ODDS", "EFFECT"),
  P      = c("P", "PVAL", "PVALUE", "P_VALUE", "P.VALUE"),
  NCAS   = c("NCA", "NCAS", "N_CAS", "NCASE", "NCASES"),
  NCON   = c("NCO", "NCON", "N_CON", "NCONTROL", "NCONTROLS"),
  N      = c("N", "NTOT", "TOTALN", "N_TOTAL")
)

pick <- function(nms, key) {
  hit <- CAND[[key]][match(toupper(CAND[[key]]), toupper(nms), nomatch = 0L) > 0L]
  if (!length(hit)) return(NA_character_)
  nms[match(toupper(hit[1]), toupper(nms))]
}

read_any <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  n_meta <- 0L
  repeat {
    l <- readLines(con, 1L)
    if (!length(l) || !startsWith(l, "##")) break
    n_meta <- n_meta + 1L
  }
  close(con); on.exit()

  con <- gzfile(path, "rt")
  hdr <- readLines(con, n_meta + 1L)[n_meta + 1L]
  close(con)
  n_tab <- length(strsplit(hdr, "\t", fixed = TRUE)[[1]])
  n_ws  <- length(strsplit(trimws(hdr), "[[:space:]]+")[[1]])
  if (n_tab > 1) { fs <- "-F'\\t'"; n_field <- n_tab }
  else           { fs <- "";        n_field <- n_ws  }
  stopifnot("cannot determine the header's field count" = n_field > 3)
  cat(sprintf("  layout   : %s-delimited, %d fields\n",
              if (n_tab > 1) "tab" else "whitespace", n_field))

  tmp <- tempfile(fileext = ".tsv")
  rc <- system2("bash", c("-c", shQuote(sprintf(
          "zcat %s | awk %s 'NF==%d' > %s",
          shQuote(path), fs, n_field, shQuote(tmp)))))
  stopifnot("bash/zcat filter failed" = rc == 0, file.exists(tmp))
  on.exit(unlink(tmp), add = TRUE)
  d <- fread(file = tmp, header = TRUE, data.table = FALSE)
  names(d)[1] <- sub("^#", "", names(d)[1])

  total <- as.integer(system2("bash", c("-c", shQuote(
      sprintf("zcat %s | wc -l", shQuote(path)))), stdout = TRUE))
  dropped <- total - n_meta - 1L - nrow(d)
  if (!is.na(dropped) && dropped > 0) {
    cat(sprintf("  MALFORMED: dropped %d line(s) whose field count != %d (of %d data lines)\n",
                dropped, n_field, total - n_meta - 1L))
    stopifnot("too many malformed lines to be junk; inspect the file by hand" =
                dropped < 0.001 * total)
  }
  return(d)
}

harmonise <- function(spec) {
  cat("\n=====", spec$id, "=====\n", spec$src, "\n")
  stopifnot("source file missing" = file.exists(spec$src))
  d <- read_any(spec$src)
  nms <- names(d)
  cat("  raw rows:", nrow(d), "\n  columns :", paste(nms, collapse = ", "), "\n")

  c_snp <- pick(nms, "SNP")
  c_a1 <- spec$a1; c_a2 <- spec$a2
  stopifnot("declared effect-allele column absent from file"  = c_a1 %in% nms,
            "declared other-allele column absent from file"   = c_a2 %in% nms)
  cat("  A1 (effect):", c_a1, " A2:", c_a2, " evidence:", spec$ev, "\n")
  c_or  <- pick(nms, "OR");  c_b  <- pick(nms, "BETA"); c_p <- pick(nms, "P")
  c_nca <- pick(nms, "NCAS"); c_nco <- pick(nms, "NCON"); c_n <- pick(nms, "N")
  stopifnot("cannot resolve SNP/A1/A2/P" =
              !any(is.na(c(c_snp, c_a1, c_a2, c_p))))
  stopifnot("no effect column (neither OR nor BETA)" = !(is.na(c_or) && is.na(c_b)))

  frq_a <- grep("^FRQ_A_[0-9]+$", nms, value = TRUE)
  frq_u <- grep("^FRQ_U_[0-9]+$", nms, value = TRUE)
  if (!is.na(c_nca) && !is.na(c_nco) &&
      anyNA(d[[c_nca]]) && length(frq_a) == 1 && length(frq_u) == 1) {
    nom_ca <- as.integer(sub("^FRQ_A_", "", frq_a))
    nom_co <- as.integer(sub("^FRQ_U_", "", frq_u))
    n_miss <- sum(is.na(d[[c_nca]]))
    cat(sprintf(paste0("  SPLIT-FILE N: %d of %d rows lack per-SNP N. Using the",
                       " nominal %d + %d = %d for EVERY SNP of this trait.\n"),
                n_miss, nrow(d), nom_ca, nom_co, nom_ca + nom_co))
    obs <- d[[c_nca]] + d[[c_nco]]
    cat(sprintf("  (where per-SNP N is present its median is %.0f, %.1f%% of nominal)\n",
                median(obs, na.rm = TRUE),
                100 * median(obs, na.rm = TRUE) / (nom_ca + nom_co)))
    d[[c_nca]] <- nom_ca; d[[c_nco]] <- nom_co
  }

  if (!is.na(c_nca) && !is.na(c_nco)) {
    N <- as.numeric(d[[c_nca]]) + as.numeric(d[[c_nco]])
    n_src <- sprintf("per-SNP %s + %s", c_nca, c_nco)
    cases    <- as.integer(median(d[[c_nca]], na.rm = TRUE))
    controls <- as.integer(median(d[[c_nco]], na.rm = TRUE))
    cases_max    <- as.integer(max(d[[c_nca]], na.rm = TRUE))
    controls_max <- as.integer(max(d[[c_nco]], na.rm = TRUE))
  } else if (!is.na(c_n)) {
    N <- as.numeric(d[[c_n]]); n_src <- sprintf("per-SNP %s", c_n)
    cases <- controls <- cases_max <- controls_max <- NA_integer_
  } else {
    stop(sprintf("%s: no per-SNP sample size columns; refusing to invent one", spec$id))
  }

  out <- data.frame(SNP = d[[c_snp]], A1 = d[[c_a1]], A2 = d[[c_a2]],
                    P = as.numeric(d[[c_p]]), N = N, stringsAsFactors = FALSE)
  if (!is.na(c_or)) { out$OR <- as.numeric(d[[c_or]]); eff <- paste0("OR=", c_or) }
  else              { out$BETA <- as.numeric(d[[c_b]]); eff <- paste0("BETA=", c_b) }

  cat(sprintf("  resolved : SNP=%s A1=%s A2=%s %s P=%s N=%s\n",
              c_snp, c_a1, c_a2, eff, c_p, n_src))

  ok <- complete.cases(out) & out$P > 0 & out$P <= 1 & out$N > 0
  if (!is.na(c_or)) ok <- ok & is.finite(out$OR) & out$OR > 0
  else              ok <- ok & is.finite(out$BETA)
  ok <- ok & out$A1 %in% c("A","C","G","T","a","c","g","t") &
             out$A2 %in% c("A","C","G","T","a","c","g","t")
  dup <- duplicated(out$SNP)
  cat(sprintf("  QC       : dropping %d failing rows, %d duplicate SNP IDs\n",
              sum(!ok), sum(ok & dup)))
  out <- out[ok & !dup, ]

  f <- file.path(OUTDIR, sprintf("%s.lava.sumstats.gz", spec$id))
  fwrite(out, f, sep = "\t")
  cat(sprintf("  wrote    : %d rows -> %s\n", nrow(out), basename(f)))
  cat(sprintf("  N        : median %.0f, range %.0f - %.0f\n",
              median(out$N), min(out$N), max(out$N)))
  if (!is.na(cases))
    cat(sprintf("  cases/controls: median %d / %d   (max %d / %d)\n",
                cases, controls, cases_max, controls_max))

  data.frame(phenotype = spec$id, cases = cases, controls = controls,
             cases_max = cases_max, controls_max = controls_max,
             includes_ukbb = spec$ukbb, n_snps = nrow(out),
             effect_allele = c_a1, effect_allele_evidence = spec$ev,
             filename = basename(f), source = basename(spec$src),
             stringsAsFactors = FALSE)
}

rows <- do.call(rbind, lapply(PANEL, harmonise))

old <- read.table(file.path(ADIR, "input.info.txt"), header = TRUE,
                  stringsAsFactors = FALSE)
stopifnot("existing input.info.txt must still hold ASD and SCZ" =
            all(c("ASD", "SCZ") %in% old$phenotype))
old_rows <- data.frame(phenotype = old$phenotype, cases = old$cases,
                       controls = old$controls, cases_max = old$cases,
                       controls_max = old$controls,
                       includes_ukbb = FALSE,
                       n_snps = NA_integer_,
                       effect_allele = "A1",
                       effect_allele_evidence = "Ricopili daner / PGC VCF: effect is for A1",
                       filename = old$filename,
                       source = "01_harmonize_sumstats.R", stringsAsFactors = FALSE)

all_rows <- rbind(old_rows, rows)
info <- all_rows[, c("phenotype", "cases", "controls", "filename")]
f_info <- file.path(ADIR, "multipair.info.txt")
write.table(info, f_info, row.names = FALSE, quote = FALSE, sep = "\t")

f_prov <- file.path(ADIR, "multipair.provenance.txt")
write.table(all_rows, f_prov, row.names = FALSE, quote = FALSE, sep = "\t")

cat("\n================ COMBINED INPUT INFO ================\n")
print(info)
cat("\nwrote", f_info, "\n     ", f_prov, "\n")

Ntot <- info$cases + info$controls
cat(sprintf("\nSmallest panel N: %s at %d -> max.K = %d, K/N at K=700 = %.4f\n",
            info$phenotype[which.min(Ntot)], min(Ntot),
            floor(0.75 * min(Ntot)), 700 / min(Ntot)))
if (min(Ntot) < 6100)
  cat("WARNING: smallest N is below ~6,100. The single-pass assumption fails.\n")
cat("\nDONE.\n")
