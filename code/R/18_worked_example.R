#!/usr/bin/env Rscript
# 18_worked_example.R -- regenerates the three worked loci of Table 12.

.args <- commandArgs(trailingOnly = FALSE)
.self <- sub("^--file=", "", grep("^--file=", .args, value = TRUE)[1])
BASE  <- normalizePath(file.path(dirname(.self), "..", ".."))
source(file.path(BASE, "code/R/lavagate.R"))

R_E   <- 0.21388
OUT   <- file.path(BASE, "results/worked_example")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

a  <- utils::read.csv(file.path(BASE, "results/panel/P1/attrition_audit.csv"))
rp <- a[!is.na(a$rho_hat), ]
o  <- lava_gate_correct(rp$rho_hat, rp$K, rp$p_SCZ, rp$p_BIP, r_e = R_E,
                        rho_lower = rp$rho_lower, rho_upper = rp$rho_upper,
                        aggregate = FALSE)

pick <- function(loc) {
  i <- which(rp$locus == loc)
  list(locus = rp$locus[i], chr = rp$chr[i], K = rp$K[i],
       p_weak = rp$p_SCZ[i], p_strong = rp$p_BIP[i],
       rho_hat = rp$rho_hat[i],
       rho_lower = rp$rho_lower[i], rho_upper = rp$rho_upper[i],
       ci_width = rp$rho_upper[i] - rp$rho_lower[i],
       theta1 = o$theta1[i], theta2 = o$theta2[i],
       slope = o$slope[i], intercept = o$intercept[i],
       at_floor = o$at_floor[i], rho_corrected = o$rho_corrected[i],
       ci_width_corrected = o$ci_width_corrected[i],
       width_inflation = o$width_inflation[i],
       naive_divide_only = rp$rho_hat[i] / o$slope[i])
}

res <- list(pair = "P1", r_e = R_E, n_reporting = nrow(rp),
            n_at_floor = sum(o$at_floor),
            pct_at_floor = round(100 * mean(o$at_floor), 1),
            A = pick(874), B = pick(38), C = pick(168))

writeLines(jsonlite::toJSON(res, auto_unbox = TRUE, pretty = TRUE, digits = 8),
           file.path(OUT, "p1_worked_example.json"))
cat(sprintf("reporting %d ; at floor %d (%.1f%%)\n",
            res$n_reporting, res$n_at_floor, res$pct_at_floor))
for (n in c("A", "B", "C")) {
  w <- res[[n]]
  cat(sprintf("%s locus %d chr%d K=%d slope=%.4f int=%.4f corr=%s naive=%s\n",
              n, w$locus, w$chr, w$K, w$slope, w$intercept,
              format(round(w$rho_corrected, 4)), format(round(w$naive_divide_only, 4))))
}
