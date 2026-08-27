#!/usr/bin/env Rscript
# lavagate.R -- gate-aware correction for LAVA local genetic correlations.

.lg_env <- new.env(parent = emptyenv())

lg_nref_scale <- function(K, n_ref = 1e5) (n_ref - K - 1) / (n_ref - 1)

lg_gate_threshold <- function(K, alpha = 0.05) qchisq(1 - alpha, K)

lg_poisson_weights <- function(lam, tol = 1e-14) {
  half <- lam / 2
  jmax <- as.integer(max(60, half + 12 * sqrt(max(half, 1))))
  j <- 0:jmax
  logw <- -half + j * log(half + 1e-300) - lgamma(j + 1)
  w <- exp(logw)
  if (max(w) > 0) {
    keep <- w > tol * max(w)
    j <- j[keep]; w <- w[keep]
  }
  list(j = j, w = w / sum(w))
}

lg_ncx2_upper_moments <- function(K, lam, t) {
  pw <- lg_poisson_weights(lam)
  m <- K + 2 * pw$j
  list(p_tail  = sum(pw$w * pchisq(t, m, lower.tail = FALSE)),
       e_trunc = sum(pw$w * m * pchisq(t, m + 2, lower.tail = FALSE)))
}

lg_delta <- function(K, theta, n_ref = 1e5, alpha = 0.05) {
  nref <- lg_nref_scale(K, n_ref)
  t <- lg_gate_threshold(K, alpha) / nref
  mm <- lg_ncx2_upper_moments(K, K * theta, t)
  mm$e_trunc / mm$p_tail - K * (1 + theta)
}

lg_grid <- function(K, n_ref = 1e5, alpha = 0.05, n = 3000, theta_max = 5) {
  key <- paste(as.integer(K), n_ref, alpha, n, sep = "|")
  if (!is.null(.lg_env[[key]])) return(.lg_env[[key]])
  th <- c(0, exp(seq(log(1e-7), log(theta_max), length.out = n)))
  D <- vapply(th, function(t) lg_delta(K, max(t, 1e-12), n_ref, alpha), numeric(1))
  g <- th + D / K
  ok <- is.finite(g) & is.finite(D)
  th <- th[ok]; g <- g[ok]; D <- D[ok]
  o <- order(g)
  out <- list(theta = th[o], g = g[o], D = D[o])
  .lg_env[[key]] <- out
  out
}


lava_theta_obs_from_p <- function(p, K, n_ref = 1e5) {
  nref <- lg_nref_scale(K, n_ref)
  qchisq(p, K, lower.tail = FALSE) / K - nref
}

lava_invert_theta <- function(theta_obs, K, n_ref = 1e5, alpha = 0.05) {
  gr <- lg_grid(K, n_ref, alpha)
  theta_obs <- as.numeric(theta_obs)
  at_floor <- theta_obs <= gr$g[1]
  th <- approx(gr$g, gr$theta, xout = theta_obs, rule = 2)$y
  th[at_floor] <- 0
  D <- approx(gr$theta, gr$D, xout = th, rule = 2)$y
  data.frame(theta = th, Delta = D, at_floor = at_floor)
}

lava_slope_intercept <- function(K, theta1, theta2, D1, D2, r_e = 0) {
  cc <- sqrt(theta1 * theta2)
  d1 <- theta1 + D1 / K
  d2 <- theta2 + D2 / K
  den <- sqrt(d1 * d2)
  s <- 1 + (D1 / (1 + 2 * theta1) + D2 / (1 + 2 * theta2)) / K
  list(slope = cc * s / den,
       intercept = (r_e / K) * (D1 * (1 + theta1) / (1 + 2 * theta1) +
                                D2 * (1 + theta2) / (1 + 2 * theta2)) / den)
}

lava_gate_correct <- function(rho_hat, K, p1, p2, r_e = 0,
                              rho_lower = NULL, rho_upper = NULL,
                              n_ref = 1e5, alpha = 0.05, aggregate = TRUE) {
  n <- length(rho_hat)
  stopifnot(length(K) %in% c(1, n), length(p1) %in% c(1, n), length(p2) %in% c(1, n))
  K <- rep_len(as.integer(K), n); p1 <- rep_len(p1, n); p2 <- rep_len(p2, n)

  t1o <- lava_theta_obs_from_p(p1, K, n_ref)
  t2o <- lava_theta_obs_from_p(p2, K, n_ref)

  out <- data.frame(rho_hat = rho_hat, K = K, p1 = p1, p2 = p2,
                    theta1 = NA_real_, theta2 = NA_real_,
                    slope = NA_real_, intercept = NA_real_,
                    at_floor = NA, rho_corrected = NA_real_,
                    rho_corrected_unclipped = NA_real_)
  for (kk in unique(K)) {
    ix <- which(K == kk)
    a <- lava_invert_theta(t1o[ix], kk, n_ref, alpha)
    b <- lava_invert_theta(t2o[ix], kk, n_ref, alpha)
    si <- lava_slope_intercept(kk, a$theta, b$theta, a$Delta, b$Delta, r_e)
    rc <- (rho_hat[ix] - si$intercept) / si$slope
    rc[!is.finite(rc)] <- NA_real_
    out$theta1[ix] <- a$theta; out$theta2[ix] <- b$theta
    out$slope[ix] <- si$slope; out$intercept[ix] <- si$intercept
    out$at_floor[ix] <- a$at_floor | b$at_floor
    out$rho_corrected_unclipped[ix] <- rc
    out$rho_corrected[ix] <- pmax(-1, pmin(1, rc))
  }
  out$rho_corrected[out$slope <= 0] <- NA_real_
  out$rho_corrected_unclipped[out$slope <= 0] <- NA_real_

  if (!is.null(rho_lower) && !is.null(rho_upper)) {
    out$ci_width <- rho_upper - rho_lower
    inv <- ifelse(out$slope > 0, 1 / out$slope, NA_real_)
    out$ci_width_corrected <- out$ci_width * inv
    out$width_inflation <- inv
  }

  if (aggregate) {
    Kb <- as.integer(round(mean(K)))
    a <- lava_invert_theta(mean(t1o), Kb, n_ref, alpha)
    b <- lava_invert_theta(mean(t2o), Kb, n_ref, alpha)
    si <- lava_slope_intercept(Kb, a$theta, b$theta, a$Delta, b$Delta, r_e)
    attr(out, "aggregate") <- list(
      K = Kb, slope = si$slope, intercept = si$intercept,
      mean_rho_hat = mean(rho_hat, na.rm = TRUE),
      mean_rho_corrected = if (si$slope > 0)
        (mean(rho_hat, na.rm = TRUE) - si$intercept) / si$slope else NA_real_)
  }
  out
}

lava_gate_diagnostic <- function(K, p1, p2, r_e = 0, n_ref = 1e5, alpha = 0.05,
                                 n_null = 40000, seed = 20260823,
                                 c1_max_fold = 1.25, c2_min_z = 5,
                                 c2_max_frac_below = 0.40, c3_max_frac_null = 0.40,
                                 p1_all = NULL, p2_all = NULL) {
  n <- length(K); K <- as.integer(K)
  t1o <- lava_theta_obs_from_p(p1, K, n_ref)
  t2o <- lava_theta_obs_from_p(p2, K, n_ref)

  strat <- function(ns) {
    b <- if (ns == 1) rep(0L, n)
         else {
           qs <- quantile(K, seq(0, 1, length.out = ns + 1), names = FALSE)
           pmin(pmax(findInterval(K, qs[-c(1, length(qs))]), 0L), ns - 1L)
         }
    sw <- 0
    for (g in unique(b)) {
      ix <- which(b == g)
      Kb <- max(25L, as.integer(round(mean(K[ix]) / 25) * 25))
      a <- lava_invert_theta(mean(t1o[ix]), Kb, n_ref, alpha)
      d <- lava_invert_theta(mean(t2o[ix]), Kb, n_ref, alpha)
      si <- lava_slope_intercept(Kb, a$theta, d$theta, a$Delta, d$Delta, r_e)
      sw <- sw + si$slope * length(ix)
    }
    sw / n
  }
  sl <- vapply(c(1, 2, 3, 4, 6, 8), strat, numeric(1))
  fold <- if (min(sl) > 0) max(sl) / min(sl) else Inf
  C1 <- list(passed = fold < c1_max_fold, slope_fold_change = fold,
             slopes = sl, threshold = c1_max_fold)

  set.seed(seed)
  kbin <- round(K / 50) * 50; kbin[kbin == 0] <- 50
  floor_one <- function(tobs) {
    zs <- ws <- fb <- numeric(0)
    for (Kb in sort(unique(kbin))) {
      nref <- lg_nref_scale(Kb, n_ref)
      t <- lg_gate_threshold(Kb, alpha) / nref
      u <- rchisq(n_null, Kb)
      u <- u[u > t]
      thn <- u / Kb - nref
      s <- tobs[kbin == Kb]
      if (!length(s)) next
      zs <- c(zs, (mean(s) - mean(thn)) / (sd(thn) / sqrt(length(s))))
      fb <- c(fb, mean(s < mean(thn)))
      ws <- c(ws, length(s))
    }
    list(pooled_z = sum(zs * ws) / sum(ws), frac_below = sum(fb * ws) / sum(ws))
  }
  f1 <- floor_one(t1o); f2 <- floor_one(t2o)
  C2 <- list(passed = (f1$pooled_z >= c2_min_z && f1$frac_below <= c2_max_frac_below &&
                       f2$pooled_z >= c2_min_z && f2$frac_below <= c2_max_frac_below),
             trait1 = f1, trait2 = f2,
             threshold_z = c2_min_z, threshold_frac_below = c2_max_frac_below)

  C3 <- list(passed = NA, note = paste(
    "Not evaluated: supply p1_all and p2_all, the univariate P-values at every",
    "locus where BOTH traits were estimable, to compute the gate pass rate."))
  if (!is.null(p1_all) && !is.null(p2_all)) {
    fr <- function(p) {
      npass <- sum(p < alpha, na.rm = TRUE)
      min(1, alpha * length(p) / max(npass, 1))
    }
    fr1 <- fr(p1_all); fr2 <- fr(p2_all)
    C3 <- list(passed = (fr1 <= c3_max_frac_null && fr2 <= c3_max_frac_null),
               frac_null_trait1 = fr1, frac_null_trait2 = fr2,
               enrichment_trait1 = mean(p1_all < alpha, na.rm = TRUE) / alpha,
               enrichment_trait2 = mean(p2_all < alpha, na.rm = TRUE) / alpha,
               threshold = c3_max_frac_null)
  }

  usable <- isTRUE(C1$passed) && isTRUE(C2$passed) && isTRUE(C3$passed)
  failed <- c(if (!isTRUE(C1$passed)) "C1 stability",
              if (!isTRUE(C2$passed)) "C2 floor distance",
              if (!isTRUE(C3$passed)) "C3 gate signal")
  list(usable = usable, C1 = C1, C2 = C2, C3 = C3, failed_criteria = failed,
       statement = if (usable)
         "A corrected magnitude MAY be reported." else
         paste0("NOT APPLICABLE. Report the uncorrected estimate with the ",
                "attenuation stated as a bound. Failed: ",
                paste(failed, collapse = ", ")))
}

lava_gate_correct_locus <- function(univ, bivar, K, r_e = 0, ...) {
  stopifnot(nrow(bivar) >= 1)
  p1 <- univ$p[match(bivar$phen1[1], univ$phen)]
  p2 <- univ$p[match(bivar$phen2[1], univ$phen)]
  if (is.na(p1) || is.na(p2))
    stop("A phenotype in `bivar` is absent from `univ`. process.locus() drops ",
         "phenotypes individually on a negative variance estimate; check for ",
         "that before correcting.")
  lava_gate_correct(bivar$rho[1], K, p1, p2, r_e,
                    rho_lower = bivar$rho.lower[1], rho_upper = bivar$rho.upper[1],
                    ...)
}
