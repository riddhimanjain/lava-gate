# p4_panel_summary.R -- compute the same PANEL_SUMMARY_<GROUP>.csv row that
r <- read.csv("../../results/applied_tierB/attrition_audit.csv", stringsAsFactors = FALSE)
r$processed <- as.logical(r$processed)
r$gate_pass <- as.logical(r$gate_pass)

ph <- c("ASD", "SCZ")
cn_p <- paste0("p_", ph)
N_SNPS_COMMON <- 5774835

pro  <- sum(r$processed)
e1   <- sum(r$processed & !is.na(r[[cn_p[1]]]))
e2   <- sum(r$processed & !is.na(r[[cn_p[2]]]))
both <- sum(r$processed & !is.na(r[[cn_p[1]]]) & !is.na(r[[cn_p[2]]]))
gp   <- sum(r$gate_pass, na.rm = TRUE)
rep_ <- sum(!is.na(r$rho_hat))
plim <- sum(r$paramlim_drop, na.rm = TRUE)
det  <- !is.na(r$rho_lower) & (r$rho_lower > 0 | r$rho_upper < 0)
w    <- r$rho_upper - r$rho_lower
bonf <- sum(r$processed & r[[cn_p[1]]] < 0.05 / pro & r[[cn_p[2]]] < 0.05 / pro, na.rm = TRUE)
weaker <- ph[which.min(c(e1, e2))]

pan <- data.frame(
  pair = "P4", traits = paste(ph, collapse = " x "), group = "P4",
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

write.csv(pan, "../../results/panel/PANEL_SUMMARY_P4.csv", row.names = FALSE)
cat("wrote PANEL_SUMMARY_P4.csv\n")
print(pan)
