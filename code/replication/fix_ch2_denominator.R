# fix_ch2_denominator.R -- one-off correction.
recompute <- function(audit_path, t1, t2) {
  r <- read.csv(audit_path, stringsAsFactors = FALSE)
  r$processed <- as.logical(r$processed)
  pro <- sum(r$processed)
  cn1 <- paste0("p_", t1); cn2 <- paste0("p_", t2)
  c1 <- round(100 * sum(r$processed & is.na(r[[cn1]])) / pro, 1)
  c2 <- round(100 * sum(r$processed & is.na(r[[cn2]])) / pro, 1)
  list(ch2_loss_1 = c1, ch2_loss_2 = c2)
}

patch_csv <- function(csv_path, pair_traits) {
  d <- read.csv(csv_path, stringsAsFactors = FALSE)
  for (pn in names(pair_traits)) {
    i <- which(d$pair == pn)
    if (!length(i)) next
    tr <- pair_traits[[pn]]
    v <- recompute(tr$audit, tr$t1, tr$t2)
    cat(sprintf("%s: ch2_loss_1 %.1f -> %.1f, ch2_loss_2 %.1f -> %.1f\n",
                pn, d$ch2_loss_1[i], v$ch2_loss_1, d$ch2_loss_2[i], v$ch2_loss_2))
    d$ch2_loss_1[i] <- v$ch2_loss_1
    d$ch2_loss_2[i] <- v$ch2_loss_2
  }
  write.csv(d, csv_path, row.names = FALSE)
}

MP <- "../../../analysis/out_R/multipair"
PANEL <- "../../results/panel"

patch_csv(file.path(PANEL, "PANEL_SUMMARY_P1.csv"),
          list(P1 = list(audit = file.path(MP, "P1", "attrition_audit.csv"),
                          t1 = "SCZ", t2 = "BIP")))

patch_csv(file.path(PANEL, "PANEL_SUMMARY_P2.csv"),
          list(P2a = list(audit = file.path(MP, "P2a", "attrition_audit.csv"),
                           t1 = "MDD", t2 = "BIP"),
               P2b = list(audit = file.path(MP, "P2b", "attrition_audit.csv"),
                           t1 = "MDD_noUKB", t2 = "BIP_noUKB")))

patch_csv(file.path(PANEL, "PANEL_SUMMARY_P3.csv"),
          list(P3 = list(audit = file.path(MP, "P3", "attrition_audit.csv"),
                          t1 = "ASD", t2 = "ADHD")))

patch_csv(file.path(PANEL, "PANEL_SUMMARY_P5.csv"),
          list(P5 = list(audit = file.path(MP, "P5", "attrition_audit.csv"),
                          t1 = "AN", t2 = "ASD")))

patch_csv(file.path(PANEL, "PANEL_SUMMARY_P4.csv"),
          list(P4 = list(audit = "../../results/applied_tierB/attrition_audit.csv",
                          t1 = "ASD", t2 = "SCZ")))

cat("\nDONE.\n")
