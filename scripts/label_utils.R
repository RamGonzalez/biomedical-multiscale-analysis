# =============================================================================
# Figure labels
#
# Builds the "<drug> - <rsid>" labels used in the figures (drug_label()).
#
# Inputs:  snp_label and rsid of each variant
# Outputs: label string
# Usage:   source("scripts/label_utils.R")
# =============================================================================

drug_label <- function(snp_label, rsid) {
  hand <- c(rs8099917   = "Interferons, Ribavirin, Telaprevir",
            rs116855232 = "Mercaptopurine, Azathioprine")
  lbl <- ifelse(rsid %in% names(hand), hand[rsid], snp_label)
  lbl <- sub(" response - .*$", "", lbl)
  lbl <- gsub(" and ", ", ", lbl, fixed = TRUE)
  lbl <- vapply(strsplit(lbl, ", ", fixed = TRUE), function(p) {
    p <- trimws(p); paste0(toupper(substring(p, 1, 1)), substring(p, 2), collapse = ", ")
  }, character(1))
  ifelse(is.na(snp_label) | snp_label == "", rsid, paste(lbl, rsid, sep = " - "))
}
