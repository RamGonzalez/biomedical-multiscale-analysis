#!/usr/bin/env Rscript
# =============================================================================
# 04 - Credible intervals
#
# Jeffreys 95% credible interval for each variant and group.
#
# Inputs:  results/freq_clst.tsv, results/freq_group.tsv
# Outputs: results/freqs_with_ci.tsv
# Usage:   Rscript 04_credible_intervals.R [config.yaml]
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(dplyr); library(yaml) })

args <- commandArgs(trailingOnly = TRUE)
cfg  <- yaml::read_yaml(if (length(args) >= 1) args[1] else "config.yaml")
level <- cfg$params$ci_level
alpha <- 1 - level

# ---- Intervals ----
jeffreys_ci <- function(c, n, a = alpha) {
  lo <- rep(NA_real_, length(c)); hi <- lo
  ok <- !is.na(n) & n > 0 & !is.na(c)
  lo[ok] <- qbeta(a/2,     c[ok] + 0.5, n[ok] - c[ok] + 0.5)
  hi[ok] <- qbeta(1 - a/2, c[ok] + 0.5, n[ok] - c[ok] + 0.5)
  list(lo = lo, hi = hi)
}

add_ci <- function(path, scale_lab) {
  df <- fread(path, colClasses = list(character = "group_code"))
  ci <- jeffreys_ci(df$c, df$n)
  df %>% mutate(group_code = as.character(group_code),
                ci_method = "jeffreys", ci_level = level,
                ci_lower = ci$lo, ci_upper = ci$hi,
                ci_width = ci$hi - ci$lo)
}

# ---- Output ----
a <- add_ci("results/freq_clst.tsv",  "clst")
b <- add_ci("results/freq_group.tsv", "group")
out <- bind_rows(a, b)
fwrite(out, "results/freqs_with_ci.tsv", sep = "\t", na = "NA")

# ---- Summary ----
w <- out %>% filter(!is.na(ci_width)) %>%
  group_by(scale) %>%
  summarise(median_n = median(n), median_ci_width = median(ci_width),
            mean_ci_width = mean(ci_width), .groups = "drop")
cat("Jeffreys", level*100, "% credible intervals written to results/freqs_with_ci.tsv\n", sep = "")
cat("Median CrI width by scale (larger groups -> tighter intervals):\n")
print(as.data.frame(w))
