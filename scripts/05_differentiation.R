#!/usr/bin/env Rscript
# =============================================================================
# 05 - Frequency differences
#
# Tests each group against all other groups at the same scale, adjusts p-values (BH-FDR), and writes the frequency table and the tables of reported and not-reported variants.
#
# Inputs:  results/freqs_with_ci.tsv, results/regional_betabinomial.tsv, results/effect_variant_table.tsv, resources/clinical_allele_sources.tsv
# Outputs: results/freqs_final.tsv, results/table_reported_snps.tsv, results/table_reported_snps_summary.tsv, results/table_excluded_snps.tsv
# Usage:   Rscript 05_differentiation.R [config.yaml]
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(dplyr); library(yaml) })

args <- commandArgs(trailingOnly = TRUE)
cfg  <- yaml::read_yaml(if (length(args) >= 1) args[1] else "config.yaml")
P    <- cfg$params
p_thr <- P$p_threshold; fdr_alpha <- P$fdr_alpha
report_scales <- unlist(P$report_scales)
min_af <- P$report_min_af

fisher_p <- function(a, b, A, B) {
  fisher.test(matrix(round(c(a, b - a, A - a, (B - b) - (A - a))), nrow = 2))$p.value
}

dat <- fread("results/freqs_with_ci.tsv", colClasses = list(character = "group_code"))
bb  <- fread("results/regional_betabinomial.tsv")[
  , .(SNP, group_name = region, n_pops, rho, bb_lo, bb_hi, p_bb)]

# ---- Tested cells ----
out <- dat %>%
  mutate(has_data = !is.na(n) & n > 0 & !is.na(c)) %>%
  group_by(scale, SNP) %>%
  mutate(C_all = sum(c[has_data]), N_all = sum(n[has_data]),
         tested = has_data & C_all > 0 & C_all < N_all,
         freq_rest = ifelse(has_data, (C_all - c) / (N_all - n), NA_real_)) %>%
  ungroup()

# ---- Fisher's exact test ----
out$p_fisher <- NA_real_
t <- out$tested
out$p_fisher[t] <- mapply(fisher_p, out$c[t], out$n[t], out$C_all[t], out$N_all[t])

# ---- Regional intervals and tests ----
out <- out %>%
  left_join(bb, by = c("SNP", "group_name")) %>%
  mutate(
    is_region       = scale == "group",
    ci_pooled_lower = ifelse(is_region, ci_lower, NA_real_),
    ci_pooled_upper = ifelse(is_region, ci_upper, NA_real_),
    use_bb    = is_region & !is.na(n_pops) & n_pops > 1,
    ci_lower  = ifelse(use_bb, bb_lo, ci_lower),
    ci_upper  = ifelse(use_bb, bb_hi, ci_upper),
    ci_width  = ci_upper - ci_lower,
    ci_method = ifelse(use_bb, "betabinomial_profile", ci_method),
    test      = ifelse(use_bb, "betabinomial_lrt", "fisher_exact"),
    p_value   = ifelse(use_bb, ifelse(tested, p_bb, NA_real_), p_fisher))

# ---- FDR and stars ----
out <- out %>%
  group_by(scale) %>%
  mutate(q_value = { v <- rep(NA_real_, n()); ok <- !is.na(p_value)
                     v[ok] <- p.adjust(p_value[ok], method = "BH"); v },
         star = !is.na(p_value) & p_value < p_thr & q_value < fdr_alpha) %>%
  ungroup()

# ---- Reported set ----
starred <- out %>% filter(scale %in% report_scales, star) %>% pull(SNP) %>% unique()
common  <- out %>% filter(scale == "clst", effect_freq > min_af) %>% pull(SNP) %>% unique()
require_star <- isTRUE(P$report_require_star)
reported <- if (require_star) intersect(starred, common) else common
manual_excl <- unlist(P$report_exclude)
reported <- setdiff(reported, out$SNP[out$rsid %in% names(manual_excl)])
out <- out %>%
  mutate(reported = SNP %in% reported) %>%
  select(-has_data, -C_all, -N_all, -tested, -is_region, -use_bb, -bb_lo, -bb_hi, -p_bb)

fwrite(out, "results/freqs_final.tsv", sep = "\t", na = "NA")

# ---- Tables of reported variants ----
source("scripts/label_utils.R")
clin <- fread(cfg$inputs$clinical_allele_sources)[, .(rsid, clinical_allele, clinical_allele_source)]
chk <- merge(unique(out %>% filter(reported) %>% select(rsid, effect_allele)), clin, by = "rsid", all.x = TRUE)
stopifnot(!anyNA(chk$clinical_allele), all(toupper(chk$clinical_allele) == toupper(chk$effect_allele)))
eff <- fread("results/effect_variant_table.tsv")[, .(SNP, gene)]
scale_lab <- c(clst = "population", group = "region")

long <- out %>% filter(reported) %>%
  left_join(eff, by = "SNP") %>%
  left_join(clin, by = "rsid") %>%
  transmute(category, rsid, SNP, gene, snp_label,
            figure_label = drug_label(snp_label, rsid), REF, effect_allele,
            clinical_allele = effect_allele,
            clinical_allele_freq = round(effect_freq, 4),
            clinical_allele_source,
            scale = scale_lab[scale], group = group_name,
            effect_copies = c, called_chromosomes = n,
            effect_freq = round(effect_freq, 4),
            ci95_lower = round(ci_lower, 4), ci95_upper = round(ci_upper, 4),
            ci_method, freq_rest_of_scale = round(freq_rest, 4),
            test, p_value = signif(p_value, 3), q_value = signif(q_value, 3), star) %>%
  arrange(category, rsid, scale, group)
fwrite(long, "results/table_reported_snps.tsv", sep = "\t", na = "NA")

starred_list <- function(g, f, r) {
  if (!length(g)) return("")
  paste(sprintf("%s (%.2f vs %.2f)", g, f, r), collapse = "; ")
}
summ <- out %>% filter(reported) %>%
  left_join(eff, by = "SNP") %>%
  left_join(clin, by = "rsid") %>%
  mutate(figure_label = drug_label(snp_label, rsid), clinical_allele = effect_allele) %>%
  group_by(category, rsid, SNP, gene, snp_label, figure_label, effect_allele,
           clinical_allele, clinical_allele_source) %>%
  summarise(
    n_pop_starred    = sum(star & scale == "clst"),
    populations      = starred_list(group_name[star & scale == "clst"],
                                    effect_freq[star & scale == "clst"],
                                    freq_rest[star & scale == "clst"]),
    n_region_starred = sum(star & scale == "group"),
    regions          = starred_list(group_name[star & scale == "group"],
                                    effect_freq[star & scale == "group"],
                                    freq_rest[star & scale == "group"]),
    .groups = "drop") %>%
  arrange(category, rsid)
fwrite(summ, "results/table_reported_snps_summary.tsv", sep = "\t", na = "NA")

# ---- Tables of not-reported variants ----
excl <- out %>% filter(!reported) %>%
  left_join(eff, by = "SNP") %>%
  mutate(figure_label = drug_label(snp_label, rsid)) %>%
  group_by(category, rsid, SNP, gene, snp_label, figure_label, effect_allele) %>%
  summarise(
    n_pop_starred    = sum(star & scale == "clst"),
    n_region_starred = sum(star & scale == "group"),
    max_pop_freq     = round(max(effect_freq[scale == "clst"], na.rm = TRUE), 4),
    max_pop          = group_name[scale == "clst"][which.max(effect_freq[scale == "clst"])],
    starred_groups   = paste(c(group_name[star & scale == "clst"],
                               group_name[star & scale == "group"]), collapse = "; "),
    .groups = "drop") %>%
  mutate(reason = ifelse(rsid %in% names(manual_excl), unname(manual_excl[rsid]),
                  ifelse(n_pop_starred + n_region_starred == 0,
                         sprintf("Frequency never exceeds %g in any population; no significant group at either scale", min_af),
                         sprintf("Frequency never exceeds %g in any population; significant in at least one group", min_af)))) %>%
  rename(alt_allele = effect_allele, max_pop_alt_freq = max_pop_freq) %>%
  arrange(reason, category, rsid)
fwrite(excl, "results/table_excluded_snps.tsv", sep = "\t", na = "NA")
cat(sprintf("Excluded SNPs (tested, not reported): %d | %s\n", nrow(excl),
            paste(names(table(excl$reason)), table(excl$reason), sep = ": ", collapse = " | ")))

# ---- Summary ----
smry <- out %>% group_by(scale, test) %>%
  summarise(tests = sum(!is.na(p_value)),
            p_below_thr = sum(p_value < p_thr, na.rm = TRUE),
            q_below_alpha = sum(q_value < fdr_alpha, na.rm = TRUE),
            starred_tiles = sum(star),
            starred_snps = n_distinct(SNP[star]),
            starred_tiles_reported = sum(star & reported), .groups = "drop")
cat(sprintf("Populations and single-population regions: Fisher's exact. Other regions: beta-binomial LRT. BH-FDR per scale; star = p < %g AND q < %g\n",
            p_thr, fdr_alpha))
print(as.data.frame(smry))
cat(sprintf("Starred SNPs (>= 1 star at %s): %d of %d\n",
            paste(report_scales, collapse = "/"), length(starred), n_distinct(out$SNP)))
cat(sprintf("Reported SNPs (effect-allele freq > %g in >= 1 population%s): %d\n",
            min_af, if (require_star) ", and >= 1 star" else ", no star required", length(reported)))
cat("Wrote results/freqs_final.tsv, table_reported_snps.tsv, table_reported_snps_summary.tsv\n")
