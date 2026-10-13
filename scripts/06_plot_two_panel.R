#!/usr/bin/env Rscript
# =============================================================================
# 06 - Two-panel figure
#
# Heatmap of effect-allele frequency for the reported variants in (a) populations and (b) regions, with gnomAD reference columns and stars for significant groups.
#
# Inputs:  results/freqs_final.tsv, results/effect_variant_table.tsv, gnomAD table (config.yaml)
# Outputs: results/two_panel_figure.pdf
# Usage:   Rscript 06_plot_two_panel.R [config.yaml]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(patchwork); library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
cfg  <- yaml::read_yaml(if (length(args) >= 1) args[1] else "config.yaml")
region_levels <- unlist(cfg$params$region_levels)
P <- cfg$params

IN  <- cfg$inputs
dat <- fread("results/freqs_final.tsv", colClasses = list(character = "group_code"))

# ---- gnomAD reference columns ----
gnomad_levels <- c("GNOMAD_EAS","GNOMAD_AFR","GNOMAD_NFE")
load_gnomad_effect <- function() {
  cohort_map <- c("African/African_American" = "GNOMAD_AFR",
                  "East_Asian"               = "GNOMAD_EAS",
                  "European_(non-Finnish)"    = "GNOMAD_NFE")
  eff <- fread("results/effect_variant_table.tsv")[
    , .(SNP, eff_allele = toupper(effect_allele), eff_ref = toupper(REF))]
  g <- fread(IN$gnomad, na.strings = c("", "NA"))
  g <- g[cohort %in% names(cohort_map)]
  g[, group_name := cohort_map[cohort]]
  g <- merge(g, eff, by = "SNP")
  g[, `:=`(ALT = toupper(ALT), REF = toupper(REF),
           AN = suppressWarnings(as.numeric(AN)),
           AF = suppressWarnings(as.numeric(AF)))]
  g[, effect_freq := fifelse(is.na(AN) | AN == 0, NA_real_,
                       fifelse(ALT == eff_allele, AF,
                        fifelse(REF == eff_allele, 1 - AF, NA_real_)))]
  g[, star := FALSE]
  g[, .(SNP, group_name, effect_freq, star)]
}
gnomad_ref <- load_gnomad_effect()

# ---- Panel a order ----
pop_order_a <- c(
  "Malagasy","Myanmar, Laos, Vietnam","Malay, Kota Kinabalu","Orang Asli","Taiwan",
  "Andaman Islands","Sulawesi and Ambon Islands","Sunda arc","Kankana-ey","Aeta","Agta",
  "East Papuan Highlands","West Papuan Highlands","Southern Papua New Guinea",
  "Madang, East Sepik","East New Britain","Kove, Mussau, Manus Island","Lavongai, Kuot",
  "Kavieng, Lihir","Northern Vanuatu","New Caledonia","Fiji","Guam","Chuuk",
  "Kapingamarangi","Pohnpei","Nauru","Palau","Majuro","Kiribati","Tuvalu Micronesia","Tuvalu Polynesia","Tonga",
  "Samoa","Atiu","Mauke","Tikehau","Raivavae","Rimatara","Rapa","South Marquesas","Mangareva")

# ---- Row labels ----
source("scripts/label_utils.R")
row_labels <- dat %>%
  distinct(SNP, rsid, snp_label) %>%
  mutate(uniq_label = drug_label(snp_label, rsid)) %>%
  distinct(SNP, uniq_label)

FILL_COLS <- c("#29CCE2","#1E77F4","#FABC3C","#FF064C","#D11149")

# ---- Panels ----
make_panel <- function(df_scale, x_levels, row_levels, title, show_y) {
  d <- df_scale %>%
    left_join(row_labels, by = "SNP") %>%
    mutate(group_name = factor(group_name, levels = x_levels),
           uniq_label = factor(uniq_label, levels = rev(row_levels))) %>%
    filter(!is.na(group_name), !is.na(uniq_label))
  d <- d %>% complete(uniq_label, group_name)
  ggplot(d, aes(group_name, uniq_label)) +
    geom_tile(aes(fill = effect_freq), colour = "white", linewidth = 0.3) +
    scale_fill_gradientn(colours = FILL_COLS,
                         values = scales::rescale(c(0,0.25,0.5,0.75,1)),
                         limits = c(0,1), na.value = "#bfc3c9",
                         name = "Effect-allele\nfrequency") +
    geom_text(data = subset(d, star %in% TRUE), aes(label = "*"),
              size = 4.5, colour = "black") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = if (show_y) element_text(size = 6.5) else element_blank(),
          plot.title  = element_text(size = 9, face = "bold"),
          panel.grid  = element_blank())
}

snp_set <- dat %>% filter(reported) %>% pull(SNP) %>% unique()
a <- dat %>% filter(scale == "clst",  SNP %in% snp_set)
b <- dat %>% filter(scale == "group", SNP %in% snp_set)
gref <- gnomad_ref %>% filter(SNP %in% snp_set)
a <- bind_rows(a, gref); b <- bind_rows(b, gref)
row_levels <- row_labels %>% filter(SNP %in% snp_set) %>%
  arrange(tolower(uniq_label)) %>% pull(uniq_label) %>% unique()
xa <- c(intersect(pop_order_a, unique(a$group_name)), gnomad_levels)
xb <- c(intersect(unlist(cfg$params$panel_b_order), unique(b$group_name)), gnomad_levels)

# ---- Figure ----
n_a <- sum(!is.na(dat$p_value[dat$scale == "clst"]))
n_b <- sum(!is.na(dat$p_value[dat$scale == "group"]))
subtitle <- sprintf(paste0(
  "* p < %g and BH-FDR q < %g, each group vs all others at its scale. a) Fisher's exact test; b) beta-binomial test that treats each population as an independent observation (Fisher's exact for Taiwan and Madagascar, one population each).\n",
  "FDR computed separately per panel: a) %s tests, b) %d tests. Shown: variants with effect-allele frequency > %g in at least one population.\n",
  "gnomAD columns are reference only."),
  P$p_threshold, P$fdr_alpha, format(n_a, big.mark = ","), n_b, P$report_min_af)

pa <- make_panel(a, xa, row_levels, "a) OGVP populations + gnomAD", show_y = TRUE)
pb <- make_panel(b, xb, row_levels, "b) Regions + gnomAD",          show_y = FALSE)
fig <- pa + pb +
  plot_layout(widths = c(length(xa), length(xb)), guides = "collect") +
  plot_annotation(
    title = "Biomedical variants: effect-allele frequencies at two scales",
    subtitle = subtitle,
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 8, colour = "#52514e")))
h <- max(4, 1.9 + 0.19 * length(row_levels))
ggsave("results/two_panel_figure.pdf", fig, width = 13, height = h, limitsize = FALSE)
cat("Wrote results/two_panel_figure.pdf", sprintf("(%d SNPs)\n", length(row_levels)))
