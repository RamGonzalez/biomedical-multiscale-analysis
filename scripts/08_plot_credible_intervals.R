#!/usr/bin/env Rscript
# =============================================================================
# 08 - Interval figure
#
# Plots the 95% interval of each estimate for the reported variants: one facet per variant, with rows for the populations grouped by region and the pooled regional estimate.
#
# Inputs:  results/freqs_final.tsv, results/biomedical_samples.tsv
# Outputs: results/ci_two_scales.pdf
# Usage:   Rscript 08_plot_credible_intervals.R [config.yaml]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr); library(ggplot2); library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
cfg  <- yaml::read_yaml(if (length(args) >= 1) args[1] else "config.yaml")
region_levels <- unlist(cfg$params$panel_b_order)
level <- cfg$params$ci_level

# ---- Colors ----
COL_POP <- "#2a78d6"; COL_REG <- "#eb6834"
INK_1 <- "#0b0b0b"; INK_2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; SURF <- "#fcfcfb"
LAB_POP <- "Population (panel a)"; LAB_REG <- "Region, pooled (panel b)"

dat <- fread("results/freqs_final.tsv", colClasses = list(character = "group_code"))
clst2group <- fread("results/biomedical_samples.tsv") %>% distinct(CLST_name, Group)

# ---- Population order ----
pop_order_a <- c(
  "Malagasy","Myanmar, Laos, Vietnam","Malay, Kota Kinabalu","Orang Asli","Taiwan",
  "Andaman Islands","Sulawesi and Ambon Islands","Sunda arc","Kankana-ey","Aeta","Agta",
  "East Papuan Highlands","West Papuan Highlands","Southern Papua New Guinea",
  "Madang, East Sepik","East New Britain","Kove, Mussau, Manus Island","Lavongai, Kuot",
  "Kavieng, Lihir","Northern Vanuatu","New Caledonia","Fiji","Guam","Chuuk",
  "Kapingamarangi","Pohnpei","Nauru","Palau","Majuro","Kiribati","Tuvalu Micronesia","Tuvalu Polynesia","Tonga",
  "Samoa","Atiu","Mauke","Tikehau","Raivavae","Rimatara","Rapa","South Marquesas","Mangareva")

source("scripts/label_utils.R")
snp_labels <- dat %>% distinct(SNP, rsid, snp_label) %>%
  mutate(facet_lab = drug_label(snp_label, rsid), sort_key = tolower(facet_lab))

# ---- Plot data ----
pops <- dat %>% filter(scale == "clst") %>%
  left_join(clst2group, by = c("group_name" = "CLST_name")) %>%
  mutate(region_i = match(Group, region_levels),
         pop_rank = match(group_name, pop_order_a),
         row_lab  = group_name, level = LAB_POP)
stopifnot(!anyNA(pops$region_i), !anyNA(pops$pop_rank))

regs <- dat %>% filter(scale == "group") %>%
  mutate(region_i = match(group_name, region_levels), pop_rank = Inf,
         row_lab = paste0(gsub("_", " ", group_name), " (pooled)"), level = LAB_REG)

plot_df <- bind_rows(pops, regs) %>%
  filter(!is.na(region_i)) %>%
  left_join(snp_labels %>% select(SNP, sort_key, facet_lab), by = "SNP") %>%
  mutate(level = factor(level, levels = c(LAB_POP, LAB_REG)))

row_order <- plot_df %>% distinct(row_lab, region_i, pop_rank) %>%
  arrange(region_i, pop_rank) %>% pull(row_lab)
n_rows <- length(row_order)
last_in_region <- cumsum(plot_df %>% distinct(row_lab, region_i) %>%
                           count(region_i) %>% arrange(region_i) %>% pull(n))
sep_y <- n_rows - last_in_region[-length(last_in_region)] + 0.5

# ---- Plot function ----
theme_ci <- theme_minimal(base_size = 8) +
  theme(plot.background  = element_rect(fill = SURF, colour = NA),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = GRID, linewidth = 0.25),
        axis.text  = element_text(colour = INK_2),
        axis.title = element_text(colour = INK_2),
        strip.text = element_text(colour = INK_1, size = 7, hjust = 0),
        plot.title = element_text(colour = INK_1, face = "bold", size = 11),
        plot.subtitle = element_text(colour = INK_2, size = 8),
        legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(colour = INK_2))

forest <- function(snp_set, title, outfile, ncol = 5) {
  d <- plot_df %>% filter(SNP %in% snp_set)
  facet_order <- snp_labels %>% filter(SNP %in% snp_set) %>% arrange(sort_key) %>% pull(facet_lab)
  d$facet_lab <- factor(d$facet_lab, levels = unique(facet_order))
  nrow_f <- ceiling(length(unique(facet_order)) / ncol)
  p <- ggplot(d, aes(effect_freq, row_lab, colour = level, shape = level)) +
    geom_hline(yintercept = sep_y, colour = GRID, linewidth = 0.3) +
    geom_linerange(data = d %>% filter(level == LAB_POP),
                   aes(xmin = ci_lower, xmax = ci_upper), orientation = "y",
                   linewidth = 0.3, alpha = 0.8) +
    geom_linerange(data = d %>% filter(level == LAB_REG),
                   aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", linewidth = 0.9) +
    geom_point(aes(size = level)) +
    scale_colour_manual(values = setNames(c(COL_POP, COL_REG), c(LAB_POP, LAB_REG))) +
    scale_shape_manual(values = setNames(c(16, 18), c(LAB_POP, LAB_REG))) +
    scale_size_manual(values = setNames(c(0.9, 2.6), c(LAB_POP, LAB_REG))) +
    scale_y_discrete(limits = rev(row_order)) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    facet_wrap(~ facet_lab, ncol = ncol) +
    labs(title = title,
         subtitle = sprintf(paste0("Effect-allele frequency with %d%% interval. Rows are grouped by region: each population (Jeffreys credible interval) and, last in each block,\n",
                                   "the pooled regional frequency (beta-binomial interval that includes variation between the region's populations; Taiwan and Madagascar have one population each, so they show its Jeffreys interval)."),
                            round(level * 100)),
         x = "Effect-allele frequency", y = NULL) +
    theme_ci +
    theme(axis.text.y = element_text(size = 5.5, colour = INK_1),
          panel.grid.major.y = element_line(colour = GRID, linewidth = 0.15),
          panel.grid.major.x = element_line(colour = GRID, linewidth = 0.25),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
  ggsave(outfile, p, width = 2.0 + 2.3 * ncol, height = 1.2 + (0.12 * n_rows + 0.5) * nrow_f,
         limitsize = FALSE)
  cat("Wrote", outfile, sprintf("(%d SNPs)\n", length(unique(facet_order))))
}

# ---- Plot ----
snp_set <- dat %>% filter(reported) %>% pull(SNP) %>% unique()
forest(snp_set, "Allele-frequency uncertainty at two scales - reported biomedical variants",
       "results/ci_two_scales.pdf")
