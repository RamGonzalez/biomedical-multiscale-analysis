#!/usr/bin/env Rscript
# =============================================================================
# 00 - Sample manifest
#
# Lists the samples used, with their population and region (Tuvalu split in two).
#
# Inputs:  genotype .fam, sample groups, population index, Tuvalu split (config.yaml)
# Outputs: results/biomedical_samples.tsv
# Usage:   Rscript 00_sample_manifest.R [config.yaml]
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(dplyr); library(yaml) })

args <- commandArgs(trailingOnly = TRUE)
cfg  <- yaml::read_yaml(if (length(args) >= 1) args[1] else "config.yaml")
IN   <- cfg$inputs
dir.create("results", showWarnings = FALSE)

# ---- Read inputs ----
fam <- fread(paste0(IN$geno_prefix, ".fam"), header = FALSE,
             col.names = c("CLST","fam_IID","V3","V4","sex","pheno")) %>%
  mutate(CLST = as.integer(CLST),
         sample_id = sub("^.*_", "", fam_IID),
         population_fam = sub("_[^_]*$", "", fam_IID))

clst_index <- fread(IN$clst_index) %>% distinct(CLST, CLST_name)

popinfo <- fread(IN$popinfo) %>%
  mutate(IID = as.character(IID)) %>%
  select(IID, Population, Ethnicity, Location, Region, Group, Country)

# ---- Build manifest ----
manifest <- fam %>%
  left_join(clst_index, by = "CLST") %>%
  left_join(popinfo, by = c("sample_id" = "IID")) %>%
  transmute(sample_id, fam_IID,
            CLST, CLST_name,
            Group,
            Region,
            Population, Country, Location, Ethnicity,
            population_fam) %>%
  arrange(CLST, sample_id)

# ---- Tuvalu split ----
TS <- cfg$params$tuvalu_split
split <- fread(IN$tuvalu_split, header = FALSE, colClasses = "character")
split <- data.table(sample_id = split$V2, split_label = split$V3)
stopifnot(!anyDuplicated(split$sample_id))
is_tuv <- manifest$CLST == TS$parent_clst
stopifnot(setequal(manifest$sample_id[is_tuv], split$sample_id))
is_mic <- is_tuv & manifest$sample_id %in% split$sample_id[split$split_label == TS$micronesia_label]
is_pol <- is_tuv & !is_mic
manifest$CLST[is_mic]      <- TS$micronesia$clst;  manifest$CLST[is_pol]      <- TS$polynesia$clst
manifest$CLST_name[is_mic] <- TS$micronesia$name;  manifest$CLST_name[is_pol] <- TS$polynesia$name
manifest$Group[is_mic]     <- TS$micronesia$group; manifest$Group[is_pol]     <- TS$polynesia$group
manifest <- manifest %>% arrange(CLST, sample_id)
cat(sprintf("Tuvalu split: %d Tuvalu Micronesia, %d Tuvalu Polynesia\n", sum(is_mic), sum(is_pol)))

# ---- Write manifest ----
fwrite(manifest, "results/biomedical_samples.tsv", sep = "\t", na = "NA")

cat(sprintf("Wrote results/biomedical_samples.tsv  (%d samples)\n", nrow(manifest)))
cat(sprintf("Unmapped to popinfo: %d\n", sum(is.na(manifest$Population))))
cat("\nSamples per Group (panel b):\n")
print(manifest %>% count(Group, sort = TRUE) %>% as.data.frame())
cat(sprintf("\nDistinct CLST populations: %d\n", n_distinct(manifest$CLST)))
