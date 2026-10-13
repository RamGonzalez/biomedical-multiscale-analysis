#!/usr/bin/env Rscript
# =============================================================================
# 01 - Effect variants
#
# Builds the variant table (GenomeAsia and ClinVar) with the allele with the clinical effect of each variant.
#
# Inputs:  annotation, ClinVar table, GenomeAsia index, genotype .bim (config.yaml)
# Outputs: results/effect_variant_table.tsv
# Usage:   Rscript 01_define_effect_variants.R [config.yaml]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr); library(yaml)
})

args    <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else "config.yaml"
cfg <- yaml::read_yaml(cfg_path)
IN  <- cfg$inputs
outdir <- "results"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out_table <- file.path(outdir, "effect_variant_table.tsv")

# ---- Variant lists and helpers ----
rsIDs <- c("rs3130690","rs1042389","rs12721649","rs2279341","rs2279344","rs2279345",
           "rs28399499","rs34826503","rs35303484","rs36079186","rs3786547","rs3826711",
           "rs4803418","rs4803419","rs707265","rs8192709","rs8192719","rs35599367",
           "rs10264272","rs776746","rs12769205","rs4244285","rs4986893","rs28399504",
           "rs12248560","rs192154563","rs1065852","rs1081003","rs1135840","rs28371725",
           "rs769258","rs3892097","rs1799853","rs1057910","rs2108622","rs3093105",
           "rs3918290","rs67376798","rs1142345","rs75543815","rs56161402","rs2306283",
           "rs4149056","rs4148323","rs35350960","rs4124874","rs17708472","rs9923231")

drop_snps <- c("19_41518221","12_21331549","10_96540410","10_96522463")

clinical_allele_is_ref <- c("rs776746")

is_palindromic <- function(a, b) {
  comp <- c(A="T", T="A", C="G", G="C")
  au <- toupper(a); bu <- toupper(b)
  ok <- au %in% names(comp) & bu %in% names(comp)
  res <- rep(NA, length(au))
  res[ok] <- comp[au[ok]] == bu[ok]
  as.logical(res)
}

# ---- OGVP genotyped variants ----
bim <- fread(paste0(IN$geno_prefix, ".bim"), header = FALSE,
             select = 2, col.names = "SNP")
ogvp_ids <- unique(bim$SNP)
message("OGVP genotyped variants (.bim): ", length(ogvp_ids))

# ---- GenomeAsia variants ----
annovar <- fread(IN$annovar_clean,
                 select = c("Chr","Start","Ref","Alt","avsnp151",
                            "Gene.refGene","Func.refGene","ExonicFunc.refGene")) %>%
  mutate(SNP = paste0(Chr, "_", Start))

ga_index <- fread(IN$genomeasia_index) %>% rename(avsnp151 = rsid)

ga <- annovar %>%
  filter(avsnp151 %in% rsIDs) %>%
  left_join(ga_index, by = "avsnp151") %>%
  transmute(
    SNP, Chr = as.integer(Chr), Position = as.integer(Start),
    rsid = avsnp151, REF = Ref, ALT = Alt,
    category = "genomeasia_pharmacogenomic",
    source   = "genomeasia",
    gene = Gene.refGene, func = Func.refGene, exonic_func = ExonicFunc.refGene,
    snp_label = Drug
  ) %>% distinct()

# ---- ClinVar variants ----
clinvar <- fread(IN$clinvar) %>%
  filter(ReviewStatus.y != "criteria provided, multiple submitters, no conflicts") %>%
  transmute(
    SNP = snpid,
    REF = V5, ALT = V6,
    clin_sig = ClinicalSignificance.y,
    review_status = ReviewStatus.y,
    review_status_sub = ReviewStatus.x,
    phenotype = ReportedPhenotypeInfo,
    rs_dbsnp = `RS# (dbSNP)`,
    gene = GeneSymbol
  ) %>%
  filter(!SNP %in% drop_snps) %>%
  filter(!SNP %in% ga$SNP)

cv <- clinvar %>%
  tidyr::separate(SNP, into = c("Chr","Position"), sep = "_", remove = FALSE) %>%
  group_by(SNP) %>%
  mutate(
    lab = if_else(review_status_sub %in% c("reviewed by expert panel","practice guideline"),
                  phenotype, NA_character_)
  ) %>%
  ungroup() %>%
  mutate(
    Chr = as.integer(Chr), Position = as.integer(Position),
    rsid = if_else(!is.na(rs_dbsnp) & rs_dbsnp != -1, paste0("rs", rs_dbsnp), NA_character_),
    snp_label = str_remove(lab, "^[^:]*:"),
    category = paste0("clinvar_", str_replace_all(tolower(clin_sig), " ", "_")),
    source = "clinvar",
    func = NA_character_, exonic_func = NA_character_
  ) %>%
  group_by(SNP, Chr, Position, REF, ALT, category, source, gene, func, exonic_func) %>%
  summarise(
    rsid = dplyr::first(na.omit(rsid)),
    snp_label = {
      v <- unique(na.omit(snp_label)); if (length(v)) str_c(v, collapse = "|") else NA_character_
    },
    .groups = "drop"
  ) %>%
  mutate(rsid = if_else(is.na(rsid), NA_character_, rsid))

# ---- Combine variants ----
common_cols <- c("SNP","Chr","Position","rsid","REF","ALT",
                 "category","source","gene","func","exonic_func","snp_label")
effect_tbl <- bind_rows(ga %>% select(all_of(common_cols)),
                        cv %>% select(all_of(common_cols))) %>%
  mutate(
    REF = toupper(REF), ALT = toupper(ALT),
    effect_allele = if_else(rsid %in% clinical_allele_is_ref, REF, ALT),
    other_allele  = if_else(rsid %in% clinical_allele_is_ref, ALT, REF),
    palindromic   = is_palindromic(REF, ALT),
    in_ogvp       = SNP %in% ogvp_ids
  ) %>%
  arrange(category, Chr, Position) %>%
  distinct()

# ---- Write table ----
n_ga  <- sum(effect_tbl$source == "genomeasia")
n_cv  <- sum(effect_tbl$source == "clinvar")
n_pal <- sum(effect_tbl$palindromic, na.rm = TRUE)
n_off <- sum(!effect_tbl$in_ogvp)

message(sprintf("Effect variants: %d total | GA %d | ClinVar %d", nrow(effect_tbl), n_ga, n_cv))
message(sprintf("  palindromic (strand-ambiguous): %d", n_pal))
message(sprintf("  NOT genotyped in OGVP (dropped downstream): %d", n_off))
message("  categories:")
print(effect_tbl %>% count(category) %>% as.data.frame())

fwrite(effect_tbl, out_table, sep = "\t", na = "NA")
message("Wrote: ", out_table)
