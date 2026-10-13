# =============================================================================
# Biomedical multi-scale allele-frequency analysis
#
# Recompute effect-allele frequencies for OGVP biomedical variants with snputils
# at two scales (panel a: 42 populations, Tuvalu split in two; panel b: 7 Group regions), add
# credible intervals (Jeffreys for populations; beta-binomial for regions, with
# populations as units), significance stars with BH-FDR, and report
# SNPs above the allele-frequency threshold (two-panel figure, CrI plots, tables).
#
# Run (from this folder):
#     snakemake --cores 1
#
# All inputs are READ-ONLY; every output stays under results/.
# =============================================================================
configfile: "config.yaml"

IN = config["inputs"]
GENO = IN["geno_prefix"]

rule all:
    input:
        "results/biomedical_samples.tsv",
        "results/freqs_final.tsv",
        "results/table_reported_snps.tsv",
        "results/table_reported_snps_summary.tsv",
        "results/two_panel_figure.pdf",
        "results/ci_two_scales.pdf",

# 00 - manifest of the exact samples used (815, from the genotype .fam); splits Tuvalu
#      into Tuvalu Micronesia / Tuvalu Polynesia
rule sample_manifest:
    input:
        fam        = GENO + ".fam",
        popinfo    = IN["popinfo"],
        clst_index = IN["clst_index"],
        split      = IN["tuvalu_split"],
        cfg        = "config.yaml",
    output:
        "results/biomedical_samples.tsv",
    shell:
        "Rscript scripts/00_sample_manifest.R {input.cfg}"

# 01 - define the biomedical effect-variant set + lock the effect allele
rule effect_variants:
    input:
        annovar   = IN["annovar_clean"],
        clinvar   = IN["clinvar"],
        ga_index  = IN["genomeasia_index"],
        bim       = GENO + ".bim",
        cfg       = "config.yaml",
    output:
        "results/effect_variant_table.tsv",
    shell:
        "Rscript scripts/01_define_effect_variants.R {input.cfg}"

# 02 - recompute effect-allele frequencies with snputils (both scales)
rule compute_freqs:
    input:
        eff  = "results/effect_variant_table.tsv",
        bed  = GENO + ".bed",
        bim  = GENO + ".bim",
        fam  = GENO + ".fam",
        samples = "results/biomedical_samples.tsv",
        cfg  = "config.yaml",
    output:
        clst  = "results/freq_clst.tsv",
        group = "results/freq_group.tsv",
        orient = "results/orientation_report.tsv",
    shell:
        "python scripts/02_compute_freqs_snputils.py {input.cfg}"

# 04 - Jeffreys credible intervals (allele-frequency uncertainty)
rule credible_intervals:
    input:
        clst  = "results/freq_clst.tsv",
        group = "results/freq_group.tsv",
        cfg   = "config.yaml",
    output:
        "results/freqs_with_ci.tsv",
    shell:
        "Rscript scripts/04_credible_intervals.R {input.cfg}"

# 04b - regional beta-binomial model (populations as units): regional intervals that
#       include between-population variation, and the regional test
rule regional_model:
    input:
        ci      = "results/freqs_with_ci.tsv",
        samples = "results/biomedical_samples.tsv",
        cfg     = "config.yaml",
    output:
        "results/regional_betabinomial.tsv",
    shell:
        "Rscript scripts/04b_regional_betabinomial.R {input.cfg}"

# 05 - tests: Fisher's exact (populations), beta-binomial LRT (regions); BH-FDR per
#      scale; star = p < 0.001 AND q < 0.05; reported SNPs and reporting tables
rule differentiation:
    input:
        ci  = "results/freqs_with_ci.tsv",
        bb  = "results/regional_betabinomial.tsv",
        eff = "results/effect_variant_table.tsv",
        clin = IN["clinical_allele_sources"],
        cfg = "config.yaml",
    output:
        final   = "results/freqs_final.tsv",
        table   = "results/table_reported_snps.tsv",
        summary = "results/table_reported_snps_summary.tsv",
    shell:
        "Rscript scripts/05_differentiation.R {input.cfg}"

# 06 - two-panel tile figure: reported SNPs, stars, gnomAD reference columns
rule figure:
    input:
        final  = "results/freqs_final.tsv",
        eff    = "results/effect_variant_table.tsv",
        gnomad = IN["gnomad"],
        cfg    = "config.yaml",
    output:
        "results/two_panel_figure.pdf",
    shell:
        "Rscript scripts/06_plot_two_panel.R {input.cfg}"

# 08 - credible-interval plots at both scales (reported SNPs)
rule ci_plots:
    input:
        final   = "results/freqs_final.tsv",
        samples = "results/biomedical_samples.tsv",
        cfg     = "config.yaml",
    output:
        main = "results/ci_two_scales.pdf",
    shell:
        "Rscript scripts/08_plot_credible_intervals.R {input.cfg}"
