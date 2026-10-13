# Allele frequencies of biomedical variants

## Inputs
- `data/`: provided separately.
- `resources/clinical_allele_sources.tsv`: clinical allele and its source for each reported variant.
- `config.yaml`: paths and thresholds.

## Pipeline steps
- `00`: sample manifest (population and region of each sample)
- `01`: variant table with the allele with the clinical effect
- `02`: copies of the clinical allele (c) and called chromosomes (n) per population and region, with snputils
- `04`: Jeffreys 95% credible intervals
- `04b`: beta-binomial regional intervals and tests
- `05`: p-values, FDR, stars and the variant tables
- `06`: two-panel figure
- `08`: interval figure

## Statistics
- Frequency = c / n of the clinical allele.
- Intervals: 95% Jeffreys credible interval for populations; 95% profile-likelihood confidence interval from a beta-binomial model (each population an independent observation) for regions with several populations; Jeffreys interval for Taiwan and Madagascar (one population each).
- Tests: each group against all other groups at the same scale. Two-sided Fisher's exact test for populations and for Taiwan and Madagascar; beta-binomial likelihood-ratio test for the other regions.
- Multiple testing: Benjamini–Hochberg FDR, separately for populations and regions.
- Star: p < 0.001 and q < 0.05.
- Reported variants: clinical-allele frequency above 0.1 in at least one population (`report_min_af`); variants listed under `report_exclude` in `config.yaml` are left out.

## Outputs (`results/`)
- `freqs_final.tsv`: all 47 variants × 42 populations and 7 regions: counts, frequency, interval, test, p, q, star
- `table_reported_snps.tsv`: the 19 reported variants, with figure label, clinical allele and source
- `table_reported_snps_summary.tsv`: one row per reported variant, with the groups that have stars
- `table_excluded_snps.tsv`: the 28 variants not reported, with the reason for each
- `two_panel_figure.pdf`, `ci_two_scales.pdf`
