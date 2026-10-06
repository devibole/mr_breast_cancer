# Cross-ancestry proteome-wide Mendelian randomization of breast cancer

Analysis code accompanying the manuscript on plasma proteins and breast cancer risk across ancestry groups.

## Analysis files

### Primary Mendelian randomization

- `analysis/main/run_primary_mr_eur.R` — MR using discovery UKB pQTL summary stats exposure with the European breast cancer GWAS summary stats as outcome data.
- `analysis/main/run_primary_mr_eas.R` — outcome data changed to East Asian breast cancer GWAS summary stats.
- `analysis/main/run_primary_mr_afr.R` — outcome data changed to African breast cancer GWAS summary stats.
- `analysis/main/compile_primary_results.R` — combines per-protein results into ancestry-level result files.
- `analysis/main/run_cross_ancestry_meta.R` — combines ancestry-specific MR estimates into cross-ancestry meta-analysis.

The primary MR scripts can be used with either the UKB-PPP discovery pQTL files or the UKB-PPP combined pQTL files by changing `MRBC_PROTEIN_ARCHIVE_DIR` and the output directories. Other input paths can also be supplied through the environment variables defined in each script.

### Follow-up analyses

- `analysis/run_subtype_mr_ivw.R` — evaluates candidate proteins across breast cancer subtypes.
- `analysis/run_reverse_mr.R` — tests breast cancer as the exposure and plasma protein levels as the outcome.
- `analysis/run_decode_mr.R` — parallels run_primary_mr_eur.R but using deCODE pQTL data.
- `analysis/run_eas_japan_top12_mr.R` — parallels run_primary_mr_eas.R but using Japanese pQTL data.
- `analysis/run_eas_bbj_bcac_sensitivity.R` — compares East Asian MR estimates using the BBJ, BCAC, and combined outcome data.

### Colocalization

- `analysis/coloc/run_protein_coloc_analysis.R` — performs regional protein–breast cancer colocalization.
- `analysis/coloc/merge_coloc_results.R` — combines the protein-specific colocalization results.

## Data requirements

The repository does not distribute the pQTL or breast cancer GWAS summary statistics, LD reference panels, annotations, or variant lookup files. Paths in the scripts must be updated for the user's computing environment.

The analyses use R, PLINK, and the R packages loaded within each script.
