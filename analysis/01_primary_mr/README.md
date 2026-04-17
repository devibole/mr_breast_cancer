# Primary Ancestry-Specific MR

This folder contains the ancestry-specific forward Mendelian randomization workflows used to pair UKB-PPP cis-pQTL instruments with breast cancer GWAS summary statistics.

## Included Workflows

- `eur/run_primary_mr_eur.R`
- `eas/run_primary_mr_eas.R`
- `afr/run_primary_mr_afr.R`

## Script Style

Each ancestry script is self-contained so it can be reviewed and run independently without relying on helper code from another file.

## Expected Inputs

Each ancestry script expects:

- UKB-PPP protein summary-statistic tar archives
- an rsID lookup RData file containing `all_rsids`
- Olink protein annotation TSV
- ancestry-specific breast cancer GWAS summary statistics
- a PLINK binary
- ancestry-matched or chosen LD reference files for clumping

## Expected Outputs

Each run writes:

- one MR-input file per protein
- one MR-result file per protein index

The exact export paths are controlled through the script-local `config` block.

## Reproducibility Note

These scripts preserve the original analysis logic, including:

- cis-window restriction to +/- 1 Mb around the coding gene,
- MAF and INFO filtering,
- LD clumping using PLINK,
- IVW, weighted median, MR-Egger, pIVW, and MR-RAPS estimates.

Before rerunning, replace environment-specific absolute paths in the `config` block.
