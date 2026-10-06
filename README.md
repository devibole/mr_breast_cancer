# Cross-ancestry proteome-wide Mendelian randomization of breast cancer

This repository contains the analysis code for the manuscript:

> **Cross-ancestry proteome-wide Mendelian randomization identifies 12 plasma protein candidates for breast cancer risk**

The study uses protein quantitative trait loci (pQTLs) as genetic instruments to test whether genetically predicted plasma protein levels are associated with breast cancer risk. The workflow begins with ancestry-specific forward Mendelian randomization (MR), combines evidence across ancestries, and then evaluates the leading proteins using subtype, reverse-MR, colocalization, replication, and shared-instrument sensitivity analyses.

## How the code maps to the analysis

The scripts are organized by their role in the manuscript rather than as a single executable pipeline:

1. **Primary ancestry-specific MR:** run the proteome-wide analysis separately using European (EUR), East Asian (EAS), and African ancestry (AFR) breast cancer GWAS outcome data.
2. **Cross-ancestry meta-analysis:** combine the ancestry-specific MR estimates and calculate heterogeneity and multiple-testing statistics.
3. **Follow-up analyses:** assess breast cancer subtypes, reverse causation, and replication in independent pQTL resources.
4. **Colocalization:** test whether protein abundance and breast cancer risk are consistent with a shared causal variant in the same locus.
5. **Overlapping-instrument sensitivity analysis:** repeat MR using only instruments shared across EUR, EAS, and AFR analyses.

```text
UKB-PPP protein GWAS
        |
        +--> EUR forward MR --+
        +--> EAS forward MR --+--> cross-ancestry meta-analysis --> candidate proteins
        +--> AFR forward MR --+                                  |
                                                                  +--> subtype MR
                                                                  +--> reverse MR
                                                                  +--> colocalization
                                                                  +--> deCODE/Japanese replication
                                                                  +--> shared-IV sensitivity analysis
```

## Repository layout

```text
analysis/
├── main/                 Primary EUR/EAS/AFR MR and cross-ancestry meta-analysis
├── coloc/                Protein–breast cancer colocalization and result merging
├── overlapping_ivs/      Sensitivity analyses restricted to shared instruments
├── run_subtype_mr_ivw.R  MR for five breast cancer molecular subtypes
├── run_reverse_mr.R      Reverse-direction MR: breast cancer liability → protein
├── run_decode_mr.R       Independent replication using deCODE pQTL data
├── run_eas_japan_top12_mr.R
│                         Replication of the 12 candidates using Japanese pQTL data
└── run_dnph1_eas_mr.R    Focused Japanese pQTL analysis of DNPH1
```

## Script-by-script guide

### 1. Primary ancestry-specific MR

The three scripts below share the same core workflow. Each accepts a **one-based protein archive index** as its command-line argument, selects that UKB-PPP protein archive, and:

- extracts and deduplicates the protein summary statistics;
- maps variants to rsIDs and GRCh37 positions;
- restricts instruments to the protein's cis region (±1 Mb around the encoding gene);
- merges the protein and ancestry-specific breast cancer summary statistics;
- aligns effect alleles and filters variants on minor allele frequency and imputation quality;
- selects protein-associated variants and performs PLINK LD clumping (`r2 = 0.01`, 2,000 kb window);
- runs inverse-variance weighted (IVW), weighted-median, MR-Egger, PIVW, and MR-RAPS analyses where the method's requirements are met; and
- writes both the harmonized instrument-level data and a protein-level MR result file.

| Script | Breast cancer outcome data | Main output role |
|---|---|---|
| [`analysis/main/run_primary_mr_eur.R`](analysis/main/run_primary_mr_eur.R) | BCAC European ancestry GWAS | Primary EUR MR estimates |
| [`analysis/main/run_primary_mr_eas.R`](analysis/main/run_primary_mr_eas.R) | BCAC–BBJ East Asian meta-analysis | Primary EAS MR estimates |
| [`analysis/main/run_primary_mr_afr.R`](analysis/main/run_primary_mr_afr.R) | African-ancestry breast cancer GWAS | Primary AFR MR estimates |

The current configuration uses a genome-wide pQTL threshold of `5e-8`. The result rows retain columns for several MR estimators so that estimates, standard errors, confidence intervals, P values, instrument counts, and instrument identifiers can be carried into the downstream summaries.

Example for the first protein archive:

```bash
Rscript analysis/main/run_primary_mr_eur.R 1
Rscript analysis/main/run_primary_mr_eas.R 1
Rscript analysis/main/run_primary_mr_afr.R 1
```

On a cluster, the index can be supplied by an array job so that one protein is processed per task. The exact protein corresponding to an index is determined by the order returned by `list.files()` in the configured protein archive directory; preserve or record that ordering when reproducing a full run.

### 2. Cross-ancestry meta-analysis

[`analysis/main/run_cross_ancestry_meta.R`](analysis/main/run_cross_ancestry_meta.R) reads the consolidated EUR, EAS, and AFR result tables and joins them by protein. For every protein, MR method, and available pQTL threshold, it performs a fixed-effect meta-analysis when at least two ancestry-specific estimates are available.

The output includes:

- the pooled effect estimate, standard error, 95% confidence interval, and P value;
- the number of contributing ancestry groups;
- Cochran's Q, its P value, I-squared, and tau-squared;
- an aggregated Cauchy association test (ACAT) P value across available MR results; and
- Bonferroni- and Benjamini–Hochberg-adjusted P values.

The script expects already consolidated ancestry-level files (`eur.csv`, `eas.csv`, and `afr.csv`). Combining the individual `result_<index>.txt` files into those three tables was performed outside this repository and is therefore a required preprocessing step before running the meta-analysis script.

### 3. Breast cancer subtype MR

[`analysis/run_subtype_mr_ivw.R`](analysis/run_subtype_mr_ivw.R) repeats cis-pQTL MR using BCAC subtype-specific outcome estimates. It evaluates five molecular subtypes:

- Luminal A;
- Luminal B;
- HER2-enriched;
- Luminal B/HER2-negative; and
- triple-negative breast cancer.

For each protein archive index, the script harmonizes and clumps the instruments, runs IVW MR separately for each available subtype, and tests heterogeneity among the subtype-specific causal estimates. Its result file contains the subtype estimates and P values together with the cross-subtype Q statistic.

```bash
Rscript analysis/run_subtype_mr_ivw.R 1
```

### 4. Reverse-direction MR

[`analysis/run_reverse_mr.R`](analysis/run_reverse_mr.R) reverses the exposure and outcome used in the primary analysis. Genome-wide significant BCAC variants are treated as instruments for breast cancer liability, and UKB-PPP protein abundance is treated as the outcome. The script processes one protein archive per index, harmonizes and clumps breast cancer instruments, runs IVW MR, and saves both the instrument-level input and protein-level result.

This analysis addresses whether the forward associations might instead reflect an effect of genetic liability to breast cancer on circulating protein levels.

```bash
Rscript analysis/run_reverse_mr.R 1
```

### 5. Independent pQTL replication

#### deCODE

[`analysis/run_decode_mr.R`](analysis/run_decode_mr.R) repeats the forward MR using deCODE protein summary statistics. In addition to the common harmonization and cis-instrument steps, it removes deCODE-designated excluded variants, merges the deCODE annotation and allele-frequency data, and lifts unmapped GRCh38 coordinates to GRCh37 before matching them to BCAC. It tests pQTL thresholds of `5e-8`, `5e-7`, and `5e-6` and reports IVW, weighted-median, MR-Egger, PIVW, and MR-RAPS results as available.

The command-line argument is the one-based index of a deCODE summary-statistic file:

```bash
Rscript analysis/run_decode_mr.R 1
```

#### Japanese pQTL data

[`analysis/run_eas_japan_top12_mr.R`](analysis/run_eas_japan_top12_mr.R) uses Japanese pQTL summary statistics to evaluate the 12 leading proteins against the EAS breast cancer outcome data. It loops over the predefined candidate list, harmonizes and clumps variants using the EAS LD reference, and runs MR at `5e-8`, `5e-7`, and `5e-6`.

[`analysis/run_dnph1_eas_mr.R`](analysis/run_dnph1_eas_mr.R) is the focused version of this workflow for **DNPH1**. It extracts DNPH1 by Ensembl gene ID and saves the threshold-specific instrument sets and MR estimates.

These two scripts contain their target protein definitions internally and therefore do not require a protein-index argument:

```bash
Rscript analysis/run_eas_japan_top12_mr.R
Rscript analysis/run_dnph1_eas_mr.R
```

### 6. Colocalization

[`analysis/coloc/run_protein_coloc_analysis.R`](analysis/coloc/run_protein_coloc_analysis.R) evaluates the candidate protein loci using approximate-Bayes-factor colocalization. The script:

- accepts a one-based index for a predefined candidate protein;
- extracts protein and BCAC summary statistics around each target signal (the configured window is ±150 kb);
- runs `coloc.abf` using a quantitative model for the protein and a case-control model for breast cancer;
- records posterior probabilities, including H3 (distinct causal variants) and H4 (a shared causal variant);
- identifies high-posterior variants and credible-set membership; and
- writes a multi-sheet Excel workbook plus regional diagnostic plots.

```bash
Rscript analysis/coloc/run_protein_coloc_analysis.R 1
```

[`analysis/coloc/run_protein_coloc_swarm.sh`](analysis/coloc/run_protein_coloc_swarm.sh) contains the original 16 protein-index calls and redirects each run's standard output and errors to separate log files. It documents the index-to-protein mapping used for the manuscript run; update its script and log paths before reuse.

After all per-protein colocalization jobs finish, [`analysis/coloc/merge_coloc_results.R`](analysis/coloc/merge_coloc_results.R) validates and combines the `All Results` and `Highly Significant SNPs` worksheets from the individual workbooks into two study-level Excel files.

### 7. Shared-instrument sensitivity analysis

The scripts in [`analysis/overlapping_ivs/`](analysis/overlapping_ivs/) identify the intersection of the EUR, EAS, and AFR instrument lists for each leading protein at `5e-8` and `5e-6`. They then rerun ancestry-specific IVW MR using exactly those shared variants. This separates differences caused by ancestry-specific instrument availability from differences in the variant–outcome associations.

| Script | pQTL and outcome analysis |
|---|---|
| [`eur_overlapping_ivs.R`](analysis/overlapping_ivs/eur_overlapping_ivs.R) | EUR shared-IV sensitivity analysis |
| [`eas_overlapping_ivs.R`](analysis/overlapping_ivs/eas_overlapping_ivs.R) | EAS shared-IV sensitivity analysis |
| [`afr_overlapping_ivs.R`](analysis/overlapping_ivs/afr_overlapping_ivs.R) | AFR shared-IV sensitivity analysis |

Each script writes a protein-level MR summary and a per-SNP table, allowing readers to see both the pooled estimate and the contribution of each shared instrument.

## Data and software requirements

The repository contains analysis code but does **not** redistribute the controlled or study-specific GWAS/pQTL summary statistics, LD reference panels, rsID map, or chain file. To rerun the analyses, users need authorized access to the corresponding resources and must update the paths in each script's configuration block.

The scripts use R and the following packages across the workflow:

```r
ACAT
MendelianRandomization
biomaRt
coloc
data.table
fs
meta
mr.pivw
mr.raps
readxl
rtracklayer
tidyverse
vroom
writexl
```

External requirements include:

- **PLINK 1.9** and the binary genotype reference files (`.bed`, `.bim`, and `.fam`) expected by each script for LD clumping;
- a GRCh38-to-GRCh37 chain file for the deCODE coordinate conversion;
- enough temporary disk space to extract per-protein pQTL archives; and
- network access to Ensembl GRCh37 through `biomaRt` when the local Olink annotation does not provide a gene region.

No R package lockfile is currently included, so users seeking an exact computational reproduction should record the R and package versions used in their environment.

## Paths and cluster assumptions

These scripts were developed in the original NIH/HPC project environment. The hardcoded absolute paths are intentionally retained as provenance, but they will not work unchanged elsewhere. Before running a script, review its path/configuration section and replace paths for:

- protein and breast cancer summary statistics;
- the Olink protein annotation and rsID lookup object (`all_rsids`);
- ancestry-specific LD reference panels;
- the PLINK executable;
- temporary/scratch storage; and
- result and instrument-output directories.

Several scripts detect `SLURM_JOB_ID` and use `/lscratch` during cluster runs. Most have a local `tempdir()` fallback, but [`analysis/run_dnph1_eas_mr.R`](analysis/run_dnph1_eas_mr.R) retains the original `/lscratch` assumption and must be edited for a non-SLURM environment.

Genome-build and LD-reference choices must be checked script by script. The primary MR and colocalization workflows match on GRCh37 coordinates, while the Japanese pQTL scripts parse the supplied `variant_id_hg19` field but currently point to an EAS LD directory labelled GRCh38. In addition, the checked-in EUR, EAS, and AFR primary MR configurations all point to the same EUR LD reference directory. These settings document the analysis environment as coded; confirm that the coordinate build and reference population are appropriate before attempting a new run.

## Recommended reproduction order

1. Obtain the required summary statistics, annotations, variant map, and LD reference panels.
2. Edit the `config` or path block in each script for the new environment.
3. Run the three primary MR scripts across every protein archive index.
4. Combine the per-protein results into one EUR, one EAS, and one AFR table with the column names expected by `run_cross_ancestry_meta.R`.
5. Run the cross-ancestry meta-analysis and identify the candidate set.
6. Run subtype MR, reverse MR, independent pQTL replication, colocalization, and overlapping-IV sensitivity analyses.
7. Compare the generated result tables with the corresponding manuscript and supplementary tables.

Because the underlying summary statistics are not included and the scripts retain environment-specific paths, this repository should currently be regarded as a transparent record of the manuscript analyses rather than a one-command reproducible package.
