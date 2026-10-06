# Cross-ancestry proteome-wide Mendelian randomization of breast cancer

This repository contains the analysis code accompanying the manuscript:

> **Cross-ancestry proteome-wide Mendelian randomization identifies 12 plasma protein candidates for breast cancer risk**

The study uses cis-protein quantitative trait loci (cis-pQTLs) as genetic instruments to test whether genetically predicted plasma protein levels are associated with breast cancer risk. Analyses are run separately with European (EUR), East Asian (EAS), and African ancestry (AFR) breast cancer GWAS data and then combined in a cross-ancestry meta-analysis. Candidate proteins are evaluated using subtype-specific MR, reverse-direction MR, colocalization, and independent pQTL replication.

## Repository scope

This is a manuscript-analysis repository, not a complete copy of the working project directory. It includes code that generates scientific results reported in the main text or Supplementary Tables 1–7. It intentionally excludes scripts used only to assemble, format, reorder, or export the final supplementary workbook.

The repository also excludes exploratory diagnostics and one-off analyses that are not reported in the current manuscript or supplement. In particular:

- the separate DNPH1-only Japanese pQTL script is omitted because DNPH1 is already analyzed by the full Japanese top-12 replication script;
- shared-instrument/overlapping-IV analyses are omitted because they are not reported in the current manuscript or Supplementary Tables 1–7; and
- later UK Biobank 10,000-participant LD-reference diagnostics are not the analysis reported in Supplementary Table 7, which uses the 1000 Genomes EUR LD reference.

## Analysis overview

```text
UKB-PPP discovery cis-pQTLs
        |
        +--> EUR forward MR --+
        +--> EAS forward MR --+--> fixed-effect cross-ancestry meta-analysis
        +--> AFR forward MR --+                  |
                                                 +--> subtype-specific MR
                                                 +--> reverse-direction MR
                                                 +--> colocalization
                                                 +--> deCODE and Japanese replication

UKB-PPP combined/full-cohort cis-pQTLs
        |
        +--> EUR/EAS/AFR forward MR with the same 1000 Genomes EUR LD reference
                                      |
                                      +--> cross-ancestry IVW sensitivity analysis
                                           (Supplementary Table 7)
```

## Mapping from results to code

| Manuscript result | Analysis code | What the code produces |
|---|---|---|
| Supplementary Table 1: cis-pQTL instruments and ancestry-specific variant associations | `analysis/main/run_primary_mr_eur.R`, `run_primary_mr_eas.R`, `run_primary_mr_afr.R` | Harmonized, clumped instrument-level files for each protein and ancestry |
| Supplementary Table 2: primary and sensitivity MR results | The three primary MR scripts plus `analysis/main/run_cross_ancestry_meta.R` | Ancestry-specific IVW and sensitivity estimates, followed by cross-ancestry meta-analysis |
| Supplementary Table 3: ancestry-specific results for the leading proteins | The three primary MR scripts | EUR, EAS, and AFR protein-level MR estimates |
| Supplementary Table 4: intrinsic-like breast cancer subtypes | `analysis/run_subtype_mr_ivw.R` | IVW estimates for five subtypes and cross-subtype heterogeneity statistics |
| Supplementary Table 5: reverse-direction MR | `analysis/run_reverse_mr.R` | Effect of genetic liability to breast cancer on candidate protein levels |
| Supplementary Table 6: colocalization | `analysis/coloc/run_protein_coloc_analysis.R` | Regional H0-H4 posterior probabilities, SNP posterior probabilities, and credible sets |
| Supplementary Table 7: combined/full-cohort pQTL sensitivity analysis | The three configurable primary MR scripts, `compile_primary_results.R`, and `run_cross_ancestry_meta.R` | Cross-ancestry IVW results using UKB-PPP Combined pQTL data and 1000 Genomes EUR LD |
| Main Table 2: independent European pQTL replication | `analysis/run_decode_mr.R` | deCODE-based MR estimates for available candidate proteins |
| Main Table 2: independent East Asian pQTL replication | `analysis/run_eas_japan_top12_mr.R` | Japanese pQTL-based MR estimates for the 12 candidates |

`analysis/main/compile_primary_results.R` only combines per-protein MR result files into the ancestry-level input needed by the meta-analysis. It does not create or format a supplementary table.

## Repository layout

```text
analysis/
├── main/
│   ├── run_primary_mr_eur.R
│   ├── run_primary_mr_eas.R
│   ├── run_primary_mr_afr.R
│   ├── compile_primary_results.R
│   └── run_cross_ancestry_meta.R
├── coloc/
│   ├── run_protein_coloc_analysis.R
│   ├── run_protein_coloc_swarm.sh
│   └── merge_coloc_results.R
├── run_subtype_mr_ivw.R
├── run_reverse_mr.R
├── run_decode_mr.R
└── run_eas_japan_top12_mr.R
```

## Primary ancestry-specific MR

The EUR, EAS, and AFR scripts use the same analysis structure but read different breast cancer outcome GWAS files. Each script accepts a one-based protein archive index and performs the following steps:

1. Select one UKB-PPP protein summary-statistic archive from a sorted archive list.
2. Extract and deduplicate the protein summary statistics.
3. Map variants to rsIDs and GRCh37 positions.
4. Define cis variants within 1 Mb of the protein-coding gene.
5. Match the pQTL data to the ancestry-specific breast cancer GWAS and align effect alleles.
6. Exclude variants with minor allele frequency below 0.01 or imputation INFO at or below 0.3.
7. Select pQTLs at `P <= 5e-8`, `5e-7`, and `5e-6` and LD-clump them with PLINK (`r2 = 0.01`, 2,000 kb window).
8. Run IVW MR and, where the required number of instruments is available, weighted-median, MR-Egger, penalized IVW, and MR-RAPS analyses.
9. Write the harmonized instrument-level data and a protein-level result file.

The primary manuscript analysis uses UKB-PPP European discovery pQTL estimates as the exposure source for all three outcome ancestries. The checked-in default paths reproduce that configuration. All important paths can also be set through environment variables, which permits the same code to run the combined/full-cohort pQTL sensitivity analysis without maintaining a second copy of the analysis logic.

Example for the first discovery-pQTL archive:

```bash
Rscript analysis/main/run_primary_mr_eur.R 1
Rscript analysis/main/run_primary_mr_eas.R 1
Rscript analysis/main/run_primary_mr_afr.R 1
```

One protein should be processed per cluster-array task. Because the scripts now sort the archive paths, a given directory has a stable index-to-protein mapping.

## Combined/full-cohort pQTL sensitivity analysis

Supplementary Table 7 repeats the ancestry-specific IVW analysis using the UKB-PPP `Combined` pQTL archives while retaining the 1000 Genomes GRCh37 EUR LD reference used in the primary analysis. Set a separate input and result directory for each ancestry so the sensitivity outputs cannot overwrite the discovery analysis.

Example EUR task:

```bash
MRBC_PROTEIN_ARCHIVE_DIR=/path/to/UKB-PPP/Combined \
MRBC_LD_REFERENCE_DIR=/path/to/1000G/GRCh37/EUR \
MRBC_P_THRESHOLDS=5e-8 \
MRBC_INPUT_EXPORT_DIR=/path/to/combined/input_files/eur \
MRBC_RESULTS_EXPORT_DIR=/path/to/combined/results/eur \
Rscript analysis/main/run_primary_mr_eur.R 1
```

Use the analogous command with `run_primary_mr_eas.R` and `run_primary_mr_afr.R`, changing the two output directories for each ancestry. The outcome paths remain ancestry specific through each script's default configuration, or they can be overridden with `MRBC_OUTCOME_PATH`.

After all array tasks finish, compile the per-protein results:

```bash
Rscript analysis/main/compile_primary_results.R \
  /path/to/combined/results/eur /path/to/combined/eur.csv

Rscript analysis/main/compile_primary_results.R \
  /path/to/combined/results/eas /path/to/combined/eas.csv

Rscript analysis/main/compile_primary_results.R \
  /path/to/combined/results/afr /path/to/combined/afr.csv
```

Then run the same fixed-effect cross-ancestry meta-analysis against the combined-pQTL results:

```bash
MRBC_EUR_RESULTS_PATH=/path/to/combined/eur.csv \
MRBC_EAS_RESULTS_PATH=/path/to/combined/eas.csv \
MRBC_AFR_RESULTS_PATH=/path/to/combined/afr.csv \
MRBC_META_OUTPUT_PATH=/path/to/combined/meta.csv \
Rscript analysis/main/run_cross_ancestry_meta.R
```

## Cross-ancestry meta-analysis

`analysis/main/run_cross_ancestry_meta.R` joins the ancestry-level result tables by protein. For each MR method and pQTL threshold represented in the input files, it performs fixed-effect meta-analysis when at least two ancestry-specific estimates are available.

The output contains:

- the pooled log-odds estimate, standard error, 95% confidence interval, and P value;
- the number of contributing ancestry groups;
- Cochran's Q, heterogeneity P value, I-squared, and tau-squared;
- an aggregated Cauchy association test P value across available MR results; and
- Bonferroni- and Benjamini-Hochberg-adjusted P values.

The manuscript's primary estimate is the fixed-effect IVW result at `5e-8`. Other columns are sensitivity results and should not be substituted for the prespecified primary result.

## Candidate-protein follow-up analyses

### Breast cancer subtypes

`analysis/run_subtype_mr_ivw.R` applies the candidate cis-pQTL instruments to five BCAC intrinsic-like subtype GWAS datasets: luminal A-like, luminal B/HER2-negative-like, luminal B-like, HER2-enriched-like, and triple-negative breast cancer. It reports IVW estimates for each subtype and a Cochran Q test comparing effects across subtypes.

### Reverse-direction MR

`analysis/run_reverse_mr.R` reverses the primary analysis. Genome-wide significant BCAC variants are instruments for breast cancer liability, and UKB-PPP protein abundance is the outcome. The analysis tests whether the forward protein-to-breast-cancer associations could instead be explained by an effect of breast cancer liability on circulating protein levels.

### Colocalization

`analysis/coloc/run_protein_coloc_analysis.R` evaluates +/-150 kb regions around candidate cis-pQTL signals with `coloc.abf`. The protein trait is modeled as quantitative and breast cancer as case-control. For every tested region, the script reports posterior probabilities for H0-H4, identifies the strongest SNP-level posterior probabilities, and records credible-set membership.

The swarm file documents the protein-index mapping used for the study run. `merge_coloc_results.R` combines the per-protein colocalization workbooks into study-level result files; it does not generate the final Supplementary Table 6 layout.

### Independent pQTL replication

`analysis/run_decode_mr.R` repeats forward MR with deCODE pQTL statistics and the EUR breast cancer GWAS. It applies the deCODE exclusion list, joins annotation and allele-frequency information, converts unmapped GRCh38 positions to GRCh37 where required, and evaluates multiple pQTL thresholds.

`analysis/run_eas_japan_top12_mr.R` evaluates all 12 candidates using Japanese pQTL statistics and the EAS breast cancer outcome. It reads the full pQTL file once, selects each candidate by gene name, harmonizes variants, performs EAS LD clumping, and reports the available MR estimates. DNPH1 is included in this candidate list, so a separate DNPH1-only replication script is unnecessary.

## Configuration and required data

The repository does not redistribute controlled or study-specific pQTL/GWAS summary statistics, individual-level data, LD panels, the rsID map, or the GRCh38-to-GRCh37 chain file. Users must obtain authorized access and update the configuration paths.

The main ancestry-specific scripts recognize these environment variables:

| Variable | Purpose |
|---|---|
| `MRBC_PROTEIN_ARCHIVE_DIR` | Discovery or Combined UKB-PPP protein archive directory |
| `MRBC_OLINK_MAP_PATH` | Olink protein/gene annotation file |
| `MRBC_RSID_LOOKUP_RDATA` | GRCh37 variant-to-rsID lookup containing `all_rsids` |
| `MRBC_OUTCOME_PATH` | Ancestry-specific breast cancer GWAS file |
| `MRBC_LD_REFERENCE_DIR` | PLINK LD-reference directory |
| `MRBC_PLINK_BINARY` | PLINK 1.9 executable |
| `MRBC_INPUT_EXPORT_DIR` | Instrument-level output directory |
| `MRBC_RESULTS_EXPORT_DIR` | Per-protein MR result directory |
| `MRBC_SCRATCH_ROOT` | Temporary working directory |
| `MRBC_P_THRESHOLDS` | Comma-separated pQTL thresholds; defaults to `5e-8,5e-7,5e-6` |

The cross-ancestry script recognizes `MRBC_EUR_RESULTS_PATH`, `MRBC_EAS_RESULTS_PATH`, `MRBC_AFR_RESULTS_PATH`, and `MRBC_META_OUTPUT_PATH`.

R packages used across the analyses include `ACAT`, `MendelianRandomization`, `biomaRt`, `coloc`, `data.table`, `dplyr`, `meta`, `mr.pivw`, `mr.raps`, `readxl`, `rtracklayer`, `tidyverse`, `vroom`, and `writexl`. PLINK 1.9 is used for the reported 1000 Genomes LD-clumping workflows.

## Reproduction order

1. Obtain the required pQTL, breast cancer GWAS, annotation, rsID, chain, and LD-reference data.
2. Run the EUR, EAS, and AFR primary scripts across all discovery protein archives.
3. Compile each ancestry's per-protein outputs with `compile_primary_results.R`.
4. Run `run_cross_ancestry_meta.R` and identify the candidate proteins from the fixed-effect IVW result at `5e-8`.
5. Run subtype MR, reverse MR, colocalization, and deCODE/Japanese replication for the candidate set.
6. Repeat steps 2-4 with the UKB-PPP Combined archive directory to reproduce Supplementary Table 7.
7. Compare analysis outputs, not workbook formatting, with the corresponding manuscript and supplementary results.

## Manuscript-code coverage note

The manuscript draft available with this project also reports polygenic-score analyses in the All of Us cohort. That analysis code was not present in the local project materials used to prepare this repository and therefore is not included here. If the All of Us results remain in the submitted manuscript, the corresponding analysis code should be obtained from the analyst who ran it and added as a separate, clearly documented workflow.

Because the source summary statistics and individual-level data are not distributed here, this repository is a transparent record of the analysis logic rather than a one-command reproduction package.
