# bc-protein-mr

Code repository for the manuscript:

`Cross-ancestry proteome-wide Mendelian randomization identifies 12 plasma protein candidates for breast cancer risk`

This repository is being assembled from the original analysis workspace. The first code family added here is the ancestry-specific forward MR workflow for:

- European ancestry (`EUR`)
- East Asian ancestry (`EAS`)
- African ancestry (`AFR`)

## Current Repository Layout

```text
bc-protein-mr/
├── docs/
├── analysis/
│   └── 01_primary_mr/
│       ├── afr/
│       ├── eas/
│       ├── eur/
│       └── README.md
├── analysis/02_cross_ancestry_meta/
├── analysis/03_subtype_specific_mr/
├── analysis/04_reverse_mr/
├── analysis/05_colocalization/
├── analysis/06_external_replication/
├── analysis/07_aou_pgs_validation/
├── analysis/08_sensitivity_analyses/
├── analysis/09_tables_and_summary_outputs/
├── figures/
└── .gitignore
```

## Important Path Note

Many scripts in this repository were developed in the original NIH/HPC project environment and therefore contain hardcoded absolute paths. These paths are preserved for provenance, but they must be updated before rerunning the code in a new local, cluster, or cloud environment.

In the cleaned scripts added here, those paths are grouped in a single `config` block near the top of each file.

## What Has Been Added So Far

- Annotated ancestry-specific MR scripts for `EUR`, `EAS`, and `AFR`
- Manuscript-based placeholder folders for all remaining analysis families
- A code-fetch checklist tying manuscript sections to likely source scripts in the current workspace
- A folder-level README describing inputs, outputs, and workflow expectations

## Notes

- The original analysis code remains untouched in the parent workspace.
- The public-facing repository contains cleaned copies intended for annotation and organization.
- Controlled-access datasets are not included in this repository.
