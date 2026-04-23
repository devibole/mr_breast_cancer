# mr_breast_cancer

Code repository for the manuscript:

`Cross-ancestry proteome-wide Mendelian randomization identifies 12 plasma protein candidates for breast cancer risk`

This repository is being assembled from the original analysis workspace. The first code family added here is the ancestry-specific forward MR workflow for:

- European ancestry (`EUR`)
- East Asian ancestry (`EAS`)
- African ancestry (`AFR`)

## Current Repository Layout

```text
mr_breast_cancer/
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

Many scripts in this repository were developed in the original NIH/HPC project environment and therefore contain hardcoded absolute paths. 

In the cleaned scripts added here, those paths are grouped in a single `config` block near the top of each file.

