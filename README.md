# mr_breast_cancer

Code repository for the manuscript:

`Cross-ancestry proteome-wide Mendelian randomization identifies 12 plasma protein candidates for breast cancer risk`

## Included Analysis Scripts

- `run_primary_mr_eur.R`
- `run_primary_mr_eas.R`
- `run_primary_mr_afr.R`
- `run_cross_ancestry_meta.R`
- `run_subtype_mr_ivw.R`
- `run_reverse_mr.R`
- `run_protein_coloc_analysis.R`
- `merge_coloc_results.R`
- `run_protein_coloc_swarm.sh`
- `run_decode_mr.R`
- `run_eas_japan_top12_mr.R`
- `run_dnph1_eas_mr.R`
- `eur_overlapping_ivs.R`
- `eas_overlapping_ivs.R`
- `afr_overlapping_ivs.R`

## Important Path Note

Many scripts in this repository were developed in the original NIH/HPC project environment and therefore contain hardcoded absolute paths. These paths are preserved here, but would need to be updated before rerunning the code in a new local, cluster, or cloud environment.
