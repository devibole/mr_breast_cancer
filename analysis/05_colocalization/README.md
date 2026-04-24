# Colocalization

This folder holds EUR colocalization analyses for candidate proteins.

## Added scripts

- `run_protein_coloc_analysis.R`
- `merge_coloc_results.R`
- `run_protein_coloc_swarm.sh`

## Manuscript links

- Methods: Colocalization Analysis
- Results: Multi-faceted Validation
- Table 2
- Supplementary Table 6

## Related source files

- `/Users/godboledd/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/haoyu/mr_bc/cispQTL/coloc/coloc_temp.R`
- `/Users/godboledd/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/haoyu/mr_bc/cispQTL/coloc/coloc/coloc_temp_JW.R`
- `/Users/godboledd/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/haoyu/mr_bc/cispQTL/coloc/coloc/coloc.abf_qtl_gwas.R`
- `/Users/godboledd/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/haoyu/mr_bc/coloc_bc_swarm.R`

## Notes

- `run_protein_coloc_analysis.R` is kept as a standalone script so this analysis can be rerun independently.
- The tarball provides a more explicit candidate-protein colocalization workflow than the older coloc scripts in the main workspace.
- The swarm submission helper was included because it is part of the original execution pattern.
- Hardcoded HPC paths are preserved for provenance and need to be replaced before use outside the original environment.
