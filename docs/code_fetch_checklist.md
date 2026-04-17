# Code Fetch Checklist

This checklist uses the manuscript text to identify all code families that still need to be added to the repository beyond the ancestry-specific primary MR scripts.

## Already Added

- `analysis/01_primary_mr/`
  - `eur/run_primary_mr_eur.R`
  - `eas/run_primary_mr_eas.R`
  - `afr/run_primary_mr_afr.R`

## Still To Fetch

| Analysis family | Manuscript evidence | Likely source files already in workspace | Status |
|---|---|---|---|
| Cross-ancestry fixed-effect meta-analysis | Primary discovery analysis; Figure 2; Figure 3; Table 2; Supplementary Tables 2-3 | `ancestries/meta.R`, `ancestries/meta_full.R`, `2026_updated/scripts/run_overlap_ivw_from_input_files.R`, `2026_updated/scripts/run_top12_overlap_ancestry_ivw.py` | Placeholder created |
| Subtype-specific MR | Figure 4; Supplementary Table 4 | `cispQTL/cis_subtypes_heterogeniety.R` | Placeholder created |
| Reverse-direction MR | Supplementary Table 5 | `bidirectional/mr_bd.R` | Placeholder created |
| Colocalization | Table 2; Supplementary Table 6 | `cispQTL/coloc/coloc_temp.R`, `cispQTL/coloc/coloc/coloc_temp_JW.R`, `cispQTL/coloc/coloc/coloc.abf_qtl_gwas.R`, `coloc_bc_swarm.R` | Placeholder created |
| deCODE replication | Table 2 | `2025_updated/decode_testing.R` | Placeholder created |
| JCTF replication | Table 2 | Not clearly identified yet in this workspace | Need to fetch |
| AoU PGS validation | Table 2; Supplementary Table 7 | Not clearly identified yet in this workspace; likely from Jacob Williams / AoU workbench code | Need to fetch |
| Sensitivity analyses: weighted median / pIVW / MR-RAPS / relaxed IV thresholds | Supplementary Figure 1; Supplementary Table 2 | partly embedded in ancestry scripts; possibly also `2026_updated/scripts/run_overlap_ivw_from_input_files.R` and `ancestries/meta_full.R` | Placeholder created |
| DNPH1 EAS-weight sensitivity analysis | Results and Methods | Not clearly isolated yet; may be embedded in one of the MR or sensitivity scripts | Need to fetch / identify |
| Supplementary table builders | Supplementary Tables 1-7 | `2026_updated/scripts/build_supp_table1_all_ivs.py`, `2026_updated/scripts/render_table_s1_like_excel.py`, `cispQTL/tables_and_misc.R` | Placeholder created |
| Figure scripts | Figures 2-4 and supplement | `2025_updated/figures/main_forest_volcano.R`; additional plotting code also embedded in `ancestries/meta.R` and `ancestries/meta_full.R` | Placeholder created |

## Likely Not Yet In Repo But Worth Checking

- AoU phenotype definition / OMOP concept extraction code
- AoU pooled and ancestry-stratified logistic regression code
- JCTF-specific instrument selection and MR script
- final figure export scripts for Figure 3 forest plot and Figure 4 heatmap if separate from analysis scripts
- any shell or swarm submission scripts needed for HPC reruns

## Suggested Next Fetch Order

1. Cross-ancestry meta-analysis
2. Subtype-specific MR
3. Reverse MR
4. Colocalization
5. deCODE replication
6. JCTF replication
7. AoU PGS validation
8. Supplementary table / figure assembly

## Likely People To Ask For Missing Code

- AoU PGS validation: Jacob Williams
- JCTF replication script if not local: Xueyao / Haoyu
- final figure assembly if not local: whichever version generated the submitted figures

