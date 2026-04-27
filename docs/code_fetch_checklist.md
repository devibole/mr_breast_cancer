# Code Fetch Checklist

This checklist is aligned to the current manuscript (`Main_Text_Final.docx`) and tracks which analysis/code families should exist in the GitHub repository.

## Already Added

- `analysis/01_primary_mr/`
  - `eur/run_primary_mr_eur.R`
  - `eas/run_primary_mr_eas.R`
  - `afr/run_primary_mr_afr.R`
- `analysis/02_cross_ancestry_meta/run_cross_ancestry_meta.R`
- `analysis/03_subtype_specific_mr/run_subtype_mr_ivw.R`
- `analysis/04_reverse_mr/run_reverse_mr.R`
- `analysis/05_colocalization/`
  - `run_protein_coloc_analysis.R`
  - `merge_coloc_results.R`
  - `run_protein_coloc_swarm.sh`
- `analysis/06_external_replication/decode_eur/run_decode_mr.R`
- `analysis/06_external_replication/jctf_eas/run_eas_japan_top12_mr.R`
- `analysis/08_sensitivity_analyses/dnph1_eas_weights/run_dnph1_eas_mr.R`
- `analysis/08_sensitivity_analyses/overlapping_ivs/`
  - `eur_overlapping_ivs.R`
  - `eas_overlapping_ivs.R`
  - `afr_overlapping_ivs.R`

## Manuscript-to-Repo Map

| Analysis family | Manuscript evidence | Repo location | Status |
|---|---|---|---|
| Ancestry-specific primary MR (EUR, EAS, AFR) | Methods: ancestry-specific two-sample MR; Results: ancestry-specific associations; Supplementary Tables 1-3 | `analysis/01_primary_mr/` | Added |
| Cross-ancestry fixed-effect meta-analysis | Methods: fixed-effects meta-analysis; Results: primary discovery; Figures 2-3; Table 2; Supplementary Tables 2-3 | `analysis/02_cross_ancestry_meta/` | Added, but figure/table export helpers may still be separate |
| Sensitivity MR methods and alternative IV thresholds | Methods: weighted median, pIVW, MR-RAPS, relaxed thresholds; Results: Supplementary Figure 1; Supplementary Table 2 | `analysis/08_sensitivity_analyses/` plus primary MR scripts | Partly present; may need additional summary/export scripts |
| DNPH1 EAS-weight sensitivity analysis | Methods and Results: targeted DNPH1 sensitivity analysis | `analysis/08_sensitivity_analyses/dnph1_eas_weights/` | Added |
| Overlapping-IV ancestry sensitivity analyses | Methods-related sensitivity support | `analysis/08_sensitivity_analyses/overlapping_ivs/` | Added |
| Subtype-specific MR | Methods and Results: intrinsic-like subtype analyses; Figure 4; Supplementary Table 4 | `analysis/03_subtype_specific_mr/` | Added |
| Reverse-direction MR | Methods and Results: reverse MR; Supplementary Table 5 | `analysis/04_reverse_mr/` | Added |
| Colocalization | Methods and Results: colocalization; Table 2; Supplementary Table 6 | `analysis/05_colocalization/` | Added |
| External replication in deCODE | Methods and Results: deCODE replication; Table 2 | `analysis/06_external_replication/decode_eur/` | Added locally; still needs cleanup/push sequence |
| External replication in JCTF | Methods and Results: JCTF replication; Table 2 | `analysis/06_external_replication/jctf_eas/` | Added locally; still needs cleanup/push sequence |
| AoU PGS validation | Methods and Results: AoU phenotype-linked validation; Table 2; Supplementary Table 7 | `analysis/07_aou_pgs_validation/` | Placeholder only; code still missing from current workspace |
| Supplementary table assembly | Supplementary Tables 1-7 | `analysis/09_tables_and_summary_outputs/` | Placeholder only |
| Figure assembly / export | Figure 1 graphical abstract; Figures 2-4; supplementary figure(s) | likely `analysis/09_tables_and_summary_outputs/` or separate figure scripts | Only partly identified |

## Missing Or Still Unclear

- AoU PGS construction code
- AoU breast cancer phenotype definition / OMOP concept extraction code
- AoU pooled multi-ancestry logistic regression code
- AoU ancestry-stratified logistic regression code
- final table-building scripts for Supplementary Tables 1-7
- final figure-export scripts for Figures 1-4 and Supplementary Figure 1

## Best Current Guess For Missing Sources

- AoU PGS validation likely lives outside this workspace, probably from Jacob Williams / AoU workbench code
- figure-export code may be separate from the analysis scripts, especially for:
  - Figure 1 graphical abstract
  - Figure 3 forest plot
  - Figure 4 heatmap
  - Supplementary Figure 1 sensitivity summary

## Recommended Next Upload Order

1. `analysis/05_colocalization/run_protein_coloc_swarm.sh`
2. `analysis/06_external_replication/decode_eur/run_decode_mr.R`
3. `analysis/06_external_replication/jctf_eas/run_eas_japan_top12_mr.R`
4. `analysis/08_sensitivity_analyses/dnph1_eas_weights/run_dnph1_eas_mr.R`
5. `analysis/08_sensitivity_analyses/overlapping_ivs/`
6. any figure / summary scripts we can locate locally
7. AoU code once obtained

## People / Sources To Ask For Missing Code

- AoU PGS validation: Jacob Williams
- any final figure assembly scripts: whoever generated the submitted figures
