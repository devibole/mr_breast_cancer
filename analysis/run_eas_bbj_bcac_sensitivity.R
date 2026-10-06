#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(MendelianRandomization)
})

base_dir <- "/data/BB_Bioinformatics/DG/MR_bc"
outcome_file <- paste0(
  "/vf/users/BB_Bioinformatics/ProjectData/breast_cancer_sum_data/",
  "BCAC_BBJ_EAS/EAS_BCAC_BBJ_meta_sumdata.txt"
)
instrument_dir <- file.path(base_dir, "2025_updated", "input_files", "eas37")
primary_results_file <- file.path(base_dir, "2025_updated", "eas37.csv")

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1]] else getwd()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

top12_symbols <- c(
  "SCAMP3", "LRRC25", "LRRC37A2", "RALB", "CASP8", "ANXA4",
  "ADM", "USP28", "RSPO3", "PARK7", "GCLM", "DNPH1"
)

primary_results <- fread(primary_results_file, check.names = FALSE)
top12 <- primary_results[sub("_.*$", "", Protein) %in% top12_symbols]
top12[, Protein_symbol := sub("_.*$", "", Protein)]
top12[, protein_order := match(Protein_symbol, top12_symbols)]
setorder(top12, protein_order)
top12[, protein_order := NULL]

if (uniqueN(top12$Protein_symbol) != 12L) {
  stop("Could not identify all 12 proteins in the primary EAS results.")
}

instrument_list <- lapply(top12$Protein, function(protein) {
  file <- file.path(instrument_dir, paste0(protein, ".txt"))
  x <- fread(file, sep = "|", check.names = FALSE)
  x[, `:=`(
    Protein = protein,
    Protein_symbol = sub("_.*$", "", protein)
  )]
  x
})

instruments <- rbindlist(instrument_list, fill = TRUE)
instruments[, pair_key := paste(
  pmin(effect_allele, other_allele),
  pmax(effect_allele, other_allele),
  sep = ":"
)]

if (anyDuplicated(instruments[, .(Protein, SNPID)])) {
  stop("Duplicate protein-SNP records in the saved EAS instruments.")
}

outcome_columns <- c(
  "unique_SNP_id",
  "effect_allele_BBJ", "non_effect_allele_BBJ", "Freq_effect_BBJ",
  "BETA_BBJ", "SE_BBJ", "P_BBJ", "N_eff_BBJ",
  "effect_allele_BCAC", "non_effect_allele_BCAC", "Freq_effect_BCAC",
  "BETA_BCAC", "SE_BCAC", "P_BCAC", "N_eff_BCAC",
  "effect_allele_meta", "non_effect_allele_meta",
  "BETA_meta", "SE_meta", "P_meta", "N_eff_meta"
)

outcome <- fread(
  outcome_file,
  select = outcome_columns,
  showProgress = TRUE
)

outcome[, c("chr", "pos", "id_a1", "id_a2") :=
          tstrsplit(unique_SNP_id, "_", fixed = TRUE)]
outcome <- outcome[chr != "X"]
outcome[, `:=`(
  chr = as.integer(chr),
  pos = as.integer(pos),
  pair_key = paste(pmin(id_a1, id_a2), pmax(id_a1, id_a2), sep = ":")
)]

instrument_keys <- unique(instruments[, .(chr, pos, pair_key)])
setkey(outcome, chr, pos, pair_key)
outcome <- outcome[instrument_keys, nomatch = 0]

if (nrow(outcome) != nrow(instrument_keys)) {
  stop("Not all EAS instrument variants were found in the outcome file.")
}
if (anyDuplicated(outcome[, .(chr, pos, pair_key)])) {
  stop("Duplicate outcome rows were found for the selected instruments.")
}

source_columns <- list(
  BBJ = list(
    effect = "effect_allele_BBJ",
    other = "non_effect_allele_BBJ",
    frequency = "Freq_effect_BBJ",
    beta = "BETA_BBJ",
    se = "SE_BBJ",
    p = "P_BBJ",
    n = "N_eff_BBJ"
  ),
  BCAC = list(
    effect = "effect_allele_BCAC",
    other = "non_effect_allele_BCAC",
    frequency = "Freq_effect_BCAC",
    beta = "BETA_BCAC",
    se = "SE_BCAC",
    p = "P_BCAC",
    n = "N_eff_BCAC"
  ),
  BBJ_BCAC_meta = list(
    effect = "effect_allele_meta",
    other = "non_effect_allele_meta",
    frequency = NULL,
    beta = "BETA_meta",
    se = "SE_meta",
    p = "P_meta",
    n = "N_eff_meta"
  )
)

get_slot <- function(object, slot_name) {
  if (is.null(object)) return(NA_real_)
  as.numeric(slot(object, slot_name))
}

run_source <- function(base, source_name, columns) {
  joined <- merge(
    base,
    outcome,
    by = c("chr", "pos", "pair_key"),
    all.x = TRUE
  )

  outcome_effect <- joined[[columns$effect]]
  outcome_other <- joined[[columns$other]]
  beta_outcome <- joined[[columns$beta]]
  se_outcome <- joined[[columns$se]]

  allele_valid <- !is.na(outcome_effect) & !is.na(outcome_other) &
    mapply(
      function(exposure_effect, outcome_effect, outcome_other) {
        exposure_effect == outcome_effect || exposure_effect == outcome_other
      },
      joined$effect_allele,
      outcome_effect,
      outcome_other
    )

  exposure_flipped <- allele_valid & joined$effect_allele != outcome_effect
  beta_exposure_aligned <- ifelse(
    exposure_flipped,
    -joined$beta_exposure,
    joined$beta_exposure
  )

  eligible <- allele_valid &
    is.finite(beta_exposure_aligned) &
    is.finite(joined$se_exposure) & joined$se_exposure > 0 &
    is.finite(beta_outcome) &
    is.finite(se_outcome) & se_outcome > 0

  fit <- NULL
  if (sum(eligible) > 0L) {
    mr_data <- mr_input(
      bx = beta_exposure_aligned[eligible],
      bxse = joined$se_exposure[eligible],
      by = beta_outcome[eligible],
      byse = se_outcome[eligible],
      exposure = unique(base$Protein),
      outcome = source_name,
      snps = joined$SNPID[eligible]
    )
    fit <- mr_ivw(mr_data)
  }

  result <- data.table(
    Protein = unique(base$Protein),
    Protein_symbol = unique(base$Protein_symbol),
    Outcome_source = source_name,
    N_primary_EAS_IVs = nrow(base),
    N_source_outcome_available = sum(eligible),
    N_source_outcome_unavailable = nrow(base) - sum(eligible),
    IVs_used = paste(joined$SNPID[eligible], collapse = ", "),
    Estimate_ivw = get_slot(fit, "Estimate"),
    SE_ivw = get_slot(fit, "StdError"),
    CILower_ivw = get_slot(fit, "CILower"),
    CIUpper_ivw = get_slot(fit, "CIUpper"),
    P_Value_ivw = get_slot(fit, "Pvalue")
  )

  frequency <- if (is.null(columns$frequency)) {
    rep(NA_real_, nrow(joined))
  } else {
    joined[[columns$frequency]]
  }

  audit <- data.table(
    Protein = joined$Protein,
    Protein_symbol = joined$Protein_symbol,
    Outcome_source = source_name,
    SNPID = joined$SNPID,
    chr = joined$chr,
    pos = joined$pos,
    effect_allele_exposure_original = joined$effect_allele,
    other_allele_exposure_original = joined$other_allele,
    beta_exposure_original = joined$beta_exposure,
    beta_exposure_aligned = beta_exposure_aligned,
    se_exposure = joined$se_exposure,
    effect_allele_outcome = outcome_effect,
    other_allele_outcome = outcome_other,
    effect_allele_frequency_outcome = frequency,
    beta_outcome = beta_outcome,
    se_outcome = se_outcome,
    p_outcome = joined[[columns$p]],
    n_effective_outcome = joined[[columns$n]],
    exposure_flipped = exposure_flipped,
    eligible = eligible,
    exclusion_reason = fifelse(
      eligible,
      "",
      fifelse(
        is.na(outcome_effect) | is.na(outcome_other) |
          is.na(beta_outcome) | is.na(se_outcome),
        "Source-specific outcome estimate unavailable",
        fifelse(
          !allele_valid,
          "Alleles could not be aligned",
          "Invalid beta or standard error"
        )
      )
    )
  )

  list(result = result, audit = audit)
}

result_list <- list()
audit_list <- list()

for (protein in top12$Protein) {
  base <- instruments[Protein == protein]

  for (source_name in c("BBJ_BCAC_meta", "BBJ", "BCAC")) {
    analysis <- run_source(base, source_name, source_columns[[source_name]])
    result_list[[length(result_list) + 1L]] <- analysis$result

    if (source_name %in% c("BBJ", "BCAC")) {
      audit_list[[length(audit_list) + 1L]] <- analysis$audit
    }
  }
}

results <- rbindlist(result_list)
results[, `:=`(
  OR = exp(Estimate_ivw),
  OR_lower = exp(CILower_ivw),
  OR_upper = exp(CIUpper_ivw),
  Direction = fifelse(
    Estimate_ivw > 0,
    "Positive",
    fifelse(Estimate_ivw < 0, "Negative", NA_character_)
  )
)]

audit <- rbindlist(audit_list)

coverage <- results[Outcome_source %in% c("BBJ", "BCAC"), .(
  Proteins = .N,
  Proteins_with_estimate = sum(is.finite(Estimate_ivw)),
  Total_primary_IV_records = sum(N_primary_EAS_IVs),
  Total_source_available_IV_records = sum(N_source_outcome_available),
  Total_source_unavailable_IV_records = sum(N_source_outcome_unavailable)
), by = Outcome_source]

if (nrow(results) != 36L) stop("Expected 36 result rows.")
if (nrow(audit) != 118L) stop("Expected 118 audit rows.")
if (sum(audit$eligible) != 94L) stop("Expected 94 eligible audit rows.")

fwrite(results, file.path(output_dir, "eas_top12_bbj_vs_bcac_ivw.csv"))
fwrite(audit, file.path(output_dir, "eas_top12_bbj_vs_bcac_instrument_audit.csv"))
fwrite(coverage, file.path(output_dir, "eas_top12_bbj_vs_bcac_coverage_summary.csv"))
