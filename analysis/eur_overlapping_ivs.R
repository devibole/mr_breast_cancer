sid <- Sys.getenv('SLURM_JOB_ID')
print(sid)

library(dplyr)
library(data.table)
library(vroom)
library(MendelianRandomization)
library(utils)

# ── Paths ──────────────────────────────────────────────────────────────────────
tar_dir     <- "/data/BB_Bioinformatics/ProjectData/UKB_protein_sumstat/UKB-PPP_pGWAS_summary_statistics/European_discovery"
base_dir    <- "/data/BB_Bioinformatics/DG/MR_bc/2025_updated"
out_dir     <- "/data/BB_Bioinformatics/DG/MR_bc_sens/overlapping_ivs/eur"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load supporting data ───────────────────────────────────────────────────────
load("/data/BB_Bioinformatics/DG/MR_bc/all_rsids.RData")
outcome_data <- vroom("/data/BB_Bioinformatics/ProjectData/breast_cancer_sum_data/BCAC_EUR/bc_summary_gwas.txt")

# ── Load top 12 results files to get overlapping IVs ──────────────────────────
message("Loading top 12 results files...")
eur_res <- fread(file.path(base_dir, "top12_eur37.csv"))
eas_res <- fread(file.path(base_dir, "top12_eas37.csv"))
afr_res <- fread(file.path(base_dir, "top12_afr37.csv"))

proteins <- eur_res$Protein
message("Proteins: ", paste(proteins, collapse = ", "))

# ── Helper: parse IVs from ivs column ─────────────────────────────────────────
parse_ivs <- function(iv_string) {
  if (is.na(iv_string) || iv_string == "") return(character(0))
  trimws(unlist(strsplit(as.character(iv_string), ",")))
}

# ── Helper: standardize alleles ───────────────────────────────────────────────
standardize_alleles <- function(allele1, allele2) {
  idx <- allele1 > allele2
  list(
    allele1 = ifelse(idx, allele2, allele1),
    allele2 = ifelse(idx, allele1, allele2)
  )
}

# ── Prepare outcome data once ──────────────────────────────────────────────────
message("Preparing outcome data...")
setDT(outcome_data)
outcome_data[, c("allele1", "allele2") := standardize_alleles(Effect.Meta, Baseline.Meta)]
outcome_data[, c("chr", "pos") := list(chr.Onco, Position.Onco)]

# ── Master results lists ───────────────────────────────────────────────────────
all_summary <- list()
all_per_snp <- list()

p_thresholds <- c("5e-08", "5e-06")

# ══════════════════════════════════════════════════════════════════════════════
# Loop over proteins
# ══════════════════════════════════════════════════════════════════════════════
for (protein in proteins) {
  message("\n╔══════════════════════════════════════")
  message("║ Processing: ", protein)
  message("╚══════════════════════════════════════")
  
  # ── Find overlapping IVs across all three ancestries ────────────────────────
  overlapping_by_thresh <- list()
  for (p_thresh in p_thresholds) {
    iv_col <- paste0(p_thresh, ".ivs")
    eur_ivs <- parse_ivs(eur_res[eur_res$Protein == protein, ][[iv_col]])
    eas_ivs <- parse_ivs(eas_res[eas_res$Protein == protein, ][[iv_col]])
    afr_ivs <- parse_ivs(afr_res[afr_res$Protein == protein, ][[iv_col]])
    overlapping <- Reduce(intersect, list(eur_ivs, eas_ivs, afr_ivs))
    message("p=", p_thresh, " | EUR: ", length(eur_ivs),
            " EAS: ", length(eas_ivs), " AFR: ", length(afr_ivs),
            " | Overlapping: ", length(overlapping),
            " (", paste(overlapping, collapse = ", "), ")")
    overlapping_by_thresh[[p_thresh]] <- overlapping
  }
  
  # Get union of all overlapping IVs across both thresholds to extract from tar once
  all_ivs_needed <- unique(unlist(overlapping_by_thresh))
  
  if (length(all_ivs_needed) == 0) {
    warning("No overlapping IVs at any threshold for ", protein, " — skipping.")
    next
  }
  
  # ── Extract exposure data from tar file ─────────────────────────────────────
  tar_file <- file.path(tar_dir, paste0(protein, ".tar"))
  if (!file.exists(tar_file)) {
    warning("Tar file not found: ", tar_file)
    next
  }
  
  temp_dir <- paste0('/lscratch/', sid, '/temp_', protein)
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  message("Extracting tar file for ", protein, "...")
  untar(tar_file, exdir = temp_dir)
  
  data_files <- list.files(temp_dir, pattern = "\\.gz$", full.names = TRUE, recursive = TRUE)
  data_files <- data_files[!grepl("chrX", data_files)]
  
  if (length(data_files) == 0) {
    warning("No data files found in tar for ", protein)
    unlink(temp_dir, recursive = TRUE)
    next
  }
  
  list_of_dfs <- lapply(data_files, fread)
  exposure_data <- rbindlist(list_of_dfs)
  setDT(exposure_data)
  
  # Deduplicate
  setkey(exposure_data, ID, N)
  exposure_data <- exposure_data[!duplicated(ID, fromLast = TRUE)]
  
  # Standardize alleles and merge rsids
  exposure_data[, c("allele1", "allele2") := standardize_alleles(ALLELE0, ALLELE1)]
  exposure_data <- merge(exposure_data, all_rsids[, c("ID", "rsid", "POS19")], by = "ID")
  exposure_data$chr <- as.numeric(exposure_data$CHROM)
  exposure_data$pos <- as.numeric(exposure_data$POS19)
  
  message("Exposure SNPs after preprocessing: ", nrow(exposure_data))
  
  # ── Filter exposure to overlapping IVs BEFORE merge ─────────────────────────
  exposure_data <- exposure_data[rsid %in% all_ivs_needed]
  message("Exposure SNPs after IV filter: ", nrow(exposure_data))
  
  if (nrow(exposure_data) == 0) {
    warning("None of the overlapping IVs found in exposure data for ", protein)
    unlink(temp_dir, recursive = TRUE)
    next
  }
  
  # Filter outcome to just the positions we need before merging
  target_pos <- exposure_data$pos
  target_chr <- exposure_data$chr
  outcome_filtered <- outcome_data[pos %in% target_pos & chr %in% target_chr]
  
  # ── Merge with outcome ───────────────────────────────────────────────────────
  merged_data <- merge(
    exposure_data,
    outcome_filtered,
    by.x = c("chr", "pos", "allele1", "allele2"),
    by.y = c("chr", "pos", "allele1", "allele2")
  )
  rm(exposure_data, outcome_filtered)
  gc()
  merged_data <- setDT(merged_data)[!duplicated(rsid)]
  message("After merging with outcome: ", nrow(merged_data))
  
  # Allele flip
  merged_data <- merged_data %>%
    mutate(
      swap    = (ALLELE1 != Effect.Meta & ALLELE1 == Baseline.Meta),
      BETA    = ifelse(swap, -BETA, BETA),
      A1FREQ  = ifelse(swap, 1 - A1FREQ, A1FREQ)
    )
  
  # MAF and INFO filter
  merged_data <- merged_data %>%
    filter(A1FREQ >= 0.01 & A1FREQ <= 0.99 & INFO > 0.3)
  
  merged_data$PVALUE <- 10^(-merged_data$LOG10P)
  
  message("After MAF/INFO filter: ", nrow(merged_data))
  
  # ── Run MR per p-threshold ───────────────────────────────────────────────────
  for (p_thresh in p_thresholds) {
    overlapping_ivs <- overlapping_by_thresh[[p_thresh]]
    
    if (length(overlapping_ivs) == 0) {
      warning("No overlapping IVs for ", protein, " at ", p_thresh, " — skipping.")
      next
    }
    
    iv_data <- merged_data %>%
      filter(rsid %in% overlapping_ivs) %>%
      dplyr::select(
        SNPID          = rsid,
        chr, pos,
        effect_allele  = ALLELE1,
        other_allele   = ALLELE0,
        beta_exposure  = BETA,
        se_exposure    = SE,
        freq_exposure  = A1FREQ,
        pvalue_exposure = PVALUE,
        beta_outcome   = Beta.meta,
        se_outcome     = sdE.meta,
        pvalue_outcome = p.meta
      ) %>%
      filter(se_outcome > 0)
    
    message("p=", p_thresh, " | IVs found in merged data: ", nrow(iv_data))
    
    if (nrow(iv_data) == 0) {
      warning("No valid IVs in merged data for ", protein, " at ", p_thresh)
      next
    }
    
    # Save input file
    write.table(iv_data,
                file      = file.path(out_dir, paste0("eur_", protein, "_p", p_thresh, "_input.txt")),
                sep       = "|", row.names = FALSE, quote = FALSE)
    
    # MR
    mr_in <- mr_input(
      bx       = iv_data$beta_exposure,
      bxse     = iv_data$se_exposure,
      by       = iv_data$beta_outcome,
      byse     = iv_data$se_outcome,
      exposure = protein,
      outcome  = "Breast Cancer (EUR)",
      snps     = iv_data$SNPID
    )
    
    ivw_res <- tryCatch(mr_ivw(mr_in), error = function(e) NA)
    
    safe <- function(res, slot_name) {
      if (is.null(res) || inherits(res, "try-error") ||
          (is.logical(res) && is.na(res))) return(NA)
      tryCatch(slot(res, slot_name), error = function(e) NA)
    }
    
    # Per-SNP estimates
    all_per_snp[[paste(protein, "EUR", p_thresh, sep = "_")]] <- data.frame(
      Protein        = protein,
      Ancestry       = "EUR",
      p_threshold    = p_thresh,
      SNPID          = iv_data$SNPID,
      beta_exposure  = iv_data$beta_exposure,
      se_exposure    = iv_data$se_exposure,
      beta_outcome   = iv_data$beta_outcome,
      se_outcome     = iv_data$se_outcome,
      wald_ratio     = iv_data$beta_outcome / iv_data$beta_exposure,
      wald_ratio_se  = abs(iv_data$se_outcome / iv_data$beta_exposure)
    )
    
    # Summary
    all_summary[[paste(protein, "EUR", p_thresh, sep = "_")]] <- data.frame(
      Protein      = protein,
      Ancestry     = "EUR",
      p_threshold  = p_thresh,
      N_IVs        = nrow(iv_data),
      IVs          = paste(iv_data$SNPID, collapse = ", "),
      Estimate_ivw = safe(ivw_res, "Estimate"),
      SE_ivw       = safe(ivw_res, "StdError"),
      CILower_ivw  = safe(ivw_res, "CILower"),
      CIUpper_ivw  = safe(ivw_res, "CIUpper"),
      P_Value_ivw  = safe(ivw_res, "Pvalue")
    )
  }
  
  # Cleanup temp dir for this protein
  unlink(temp_dir, recursive = TRUE)
  message("Done: ", protein)
}

# ── Save results ───────────────────────────────────────────────────────────────
summary_out <- do.call(rbind, all_summary)
per_snp_out <- do.call(rbind, all_per_snp)

write.table(summary_out,
            file      = file.path(out_dir, "eur_overlapping_iv_mr_summary.txt"),
            sep       = "|", row.names = FALSE, quote = FALSE)

write.table(per_snp_out,
            file      = file.path(out_dir, "eur_overlapping_iv_mr_persnp.txt"),
            sep       = "|", row.names = FALSE, quote = FALSE)

message("\nAll done! Results saved to: ", out_dir)
rm(list = ls())