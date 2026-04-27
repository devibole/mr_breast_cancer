sid <- Sys.getenv("SLURM_JOB_ID")
if (sid == "") sid <- as.character(Sys.getpid())
print(sid)

suppressPackageStartupMessages({
  library(data.table)
  library(MendelianRandomization)
  library(utils)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
base_dir <- "/data/BB_Bioinformatics/DG/MR_bc/2025_updated"

# EAS top12 pQTL tar files
tar_dir <- "/vf/users/BB_Bioinformatics/ProjectData/UKB_protein_sumstat/UKB-PPP_pGWAS_summary_statistics/top12/East_Asian"

# EAS outcome input files (already ancestry-specific outcomes)
outcome_input_dir <- file.path(base_dir, "input_files/eas37")

# top12 overlap definition files
eur_top12 <- file.path(base_dir, "top12_eur37.csv")
eas_top12 <- file.path(base_dir, "top12_eas37.csv")
afr_top12 <- file.path(base_dir, "top12_afr37.csv")

# rsid map
all_rsids_file <- "/data/BB_Bioinformatics/DG/MR_bc/all_rsids.RData"

# output
out_dir <- "/data/BB_Bioinformatics/DG/MR_bc_sens/overlapping_ivs/eas"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# keep same style as your prior script
out_prefix <- "eas"

# ── Settings ───────────────────────────────────────────────────────────────────
p_thresholds <- c("5e-08", "5e-06")
p_threshold_num <- c("5e-08" = 5e-8, "5e-06" = 5e-6)

# ── Helpers ────────────────────────────────────────────────────────────────────
parse_ivs <- function(iv_string) {
  if (length(iv_string) == 0 || is.na(iv_string) || iv_string == "") return(character(0))
  unique(trimws(unlist(strsplit(as.character(iv_string), ","))))
}

safe_slot <- function(res, slot_name) {
  if (is.null(res) || inherits(res, "try-error") || (is.atomic(res) && length(res) == 1 && is.na(res))) return(NA_real_)
  tryCatch(slot(res, slot_name), error = function(e) NA_real_)
}

comp <- c(A = "T", T = "A", C = "G", G = "C")
allele_orient <- function(ea_exp, oa_exp, ea_out, oa_out) {
  ea_exp <- toupper(ea_exp); oa_exp <- toupper(oa_exp)
  ea_out <- toupper(ea_out); oa_out <- toupper(oa_out)
  
  if (ea_exp == ea_out && oa_exp == oa_out) return(1L)
  if (ea_exp == oa_out && oa_exp == ea_out) return(-1L)
  
  cea <- comp[ea_out]; coa <- comp[oa_out]
  if (!is.na(cea) && !is.na(coa)) {
    if (ea_exp == cea && oa_exp == coa) return(1L)
    if (ea_exp == coa && oa_exp == cea) return(-1L)
  }
  NA_integer_
}

# ── Load support data ──────────────────────────────────────────────────────────
message("Loading top12 files...")
eur_res <- fread(eur_top12)
eas_res <- fread(eas_top12)
afr_res <- fread(afr_top12)

proteins <- eas_res$Protein
message("Proteins: ", paste(proteins, collapse = ", "))

message("Loading rsid map...")
load(all_rsids_file)  # expects object all_rsids
if (!exists("all_rsids")) stop("all_rsids object not found in: ", all_rsids_file)

rs_map <- as.data.table(all_rsids)[, .(ID, SNPID = as.character(rsid), POS19 = as.integer(POS19))]
rs_map <- rs_map[!is.na(ID) & !is.na(SNPID)]
setkey(rs_map, ID)

# ── Master results lists ───────────────────────────────────────────────────────
all_summary <- list()
all_per_snp <- list()

# ═══════════════════════════════════════════════════════════════════════════════
# Main loop
# ═══════════════════════════════════════════════════════════════════════════════
for (protein in proteins) {
  message("\n╔══════════════════════════════════════")
  message("║ Processing: ", protein)
  message("╚══════════════════════════════════════")
  
  # ── Overlapping IVs across EUR/EAS/AFR by threshold ─────────────────────────
  overlapping_by_thresh <- list()
  for (p_thresh in p_thresholds) {
    iv_col <- paste0(p_thresh, ".ivs")
    
    eur_ivs <- parse_ivs(eur_res[Protein == protein][[iv_col]])
    eas_ivs <- parse_ivs(eas_res[Protein == protein][[iv_col]])
    afr_ivs <- parse_ivs(afr_res[Protein == protein][[iv_col]])
    
    overlapping <- Reduce(intersect, list(eur_ivs, eas_ivs, afr_ivs))
    overlapping_by_thresh[[p_thresh]] <- overlapping
    
    message("p=", p_thresh, " | EUR: ", length(eur_ivs),
            " EAS: ", length(eas_ivs), " AFR: ", length(afr_ivs),
            " | Overlapping: ", length(overlapping),
            " (", paste(overlapping, collapse = ", "), ")")
  }
  
  all_ivs_needed <- unique(unlist(overlapping_by_thresh, use.names = FALSE))
  if (length(all_ivs_needed) == 0) {
    warning("No overlapping IVs at any threshold for ", protein, " — skipping.")
    next
  }
  
  # map rsids -> IDs so we only read needed records from tar chunks
  needed_ids <- unique(rs_map[SNPID %in% all_ivs_needed, ID])
  if (length(needed_ids) == 0) {
    warning("No rs_map IDs found for overlapping IVs in ", protein)
    next
  }
  
  # ── Extract EAS exposure data from tar ───────────────────────────────────────
  tar_file <- file.path(tar_dir, paste0(protein, ".tar"))
  if (!file.exists(tar_file)) {
    warning("Tar file not found: ", tar_file)
    next
  }
  
  temp_dir <- file.path("/lscratch", sid, paste0("temp_", protein))
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
  
  cols_needed <- c("ID", "CHROM", "ALLELE1", "ALLELE0", "BETA", "SE", "A1FREQ", "INFO", "LOG10P", "N")
  chunks <- vector("list", length(data_files))
  k <- 0L
  
  for (f in data_files) {
    hdr <- names(fread(f, nrows = 0))
    use_cols <- intersect(cols_needed, hdr)
    if (length(use_cols) == 0) next
    
    dt <- fread(f, select = use_cols, showProgress = FALSE)
    if (!("ID" %in% names(dt))) next
    dt <- dt[ID %in% needed_ids]
    if (nrow(dt) == 0) next
    
    k <- k + 1L
    chunks[[k]] <- dt
  }
  
  if (k == 0L) {
    warning("No needed IV IDs found in tar for ", protein)
    unlink(temp_dir, recursive = TRUE)
    next
  }
  
  exposure_data <- rbindlist(chunks[seq_len(k)], fill = TRUE, use.names = TRUE)
  setDT(exposure_data)
  
  # Deduplicate by ID similar to your previous logic
  if ("N" %in% names(exposure_data)) {
    setorderv(exposure_data, c("ID", "N"), order = c(1L, -1L), na.last = TRUE)
  } else {
    setorder(exposure_data, ID)
  }
  exposure_data <- exposure_data[!duplicated(ID)]
  
  # Merge rsid map
  exposure_data <- merge(exposure_data, rs_map[, .(ID, SNPID, POS19)], by = "ID", all.x = TRUE, all.y = FALSE)
  exposure_data <- exposure_data[!is.na(SNPID)]
  
  # Build exposure fields
  exposure_data[, `:=`(
    chr = as.numeric(CHROM),
    pos = as.numeric(POS19),
    effect_allele_exp = toupper(ALLELE1),
    other_allele_exp  = toupper(ALLELE0),
    beta_exposure = as.numeric(BETA),
    se_exposure = as.numeric(SE),
    freq_exposure = as.numeric(A1FREQ),
    info = as.numeric(INFO),
    pvalue_exposure = 10^(-as.numeric(LOG10P))
  )]
  
  # Same QC style
  exposure_data <- exposure_data[(is.na(freq_exposure) | (freq_exposure >= 0.01 & freq_exposure <= 0.99)) &
                                   (is.na(info) | info > 0.3)]
  exposure_data <- exposure_data[!is.na(beta_exposure) & !is.na(se_exposure) & !is.na(pvalue_exposure)]
  exposure_data <- exposure_data[SNPID %in% all_ivs_needed]
  
  # dedup by SNPID
  setorder(exposure_data, SNPID, pvalue_exposure)
  exposure_data <- exposure_data[!duplicated(SNPID)]
  
  message("Exposure SNPs after preprocessing/filter: ", nrow(exposure_data))
  if (nrow(exposure_data) == 0) {
    warning("No overlapping IV exposure SNPs found for ", protein)
    unlink(temp_dir, recursive = TRUE)
    next
  }
  
  # ── Load EAS outcome info from existing input file ───────────────────────────
  outcome_file <- file.path(outcome_input_dir, paste0(protein, ".txt"))
  if (!file.exists(outcome_file)) {
    warning("Outcome input file not found: ", outcome_file)
    unlink(temp_dir, recursive = TRUE)
    next
  }
  
  outcome_data <- fread(outcome_file, sep = "|", showProgress = FALSE)
  required_out_cols <- c("SNPID", "effect_allele", "other_allele", "beta_outcome", "se_outcome", "pvalue_outcome")
  miss <- setdiff(required_out_cols, names(outcome_data))
  if (length(miss) > 0) {
    warning("Missing columns in outcome input for ", protein, ": ", paste(miss, collapse = ", "))
    unlink(temp_dir, recursive = TRUE)
    next
  }
  
  outcome_data <- outcome_data[, .(
    SNPID = as.character(SNPID),
    effect_allele_out = toupper(effect_allele),
    other_allele_out  = toupper(other_allele),
    beta_outcome = as.numeric(beta_outcome),
    se_outcome = as.numeric(se_outcome),
    pvalue_outcome = as.numeric(pvalue_outcome)
  )]
  outcome_data <- outcome_data[SNPID %in% all_ivs_needed]
  outcome_data <- outcome_data[!is.na(beta_outcome) & !is.na(se_outcome) & se_outcome > 0]
  setorder(outcome_data, SNPID)
  outcome_data <- outcome_data[!duplicated(SNPID)]
  
  # ── Merge by SNPID and harmonize alleles ─────────────────────────────────────
  merged_data <- merge(
    exposure_data[, .(SNPID, chr, pos,
                      effect_allele_exp, other_allele_exp,
                      beta_exposure, se_exposure, freq_exposure, pvalue_exposure)],
    outcome_data,
    by = "SNPID",
    all = FALSE
  )
  
  if (nrow(merged_data) == 0) {
    warning("No SNP overlap after merge for ", protein)
    unlink(temp_dir, recursive = TRUE)
    next
  }
  
  merged_data[, orient := mapply(
    allele_orient,
    effect_allele_exp, other_allele_exp,
    effect_allele_out, other_allele_out
  )]
  merged_data <- merged_data[!is.na(orient)]
  merged_data[orient == -1L, beta_outcome := -beta_outcome]
  
  message("After SNP merge + allele harmonization: ", nrow(merged_data))
  if (nrow(merged_data) == 0) {
    warning("No SNPs left after harmonization for ", protein)
    unlink(temp_dir, recursive = TRUE)
    next
  }
  
  # ── Run MR per threshold ─────────────────────────────────────────────────────
  for (p_thresh in p_thresholds) {
    overlapping_ivs <- overlapping_by_thresh[[p_thresh]]
    p_cut <- p_threshold_num[[p_thresh]]
    
    if (length(overlapping_ivs) == 0) {
      warning("No overlapping IVs for ", protein, " at ", p_thresh, " — skipping.")
      next
    }
    
    iv_data <- merged_data[
      SNPID %in% overlapping_ivs &
        pvalue_exposure <= p_cut &
        se_outcome > 0 &
        beta_exposure != 0
    ][, .(
      SNPID,
      chr, pos,
      effect_allele = effect_allele_exp,
      other_allele  = other_allele_exp,
      beta_exposure,
      se_exposure,
      freq_exposure,
      pvalue_exposure,
      beta_outcome,
      se_outcome,
      pvalue_outcome
    )]
    
    message("p=", p_thresh, " | IVs found in merged data: ", nrow(iv_data))
    
    if (nrow(iv_data) == 0) {
      warning("No valid IVs for ", protein, " at ", p_thresh)
      next
    }
    
    # Save input file (same style)
    write.table(
      iv_data,
      file = file.path(out_dir, paste0(out_prefix, "_", protein, "_p", p_thresh, "_input.txt")),
      sep = "|", row.names = FALSE, quote = FALSE
    )
    
    # MR IVW
    mr_in <- mr_input(
      bx = iv_data$beta_exposure,
      bxse = iv_data$se_exposure,
      by = iv_data$beta_outcome,
      byse = iv_data$se_outcome,
      exposure = protein,
      outcome = "Breast Cancer (EAS)",
      snps = iv_data$SNPID
    )
    
    ivw_res <- tryCatch(mr_ivw(mr_in), error = function(e) NA)
    
    # Per-SNP estimates (same columns as your old script)
    all_per_snp[[paste(protein, "EAS", p_thresh, sep = "_")]] <- data.frame(
      Protein       = protein,
      Ancestry      = "EAS",
      p_threshold   = p_thresh,
      SNPID         = iv_data$SNPID,
      beta_exposure = iv_data$beta_exposure,
      se_exposure   = iv_data$se_exposure,
      beta_outcome  = iv_data$beta_outcome,
      se_outcome    = iv_data$se_outcome,
      wald_ratio    = iv_data$beta_outcome / iv_data$beta_exposure,
      wald_ratio_se = abs(iv_data$se_outcome / iv_data$beta_exposure)
    )
    
    # Summary (same columns as your old script)
    all_summary[[paste(protein, "EAS", p_thresh, sep = "_")]] <- data.frame(
      Protein      = protein,
      Ancestry     = "EAS",
      p_threshold  = p_thresh,
      N_IVs        = nrow(iv_data),
      IVs          = paste(iv_data$SNPID, collapse = ", "),
      Estimate_ivw = safe_slot(ivw_res, "Estimate"),
      SE_ivw       = safe_slot(ivw_res, "StdError"),
      CILower_ivw  = safe_slot(ivw_res, "CILower"),
      CIUpper_ivw  = safe_slot(ivw_res, "CIUpper"),
      P_Value_ivw  = safe_slot(ivw_res, "Pvalue")
    )
  }
  
  unlink(temp_dir, recursive = TRUE)
  gc(verbose = FALSE)
  message("Done: ", protein)
}

# ── Save results (same format/style) ──────────────────────────────────────────
summary_out <- if (length(all_summary) > 0) do.call(rbind, all_summary) else data.frame()
per_snp_out <- if (length(all_per_snp) > 0) do.call(rbind, all_per_snp) else data.frame()

write.table(
  summary_out,
  file = file.path(out_dir, paste0(out_prefix, "_overlapping_iv_mr_summary.txt")),
  sep = "|", row.names = FALSE, quote = FALSE
)

write.table(
  per_snp_out,
  file = file.path(out_dir, paste0(out_prefix, "_overlapping_iv_mr_persnp.txt")),
  sep = "|", row.names = FALSE, quote = FALSE
)

message("\nAll done! Results saved to: ", out_dir)
rm(list = ls())
