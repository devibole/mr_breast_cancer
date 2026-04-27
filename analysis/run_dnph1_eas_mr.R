sid <- Sys.getenv('SLURM_JOB_ID')
print(sid)
dir.create(paste0('/lscratch/', sid, '/temp_dnph1'), showWarnings = FALSE)
temp_dir <- paste0('/lscratch/', sid, '/temp_dnph1')
print(temp_dir)

library(tidyverse)
library(data.table)
library(MendelianRandomization)
library(vroom)
library(mr.pivw)
library(mr.raps)

# ── Paths ──────────────────────────────────────────────────────────────────────
pqtl_path    <- "/data/BB_Bioinformatics/DG/MR_bc_sens/eas_dnph1/hum0343.v3.qtl.v1/hum0343.v3.pqtl.v1/pqtl_releasedata_v0.1.tsv"
ld_dir       <- "/data/BB_Bioinformatics/ProjectData/1000G_unrelated/1000G_full_data/GRCh38/EAS"
out_dir      <- "/data/BB_Bioinformatics/DG/MR_bc_sens/eas_dnph1"
input_dir    <- file.path(out_dir, "input_files")
dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load outcome data ──────────────────────────────────────────────────────────
load("/data/BB_Bioinformatics/DG/MR_bc_sens/eas/outcome_data_file.RData")

# ── Helper functions ───────────────────────────────────────────────────────────
standardize_alleles <- function(allele1, allele2) {
  idx <- allele1 > allele2
  list(
    allele1 = ifelse(idx, allele2, allele1),
    allele2 = ifelse(idx, allele1, allele2)
  )
}

# ── Step 1: Define DNPH1 Ensembl ID and extract from pQTL file ────────────────
dnph1_ensembl_id <- "ENSG00000112667"
message("DNPH1 Ensembl ID: ", dnph1_ensembl_id)

# Read only DNPH1 rows (gene_id column has versioned IDs like ENSG00000136273.10)
message("Reading pQTL file and filtering for DNPH1...")
pqtl_raw <- vroom(pqtl_path, delim = "\t", show_col_types = FALSE)
exposure_data <- pqtl_raw %>%
  filter(grepl(dnph1_ensembl_id, gene_id, fixed = TRUE))

message("DNPH1 pQTL rows found: ", nrow(exposure_data))
if (nrow(exposure_data) == 0) stop("No DNPH1 rows found — check Ensembl ID or gene_id format.")

# ── Step 2: Parse hg19 position from variant_id_hg19 (chr:pos:ref:alt) ────────
# Format: 1:26963434:G:A
exposure_data <- exposure_data %>%
  mutate(
    chr   = as.numeric(sapply(strsplit(variant_id_hg19, ":"), `[`, 1)),
    pos   = as.numeric(sapply(strsplit(variant_id_hg19, ":"), `[`, 2)),
    REF   = sapply(strsplit(variant_id_hg19, ":"), `[`, 3),
    ALT   = sapply(strsplit(variant_id_hg19, ":"), `[`, 4)
  )

# Standardize alleles
std <- standardize_alleles(exposure_data$REF, exposure_data$ALT)
exposure_data$allele1 <- std$allele1
exposure_data$allele2 <- std$allele2

# Use rsid from file; flag missing rsids
exposure_data <- exposure_data %>%
  filter(!is.na(chr) & !is.na(pos))

message("After position parsing: ", nrow(exposure_data))

# ── Step 3: Prepare outcome data ───────────────────────────────────────────────
outcome_data <- setDT(outcome_data_file)
outcome_data[, c("chr", "pos", "A1", "A2") := tstrsplit(unique_SNP_id, "_", fixed = TRUE)]
outcome_data[, c("allele1", "allele2") := standardize_alleles(A1, A2)]
outcome_data <- outcome_data[chr != "X"]
outcome_data$chr <- as.numeric(outcome_data$chr)
outcome_data$pos <- as.numeric(outcome_data$pos)

# ── Step 4: Merge exposure and outcome ────────────────────────────────────────
merged_data <- merge(
  setDT(exposure_data),
  outcome_data,
  by = c("chr", "pos", "allele1", "allele2")
)
message("After merging with outcome: ", nrow(merged_data))
merged_data <- merged_data[!duplicated(rsid)]
message("After deduplication: ", nrow(merged_data))

# ── Step 5: Allele flip to align effect directions ────────────────────────────
# effect allele in pQTL file is ALT; align with A1 in outcome
merged_data <- merged_data %>%
  mutate(
    swap     = (ALT != A1 & ALT == A2),
    slope    = ifelse(swap, -slope, slope),
    maf      = ifelse(swap, 1 - maf,  maf)
  )

# ── Step 7: MAF filter ─────────────────────────────────────────────────────────
merged_data <- merged_data %>% filter(maf >= 0.01 & maf <= 0.99)
message("After MAF filter: ", nrow(merged_data))

# Rename for clarity downstream
merged_data$BETA      <- merged_data$slope
merged_data$SE        <- merged_data$slope_se
merged_data$PVALUE    <- merged_data$pval_nominal
merged_data$beta_out  <- merged_data$BETA_meta
merged_data$se_out    <- merged_data$SE_meta
merged_data$p_out     <- merged_data$P_meta

# ── Step 8: MR loop over p-value thresholds ───────────────────────────────────
p_thresholds <- c(5e-08, 5e-07, 5e-06)
results_list <- list()

for (p_thresh in p_thresholds) {
  message("\n── Processing p_thresh: ", p_thresh, " ──")
  sig_snps <- merged_data[merged_data$PVALUE <= p_thresh, ]
  
  tempWrite <- data.frame(
    SNP  = sig_snps$rsid,
    CHR  = sig_snps$chr,
    BP   = sig_snps$pos,
    P    = sig_snps$PVALUE,
    A1   = sig_snps$ALT,
    BETA = sig_snps$BETA
  )
  
  # Drop rows with missing rsid (needed for clumping)
  tempWrite <- tempWrite[!is.na(tempWrite$SNP) & tempWrite$SNP != "", ]
  
  if (nrow(tempWrite) == 0) {
    warning("No SNPs at p-threshold: ", p_thresh)
    next
  }
  message("SNPs going into clumping: ", nrow(tempWrite))
  
  temp_write_path <- paste0(temp_dir, "/clump_in.txt")
  write.table(tempWrite, file = temp_write_path, row.names = FALSE, quote = FALSE)
  
  # Clumping with EAS LD panel
  run_plink_clumping <- function(chromosome_file, temp_write_path) {
    out_prefix <- paste0(temp_dir, "/clump_out_", tools::file_path_sans_ext(basename(chromosome_file)))
    cmd <- paste(
      "module load plink/1.9.0-beta4.4;",
      "plink --bfile", tools::file_path_sans_ext(chromosome_file),
      "--clump", temp_write_path,
      "--clump-p1 0.5",
      "--clump-r2 0.01",
      "--clump-kb 2000",
      "--out", out_prefix
    )
    system(cmd)
  }
  
  chr_files <- list.files(ld_dir, pattern = "\\.bed$", full.names = TRUE)
  for (f in chr_files) run_plink_clumping(f, temp_write_path)
  message("Clumping done for p_thresh: ", p_thresh)
  
  # Combine clumped files
  output_clump <- paste0(temp_dir, "/combined_clumps.txt")
  combine_cmd <- sprintf('
    > "%s"
    CLUMPED_FILES=(%s/clump_out_chr*.clumped)
    if [ ${#CLUMPED_FILES[@]} -eq 0 ]; then echo "No clumped files found"; exit 1; fi
    awk "NF > 0" "${CLUMPED_FILES[0]}" > "%s"
    for ((i=1; i<${#CLUMPED_FILES[@]}; i++)); do
      awk "NR > 1 && NF > 0" "${CLUMPED_FILES[$i]}" >> "%s"
    done
  ', output_clump, temp_dir, output_clump, output_clump)
  system(combine_cmd)
  
  clumps <- fread(output_clump, header = TRUE)
  message("IVs after clumping: ", nrow(clumps))
  
  iv_df <- sig_snps %>%
    filter(rsid %in% clumps$SNP) %>%
    select(rsid, BETA, SE, beta_out, se_out) %>%
    rename(SNPID = rsid, beta_exposure = BETA, se_exposure = SE,
           beta_outcome = beta_out, se_outcome = se_out) %>%
    as.data.frame()
  
  valid_idx <- which(iv_df$se_outcome > 0)
  
  if (length(valid_idx) == 0) {
    warning("No valid SNPs for MR at p-threshold: ", p_thresh)
    next
  }
  
  # Save IV input file
  tempInput_df <- sig_snps %>%
    filter(rsid %in% iv_df$SNPID) %>%
    select(
      SNPID            = rsid,
      chr, pos,
      effect_allele    = ALT,
      other_allele     = REF,
      beta_exposure    = BETA,
      se_exposure      = SE,
      pvalue_exposure  = PVALUE,
      beta_outcome     = beta_out,
      se_outcome       = se_out,
      pvalue_outcome   = p_out
    )
  write.table(tempInput_df,
              file      = file.path(input_dir, paste0("eas_dnph1_p", p_thresh, ".txt")),
              sep       = "|", row.names = FALSE, quote = FALSE)
  
  # MR input object
  tempInput <- mr_input(
    bx       = iv_df[valid_idx, "beta_exposure"],
    bxse     = iv_df[valid_idx, "se_exposure"],
    by       = iv_df[valid_idx, "beta_outcome"],
    byse     = iv_df[valid_idx, "se_outcome"],
    exposure = "DNPH1",
    outcome  = "Breast Cancer (EAS)",
    snps     = iv_df[valid_idx, "SNPID"]
  )
  
  safe_extract <- function(res, slot_name) {
    if (is.null(res) || inherits(res, "try-error") ||
        (is.logical(res) && is.na(res))) return(NA)
    tryCatch(slot(res, slot_name), error = function(e) NA)
  }
  
  n_valid <- length(valid_idx)
  tempOutput_ivw    <- if (n_valid >= 1) tryCatch(mr_ivw(tempInput),    error = function(e) NA) else NA
  tempOutput_median <- if (n_valid >= 3) tryCatch(mr_median(tempInput),  error = function(e) NA) else NA
  tempOutput_egger  <- if (n_valid >= 3) tryCatch(mr_egger(tempInput),   error = function(e) NA) else NA
  
  tempOutput_pivw <- tryCatch(
    mr_pivw(Bx    = iv_df[valid_idx, "beta_exposure"],
            Bxse  = iv_df[valid_idx, "se_exposure"],
            By    = iv_df[valid_idx, "beta_outcome"],
            Byse  = iv_df[valid_idx, "se_outcome"],
            n.boot = 10000),
    error = function(e) NA
  )
  
  tempOutput_raps <- tryCatch(
    mr.raps(data = data.frame(
      beta.exposure = iv_df[valid_idx, "beta_exposure"],
      beta.outcome  = iv_df[valid_idx, "beta_outcome"],
      se.exposure   = iv_df[valid_idx, "se_exposure"],
      se.outcome    = iv_df[valid_idx, "se_outcome"])),
    error = function(e) NA
  )
  
  if (!is.null(tempOutput_raps) && !inherits(tempOutput_raps, "try-error") &&
      !(is.logical(tempOutput_raps) && is.na(tempOutput_raps))) {
    CI_bounds <- c(tempOutput_raps$beta.hat - 1.96 * tempOutput_raps$beta.se,
                   tempOutput_raps$beta.hat + 1.96 * tempOutput_raps$beta.se)
    raps_p    <- 2 * pnorm(abs(tempOutput_raps$beta.hat / tempOutput_raps$beta.se), lower.tail = FALSE)
  } else {
    CI_bounds <- c(NA, NA); raps_p <- NA
  }
  
  results_list[[as.character(p_thresh)]] <- data.frame(
    Protein              = "DNPH1",
    p_threshold          = p_thresh,
    SNP_Count            = nrow(sig_snps),
    N_IVs                = length(valid_idx),
    IVs                  = paste(iv_df[valid_idx, "SNPID"], collapse = ", "),
    Estimate_ivw         = safe_extract(tempOutput_ivw,    "Estimate"),
    SE_ivw               = safe_extract(tempOutput_ivw,    "StdError"),
    CILower_ivw          = safe_extract(tempOutput_ivw,    "CILower"),
    CIUpper_ivw          = safe_extract(tempOutput_ivw,    "CIUpper"),
    P_Value_ivw          = safe_extract(tempOutput_ivw,    "Pvalue"),
    Estimate_median      = safe_extract(tempOutput_median, "Estimate"),
    SE_median            = safe_extract(tempOutput_median, "StdError"),
    CILower_median       = safe_extract(tempOutput_median, "CILower"),
    CIUpper_median       = safe_extract(tempOutput_median, "CIUpper"),
    P_Value_median       = safe_extract(tempOutput_median, "Pvalue"),
    Estimate_egger       = safe_extract(tempOutput_egger,  "Estimate"),
    SE_egger             = safe_extract(tempOutput_egger,  "StdError.Est"),
    CILower_egger        = safe_extract(tempOutput_egger,  "CILower.Est"),
    CIUpper_egger        = safe_extract(tempOutput_egger,  "CIUpper.Est"),
    P_Value_egger        = safe_extract(tempOutput_egger,  "Pvalue.Est"),
    Intercept_P_egger    = safe_extract(tempOutput_egger,  "Pvalue.Int"),
    Estimate_pivw        = safe_extract(tempOutput_pivw,   "Estimate"),
    SE_pivw              = safe_extract(tempOutput_pivw,   "StdError"),
    CILower_pivw         = safe_extract(tempOutput_pivw,   "CILower"),
    CIUpper_pivw         = safe_extract(tempOutput_pivw,   "CIUpper"),
    P_Value_pivw         = safe_extract(tempOutput_pivw,   "Pvalue"),
    Estimate_raps        = if (!is.null(tempOutput_raps) && !inherits(tempOutput_raps, "try-error") &&
                               !(is.logical(tempOutput_raps) && is.na(tempOutput_raps))) tempOutput_raps$beta.hat else NA,
    SE_raps              = if (!is.null(tempOutput_raps) && !inherits(tempOutput_raps, "try-error") &&
                               !(is.logical(tempOutput_raps) && is.na(tempOutput_raps))) tempOutput_raps$beta.se  else NA,
    CILower_raps         = CI_bounds[1],
    CIUpper_raps         = CI_bounds[2],
    P_Value_raps         = raps_p
  )
}

# ── Save results ───────────────────────────────────────────────────────────────
combined_results <- do.call(rbind, results_list)
print(combined_results)
write.table(combined_results,
            file      = file.path(out_dir, "result_dnph1_eas.txt"),
            sep       = "|", row.names = FALSE, quote = FALSE)
message("Results saved to: ", file.path(out_dir, "result_dnph1_eas.txt"))

# ── Cleanup ────────────────────────────────────────────────────────────────────
unlink(temp_dir, recursive = TRUE)
rm(list = ls())
