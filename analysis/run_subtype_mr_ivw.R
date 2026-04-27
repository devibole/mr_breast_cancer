
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Provide the tar-file index as the first command-line argument.")
}

job_index <- suppressWarnings(as.integer(args[1]))
if (is.na(job_index) || job_index < 1) {
  stop("The tar-file index must be a positive integer.")
}

slurm_job_id <- Sys.getenv("SLURM_JOB_ID", unset = "")
if (nzchar(slurm_job_id)) {
  temp_dir <- file.path("/lscratch", slurm_job_id, paste0("temp_", job_index))
} else {
  temp_dir <- file.path(tempdir(), paste0("temp_", job_index))
}

dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
message("Temporary directory: ", temp_dir)

library(tidyverse)
library(data.table)
library(vroom)
library(MendelianRandomization)
library(biomaRt)
library(meta)

# Paths
config <- list(
  tar_files_dir = "/data/BB_Bioinformatics/ProjectData/UKB_protein_sumstat/UKB-PPP_pGWAS_summary_statistics/European_discovery",
  olink_tsv_path = "/data/BB_Bioinformatics/ProjectData/UKB_protein_sumstat/Metadata/Protein_annotation/olink_protein_map_3k_v1.tsv",
  all_rsids_rdata = "/data/BB_Bioinformatics/DG/MR_bc/all_rsids.RData",
  outcome_data_path = "/data/BB_Bioinformatics/ProjectData/breast_cancer_sum_data/BCAC_EUR/bc_meta_subtype.txt",
  chromosome_dir = "/data/BB_Bioinformatics/ProjectData/1000G_unrelated/1000G_full_data/GRCh37/EUR",
  plink_binary = "/usr/local/apps/plink/1.9.0-beta4.4/plink",
  results_dir = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/subtypes/ivw/results37",
  input_files_dir = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/subtypes/ivw/input_files37"
)

# Load necessary data
load(config$all_rsids_rdata)
outcome_data_file <- vroom(config$outcome_data_path)
dir.create(config$results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$input_files_dir, recursive = TRUE, showWarnings = FALSE)

# Function to standardize alleles
standardize_alleles <- function(allele1, allele2) {
  idx <- allele1 > allele2
  allele1_std <- ifelse(idx, allele2, allele1)
  allele2_std <- ifelse(idx, allele1, allele2)
  list(
    allele1 = allele1_std,
    allele2 = allele2_std
  )
}

# Function to get cis-regions using biomaRt (fallback method)
get_cis_regions <- function(protein_names, upstream = 1e6, downstream = 1e6) {
  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", host="https://grch37.ensembl.org")
  
  genes <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"),
                 filters = "hgnc_symbol",
                 values = protein_names,
                 mart = ensembl)
  
  genes$cis_start <- ifelse(genes$strand == 1,
                            genes$start_position - upstream,
                            genes$end_position - downstream)
  genes$cis_end <- ifelse(genes$strand == 1,
                          genes$end_position + downstream,
                          genes$start_position + upstream)
  
  genes$cis_start <- pmax(genes$cis_start, 1)
  
  return(genes)
}

# Function to convert protein name to OID
convert_protein_name <- function(protein_name) {
  oid <- sub(".*_(OID\\d+).*", "\\1", protein_name)
  return(oid)
}

# Heterogeneity testing function
test_heterogeneity <- function(estimates, ses) {
  meta_analysis <- metagen(estimates, ses, comb.random = FALSE)
  return(list(
    Q = meta_analysis$Q,
    df = meta_analysis$df.Q,
    p_value = meta_analysis$pval.Q
  ))
}

# Define subtypes
subtypes <- list(
  Luminal_A = c("Luminal_A_log_or_meta", "Luminal_A_se_meta"),
  Luminal_B = c("Luminal_B_log_or_meta", "Luminal_B_se_meta"),
  HER2_Enriched = c("HER2_Enriched_log_or_meta", "HER2_Enriched_se_meta"),
  Luminal_B_HER2Neg = c("Luminal_B_HER2Neg_log_or_meta", "Luminal_B_HER2Neg_se_meta"),
  Triple_Neg = c("Triple_Neg_log_or_meta", "Triple_Neg_se_meta")
)

tar_files <- list.files(config$tar_files_dir, pattern = "\\.tar$", full.names = TRUE)
if (length(tar_files) == 0) {
  stop("No tar files were found in the exposure summary-statistics directory.")
}
if (job_index > length(tar_files)) {
  stop(sprintf("Requested tar-file index %s, but only %s tar files are available.", job_index, length(tar_files)))
}

tar_file_to_process <- tar_files[job_index]
print(tar_file_to_process)

# Extract and process exposure data
untar(tar_file_to_process, exdir = temp_dir)
data_files <- list.files(temp_dir, pattern = "\\.gz$", full.names = TRUE, recursive = TRUE)

if (length(data_files) == 0) {
  stop("No files found in the tar archive.")
}

data_files <- data_files[!grepl("chrX", data_files)]
list_of_data_frames <- lapply(data_files, fread)
exp_data <- rbindlist(list_of_data_frames)
exposure_data <- setDT(exp_data)

print(paste("Initial exposure SNPs:", nrow(exposure_data)))

# Remove duplicate IDs, keeping the one with highest N
setkey(exposure_data, ID, N)
exposure_data <- exposure_data[!duplicated(ID, fromLast = TRUE)]
print(paste("After removed repeated IDs:", nrow(exposure_data)))

# Standardize alleles and merge with rsids
exposure_data[, c("allele1", "allele2") := standardize_alleles(ALLELE0, ALLELE1)]
exposure_data <- merge(exposure_data, all_rsids[, c("ID", "rsid", "POS19")], by = "ID")
exposure_data$chr <- as.numeric(exposure_data$CHROM)
exposure_data$pos <- as.numeric(exposure_data$POS19)
print(paste("After rsid merge:", nrow(exposure_data)))

# Get cis-regions 
protein_name <- sub(".*/([^/]+)\\.tar$", "\\1", tar_file_to_process)
print(protein_name)

# Try Olink method first
olink_data <- fread(config$olink_tsv_path)
converted_protein_name <- convert_protein_name(protein_name)
print(converted_protein_name)
cis_regions <- olink_data[grepl(converted_protein_name, UKBPPP_ProteinID), ]

if (nrow(cis_regions) == 0) {
  # Fallback to biomaRt method
  print("Olink method failed, trying biomaRt...")
  protein_symbol <- strsplit(protein_name, "_")[[1]][1]
  cis_regions <- get_cis_regions(protein_symbol)
  
  if (nrow(cis_regions) == 0) {
    stop(paste("No cis-region found for protein:", protein_name))
  }
  
  # Convert biomaRt format to match expected format
  cis_regions <- cis_regions %>%
    rename(chr = chromosome_name) %>%
    mutate(
      chr = as.numeric(chr),
      cis_start = ifelse(strand == 1, start_position - 1e6, end_position - 1e6),
      cis_end = ifelse(strand == 1, end_position + 1e6, start_position + 1e6),
      cis_start = pmax(cis_start, 1)
    )
} else {
  # Use Olink format
  cis_regions <- cis_regions %>%
    mutate(
      cis_start = ifelse(Strand == 1, gene_start - 1e6, gene_end - 1e6),
      cis_end = ifelse(Strand == 1, gene_end + 1e6, gene_start + 1e6),
      cis_start = pmax(cis_start, 1)
    )
}

print("total PQTLs before cis filtering")
print(dim(exposure_data))

# Apply cis-region filtering
exposure_data <- setDT(exposure_data)

# Use appropriate position column based on data structure
pos_col <- ifelse("GENPOS" %in% colnames(exposure_data), "GENPOS", "pos")

exposure_data <- exposure_data[apply(cis_regions, 1, function(cis_row) {
  exposure_data[[pos_col]] >= as.numeric(cis_row["cis_start"]) & 
    exposure_data[[pos_col]] <= as.numeric(cis_row["cis_end"]) &
    exposure_data$chr == as.numeric(cis_row["chr"])
}) %>% rowSums() > 0]

print(paste("After cis-region filtering:", nrow(exposure_data)))

# Prepare outcome data (subtype specific format)
outcome_data <- setDT(outcome_data_file)
outcome_data[, c("chr", "pos", "A1", "A2") := tstrsplit(var_name, "_", fixed=TRUE)]
outcome_data[, c("allele1", "allele2") := standardize_alleles(Effect.Meta, Baseline.Meta)]
outcome_data$chr <- as.numeric(outcome_data$chr)
outcome_data$pos <- as.numeric(outcome_data$pos)

# Merge exposure and outcome data
merged_data <- merge(
  exposure_data,
  outcome_data,
  by.x = c("chr", "pos", "allele1", "allele2"),
  by.y = c("chr", "pos", "allele1", "allele2")
)
print(paste("After merging with outcome:", nrow(merged_data)))

setDT(merged_data)
merged_data <- merged_data[!duplicated(rsid)]

print(paste("After duplicate resolution:", nrow(merged_data)))
print("total SNPs")
print(dim(merged_data))

# Allele flipping - align exposure alleles with outcome alleles
merged_data <- merged_data %>%
  mutate(
    swap = (ALLELE1 != Effect.Meta),
    ALLELE1_new = ifelse(swap, ALLELE0, ALLELE1),
    ALLELE0_new = ifelse(swap, ALLELE1, ALLELE0),
    BETA = ifelse(swap, -BETA, BETA),
    A1FREQ = ifelse(swap, 1 - A1FREQ, A1FREQ)
  ) %>%
  select(-ALLELE1, -ALLELE0) %>%
  rename(ALLELE1 = ALLELE1_new, ALLELE0 = ALLELE0_new)

# Calculate p-values and apply filters
merged_data$PVALUE <- 10^(-merged_data$LOG10P)
merged_data <- merged_data %>%
  filter(A1FREQ >= 0.01 & A1FREQ <= 0.99 & INFO > 0.3)

print(paste("After MAF/INFO filtering:", nrow(merged_data)))

# Get significant SNPs
significant_snps <- merged_data[merged_data$PVALUE <= 5e-08, ]

tempWrite <- data.frame(
  SNP = significant_snps$rsid,
  CHR = significant_snps$chr, 
  BP = significant_snps$pos, 
  P = significant_snps$PVALUE, 
  A1 = significant_snps$ALLELE1, 
  BETA = significant_snps$BETA
)

if (nrow(tempWrite) == 0) {
  warning("No SNPs made it past preprocessing")
  # Create empty results and exit gracefully
  final_results <- data.frame(
    Protein = sub(".*/([^/]+)\\.tar$", "\\1", tar_file_to_process),
    SNP_Count = 0,
    N = 0,
    ivs = NA
  )
  
  for (subtype in names(subtypes)) {
    final_results[paste0("Estimate_", subtype)] <- NA
    final_results[paste0("SE_", subtype)] <- NA
    final_results[paste0("CILower_ivw_", subtype)] <- NA
    final_results[paste0("CIUpper_ivw_", subtype)] <- NA
    final_results[paste0("P_Value_ivw_", subtype)] <- NA
  }
  
  output_file_path <- file.path(config$results_dir, paste0("result_", job_index, ".csv"))
  write.csv(final_results, file = output_file_path, row.names = FALSE)
  
  # Cleanup and exit
  unlink(temp_dir, recursive = TRUE)
  quit(save = "no")
} else { 
  print("TempWrite with nrow")
  print(nrow(tempWrite))
  print(head(tempWrite))
}

temp_write_path <- paste0(temp_dir, "/clump_in.txt")
write.table(tempWrite, file = temp_write_path, row.names = FALSE, quote = FALSE)

run_plink_clumping <- function(chromosome_file, temp_write_path) {
  chromosome_file_name <- tools::file_path_sans_ext(basename(chromosome_file))
  output_prefix <- paste0(temp_dir, "/clump_out_", chromosome_file_name)
  plink_command <- paste(config$plink_binary,
                         "--bfile", tools::file_path_sans_ext(chromosome_file),
                         "--clump", temp_write_path,
                         "--clump-p1 0.5",
                         "--clump-r2 0.01",
                         "--clump-kb 2000",
                         "--out", output_prefix, sep=" ")
  
  print(paste("Running:", plink_command))
  system(plink_command)
}

chromosome_dir <- config$chromosome_dir
chromosome_files <- list.files(path = chromosome_dir, pattern = "\\.bed$", full.names = TRUE)

for (file in chromosome_files) {
  run_plink_clumping(file, temp_write_path)
}

message("Clumping finished")

combine_clumped_files <- function(clumped_dir, output_file) {
  command <- sprintf("
  # Ensure the output file is empty or does not exist
  > '%s'

  # Find all .clumped files and store them in an array
  CLUMPED_FILES=(%s/clump_out_chr*.clumped)

  # Check if there are any clumped files to process
  if [ ${#CLUMPED_FILES[@]} -eq 0 ]; then
    echo 'No .clumped files found in %s'
    exit 1
  fi

  # Process the first file: include the header and exclude empty lines
  awk 'NF > 0' \"${CLUMPED_FILES[0]}\" > '%s'

  # Process the rest of the files: skip the header and exclude empty lines
  for ((i=1; i<${#CLUMPED_FILES[@]}; i++)); do
    awk 'NR > 1 && NF > 0' \"${CLUMPED_FILES[i]}\" >> '%s'
  done

  echo 'Combined clumps are in %s'
  ", output_file, clumped_dir, clumped_dir, output_file, output_file, output_file)
  
  system(command, intern = FALSE)
}

clumped_dir <- temp_dir
output_file <- paste0(temp_dir, "/combined_clumps.txt")
combine_clumped_files(clumped_dir, output_file)

clumps <- fread(output_file, header = TRUE)
head(clumps)

# Process each subtype
subtype_results <- list()
final_results <- data.frame()

tempLabel <- sub(".*/([^/]+)\\.tar$", "\\1", tar_file_to_process)

for (subtype in names(subtypes)) {
  message("Processing subtype: ", subtype)
  log_or_column <- subtypes[[subtype]][1]
  se_column <- subtypes[[subtype]][2]
  
  iv_df <- significant_snps %>% 
    filter(significant_snps$rsid %in% clumps$SNP) %>%
    select(rsid, BETA, SE, log_or = all_of(log_or_column), se = all_of(se_column))
  
  colnames(iv_df) <- c("SNPID", "beta_exposure", "se_exposure", "beta_outcome", "se_outcome")
  iv_df <- as.data.frame(iv_df)
  
  tempSNPVec <- clumps$SNP
  tempIdx <- which(iv_df$SNPID %in% tempSNPVec)
  valid_idx <- tempIdx[iv_df[tempIdx, "se_outcome"] > 0]
  
  if (length(valid_idx) == 0) {
    warning(paste("No valid SNPs for MR analysis for subtype:", subtype))
    next
  }
  
  # Create extended input data for saving
  tempInput_df <- significant_snps %>%
    filter(rsid %in% iv_df$SNPID) %>%  # Only keep the SNPs that made it to iv_df
    select(
      SNPID = rsid,
      chr,
      pos, 
      effect_allele = ALLELE1,
      other_allele = ALLELE0,
      beta_exposure = BETA,
      se_exposure = SE,
      freq_exposure = A1FREQ,
      pvalue_exposure = PVALUE,
      beta_outcome = all_of(log_or_column),
      se_outcome = all_of(se_column)
    )
  
  # Save the extended input data
  protein_full_name <- paste0(tempLabel, "_", subtype)
  input_output_path <- file.path(config$input_files_dir, paste0(protein_full_name, ".txt"))
  write.table(tempInput_df, file = input_output_path, sep = "|", row.names = FALSE, quote = FALSE)
  
  tempInput <- mr_input(
    bx = iv_df[valid_idx, "beta_exposure"],
    bxse = iv_df[valid_idx, "se_exposure"],
    by = iv_df[valid_idx, "beta_outcome"],
    byse = iv_df[valid_idx, "se_outcome"],
    exposure = tempLabel,
    outcome = "Breast Cancer",
    snps = iv_df[valid_idx, "SNPID"]
  )
  
  # IVW analysis
  tempOutput_ivw <- tryCatch({
    mr_ivw(tempInput)
  }, error = function(e) {
    warning(paste("IVW analysis failed for subtype", subtype, ":", e$message))
    return(NA)
  })
  
  # Check if IVW analysis was successful
  if (!is.null(tempOutput_ivw) && !is.na(tempOutput_ivw)) {
    resultDF <- data.frame(
      Protein = tempLabel,
      SNP_Count = nrow(significant_snps),
      N = nrow(iv_df[valid_idx,]),
      ivs = paste(iv_df[valid_idx,]$SNPID, collapse = ", "),
      Subtype = subtype,
      Estimate = tempOutput_ivw@Estimate,
      SE = tempOutput_ivw@StdError,
      CILower_ivw = tempOutput_ivw@CILower,
      CIUpper_ivw = tempOutput_ivw@CIUpper,
      P_Value_ivw = tempOutput_ivw@Pvalue
    )
    
    subtype_results[[subtype]] <- resultDF
  } else {
    warning(paste("Failed to generate results for subtype:", subtype))
    
    # Create empty result for this subtype
    resultDF <- data.frame(
      Protein = tempLabel,
      SNP_Count = nrow(significant_snps),
      N = nrow(iv_df[valid_idx,]),
      ivs = paste(iv_df[valid_idx,]$SNPID, collapse = ", "),
      Subtype = subtype,
      Estimate = NA,
      SE = NA,
      CILower_ivw = NA,
      CIUpper_ivw = NA,
      P_Value_ivw = NA
    )
    
    subtype_results[[subtype]] <- resultDF
  }
}

# Perform heterogeneity test
subtype_estimates <- numeric()
subtype_ses <- numeric()
subtype_names <- character()

for (subtype in names(subtypes)) {
  if (!is.null(subtype_results[[subtype]])) {
    subtype_estimates <- c(subtype_estimates, subtype_results[[subtype]]$Estimate)
    subtype_ses <- c(subtype_ses, subtype_results[[subtype]]$SE)
    subtype_names <- c(subtype_names, subtype)
  }
}

# Initialize final_results with basic protein info
final_results <- data.frame(
  Protein = tempLabel,
  SNP_Count = nrow(significant_snps),
  N = length(unique(clumps$SNP)),
  ivs = paste(clumps$SNP, collapse = ", ")
)

# Add subtype results to final_results
for (subtype in names(subtypes)) {
  if (!is.null(subtype_results[[subtype]])) {
    subtype_data <- subtype_results[[subtype]]
    final_results[paste0("Estimate_", subtype)] <- subtype_data$Estimate
    final_results[paste0("SE_", subtype)] <- subtype_data$SE
    final_results[paste0("CILower_ivw_", subtype)] <- subtype_data$CILower_ivw
    final_results[paste0("CIUpper_ivw_", subtype)] <- subtype_data$CIUpper_ivw
    final_results[paste0("P_Value_ivw_", subtype)] <- subtype_data$P_Value_ivw
  } else {
    final_results[paste0("Estimate_", subtype)] <- NA
    final_results[paste0("SE_", subtype)] <- NA
    final_results[paste0("CILower_ivw_", subtype)] <- NA
    final_results[paste0("CIUpper_ivw_", subtype)] <- NA
    final_results[paste0("P_Value_ivw_", subtype)] <- NA
  }
}

# Add heterogeneity test results if applicable
if (length(subtype_estimates) > 1) {
  het_test <- tryCatch({
    test_heterogeneity(subtype_estimates, subtype_ses)
  }, error = function(e) {
    warning(paste("Heterogeneity test failed:", e$message))
    return(list(Q = NA, df = NA, p_value = NA))
  })
  
  final_results$Q_statistic <- het_test$Q
  final_results$Het_df <- het_test$df
  final_results$Het_P_Value <- het_test$p_value
  final_results$N_subtypes <- length(subtype_estimates)
  final_results$Subtypes_tested <- paste(subtype_names, collapse=";")
} else {
  final_results$Q_statistic <- NA
  final_results$Het_df <- NA
  final_results$Het_P_Value <- NA
  final_results$N_subtypes <- length(subtype_estimates)
  final_results$Subtypes_tested <- paste(subtype_names, collapse=";")
}

print("Final results:")
print(final_results)

# Save final combined results
output_file_path <- file.path(config$results_dir, paste0("result_", job_index, ".csv"))
write.csv(final_results, file = output_file_path, row.names = FALSE)

# Cleanup
unlink(temp_dir, recursive = TRUE)
rm(list = ls())
