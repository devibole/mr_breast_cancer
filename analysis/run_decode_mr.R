# Standalone deCODE replication workflow.
# Hardcoded HPC paths are preserved for provenance and must be replaced for local reruns.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Provide the deCODE file index as the first command-line argument.")
}

job_index <- suppressWarnings(as.integer(args[1]))
if (is.na(job_index) || job_index < 1) {
  stop("The deCODE file index must be a positive integer.")
}

slurm_job_id <- Sys.getenv("SLURM_JOB_ID", unset = "")
if (nzchar(slurm_job_id)) {
  temp_dir <- file.path("/lscratch", slurm_job_id, paste0("temp_", job_index))
} else {
  temp_dir <- file.path(tempdir(), paste0("temp_", job_index))
}
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
message("Temporary directory: ", temp_dir)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(MendelianRandomization)
  library(vroom)
  library(biomaRt)
  library(mr.raps)
  library(rtracklayer)
})

if (!requireNamespace("mr.pivw", quietly = TRUE)) {
  stop("Install the 'mr.pivw' package before running this script.")
}
library(mr.pivw)

config <- list(
  excluded_variants_path = "/data/BB_Bioinformatics/ProjectData/DeCode_protein_sumstat/info/assocvariants.excluded.txt.gz",
  decode_files_dir = "/data/BB_Bioinformatics/ProjectData/DeCode_protein_sumstat/sumstats",
  all_rsids_rdata = "/data/BB_Bioinformatics/DG/MR_bc/all_rsids.RData",
  outcome_data_path = "/data/BB_Bioinformatics/ProjectData/breast_cancer_sum_data/BCAC_EUR/bc_summary_gwas.txt",
  info_file_path = "/data/BB_Bioinformatics/ProjectData/DeCode_protein_sumstat/info/assocvariants.annotated.txt.gz",
  chain_file_path = "/data/BB_Bioinformatics/Kevin/tools/liftover/hg38ToHg19.over.chain",
  chromosome_dir = "/data/BB_Bioinformatics/ProjectData/1000G_unrelated/1000G_full_data/GRCh37/EUR",
  input_files_dir = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/decode/input_files",
  results_dir = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/decode/results"
)

# Load excluded variants (Decode-specific)
excluded_variants <- fread(config$excluded_variants_path, fill = TRUE)
load(config$all_rsids_rdata)
outcome_data <- vroom(config$outcome_data_path)
dir.create(config$input_files_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$results_dir, recursive = TRUE, showWarnings = FALSE)

# Function to standardize alleles
standardize_alleles <- function(allele1, allele2) {
  idx <- allele1 > allele2
  list(allele1 = ifelse(idx, allele2, allele1), 
       allele2 = ifelse(idx, allele1, allele2))
}

# Load Decode protein data
data_files <- list.files(config$decode_files_dir, pattern = "\\.gz$", full.names = TRUE)
if (length(data_files) == 0) {
  stop("No deCODE summary-statistic files were found.")
}
if (job_index > length(data_files)) {
  stop(sprintf("Requested file index %s, but only %s deCODE files are available.", job_index, length(data_files)))
}

data_file_to_process <- data_files[job_index]
print(data_file_to_process)
exp <- vroom(data_file_to_process)

exposure_data <- setDT(exp)
print(paste("Initial exposure SNPs:", nrow(exposure_data)))

# Merge frequency data (Decode-specific)
info_data <- fread(config$info_file_path)

# Merge frequency data using the Name column
exposure_data <- merge(exposure_data, 
                      info_data[, .(Name, effectAlleleFreq)], 
                      by = "Name", 
                      all.x = TRUE)

setnames(exposure_data, "effectAlleleFreq", "decode_A1_freq")

print(paste("Added A1 frequency data. Match rate:", 
           round(sum(!is.na(exposure_data$decode_A1_freq))/nrow(exposure_data)*100, 1), "%"))

# Remove exact row duplicates
exposure_data <- unique(exposure_data)
print(paste("After removing exact duplicates:", nrow(exposure_data)))

# Remove chrX variants
exposure_data <- exposure_data[Chrom != "chrX"]
print(paste("After removing chrX:", nrow(exposure_data)))

# Exclude problematic variants (Decode-specific)
setkey(excluded_variants, Name) 
exposure_data <- exposure_data[!Name %in% excluded_variants$Name]
print(paste("After excluding variants:", nrow(exposure_data)))

# Standardize alleles
exposure_data[, c("allele1", "allele2") := standardize_alleles(effectAllele, otherAllele)]


# Create rsids for unmapped variants
exposure_data[is.na(rsids), rsids := paste0(
  gsub("chr", "", Chrom), ":",
  Pos, "_",
  otherAllele, "_",
  effectAllele
)]

print(paste("Total SNPs in raw exposure:", nrow(exposure_data)))

# Map to rsids using all_rsids
all_rsids[, first_char := substr(rsid, 1, 1)]
all_rsids[, is_pos_format := first_char %in% c("1","2","3","4","5","6","7","8","9")]

all_rsids[is_pos_format == TRUE, 
          rsids := paste0(
            substr(rsid, 1, regexpr(":", rsid) - 1),
            ":",
            POS38,                 
            "_",
            REF,
            "_",
            ALT
          )]
all_rsids[is_pos_format == FALSE, rsids := rsid]

exposure_data_mapped <- merge(exposure_data, all_rsids[,c("ID", "rsids", "POS19")], by="rsids", all.x=FALSE)
exposure_data_mapped <- exposure_data_mapped[!duplicated(exposure_data_mapped[, .(Chrom, Pos, effectAllele, otherAllele, Beta)])]

print(paste("SNPs successfully mapped to hg19:", nrow(exposure_data_mapped)))

# Liftover for unmapped SNPs
exposure_data_unmapped <- exposure_data[!rsids %in% all_rsids$rsids]
print(paste("SNPs not mapped - attempting liftover:", nrow(exposure_data_unmapped)))

if(nrow(exposure_data_unmapped) > 0) {
  # Apply liftover to unmapped SNPs
  chain <- import.chain(config$chain_file_path)
  
  # Prepare data for liftover function
  liftover_input <- exposure_data_unmapped[, .(CHR = as.numeric(gsub("chr", "", Chrom)), BP = Pos)]
  liftover_input <- cbind(liftover_input, exposure_data_unmapped)
  
  # Liftover function
  liftover38_19 <- function(dat) {
    gr_dat <- GRanges(seqnames = paste0("chr", dat$CHR), ranges = IRanges(start = dat$BP, width = 1))
    liftres <- liftOver(gr_dat, chain)
    liftres <- as.data.frame(liftres)
    dat$hg19_bp <- NA
    dat$hg19_chr <- NA
    if(nrow(liftres) > 0) {
      dat$hg19_bp[liftres$group] <- liftres$start
      dat$hg19_chr[liftres$group] <- as.character(liftres$seqnames)
      dat$hg19_chr <- gsub("chr", "", dat$hg19_chr)
    }
    return(dat)
  }
  
  # Apply liftover
  exposure_data_unmapped_lifted <- liftover38_19(liftover_input)
  
  # Keep only successfully lifted SNPs
  exposure_data_unmapped_lifted <- exposure_data_unmapped_lifted[!is.na(hg19_bp)]
  exposure_data_unmapped_lifted[, POS19 := hg19_bp]
  exposure_data_unmapped_lifted[, ID := NA]
  
  print(paste("SNPs successfully lifted from hg38 to hg19:", nrow(exposure_data_unmapped_lifted)))
  
  # Combine mapped and lifted data
  exposure_data <- rbindlist(list(exposure_data_mapped, exposure_data_unmapped_lifted), fill=TRUE)
} else {
  exposure_data <- exposure_data_mapped
}

print(paste("Total SNPs after rsid mapping + liftover:", nrow(exposure_data)))

exposure_data$chr <- as.numeric(gsub("chr", "", exposure_data$Chrom))
exposure_data$pos <- as.numeric(exposure_data$POS19)

print(paste("After rsids/pos19 conversion:", nrow(exposure_data)))

exposure_data <- setDT(exposure_data)
exposure_data[, rsid := rsids] 

print(paste("Rsid duplicates (multi-allelic variants - keeping all):", nrow(exposure_data) - length(unique(exposure_data$rsid))))

outcome_data_prepared <- outcome_data %>%
  rename(
    beta_outcome = Beta.meta,
    se_outcome = sdE.meta,
    pvalue_outcome = p.meta

  )

outcome_data_prepared <- setDT(outcome_data_prepared)
outcome_data_prepared <- unique(outcome_data_prepared)
outcome_data_prepared[, c("allele1", "allele2") := standardize_alleles(Effect.Meta, Baseline.Meta)]
outcome_data_prepared[, c("chr", "pos") := list(chr.Onco, Position.Onco)]

# Merge exposure and outcome data
merged_data <- merge(
  exposure_data,
  outcome_data_prepared,
  by.x = c("chr", "pos", "allele1", "allele2"),
  by.y = c("chr", "pos", "allele1", "allele2")
)
print(paste("After merging with breast cancer outcome:", nrow(merged_data)))

setDT(merged_data)

print(paste("After duplicate resolution:", nrow(merged_data)))
print("Total SNPs")
print(dim(merged_data))

# Get cis-regions function
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

# Extract protein information
filename <- sub(".*/([^/]+)\\.txt\\.gz$", "\\1", data_file_to_process)
protein_full_name <- filename
protein_symbol <- strsplit(filename, "_")[[1]][3]  # Extract protein symbol for biomart
protein_symbol

# Get cis-regions
cis_regions <- get_cis_regions(protein_symbol)
cis_regions <- cis_regions[cis_regions$chromosome_name %in% 1:22,]

if (nrow(cis_regions) == 0) {
  stop(paste("No cis-region found for protein symbol:", protein_symbol))
}

print("Total PQTLs before cis filtering")
print(dim(merged_data))

# Apply cis-region filtering
merged_data <- setDT(merged_data)

merged_data <- merged_data[apply(cis_regions, 1, function(cis_row) {
  merged_data$chr == as.numeric(cis_row["chromosome_name"]) &
    merged_data$pos >= cis_row["cis_start"] &
    merged_data$pos <= cis_row["cis_end"]
}) %>% rowSums() > 0]

print(paste("After cis-region filtering:", nrow(merged_data)))

if (nrow(merged_data) == 0) {
  stop("No SNPs in cis-region")
}

# Harmonize effect alleles for breast cancer outcome
merged_data <- merged_data %>%
  mutate(
    swap = (effectAllele != Effect.Meta),
    ALLELE1 = ifelse(swap, otherAllele, effectAllele),
    ALLELE0 = ifelse(swap, effectAllele, otherAllele),
    BETA = ifelse(swap, -Beta, Beta),
    rsid = rsid, 
    PVALUE = Pval
  )

# Apply MAF filtering
merged_data <- merged_data %>%
  filter(ImpMAF >= 0.01)

print(paste("After MAF filtering:", nrow(merged_data)))

# MR Analysis with multiple p-value thresholds (like UKB breast cancer script)
p_thresholds <- c(5e-08, 5e-07, 5e-06)
results_list <- list()

for (p_thresh in p_thresholds) {
  message("Processing p_thresh: ", p_thresh)
  significant_snps <- merged_data[merged_data$PVALUE <= p_thresh, ]
  
  tempWrite <- data.frame(
    SNP = significant_snps$rsid,
    CHR = significant_snps$chr, 
    BP = significant_snps$pos, 
    P = significant_snps$PVALUE, 
    A1 = significant_snps$ALLELE1, 
    BETA = significant_snps$BETA
  )
  
  if (nrow(tempWrite) == 0) {
    warning(paste("No SNPs made it past preprocessing for p-threshold:", p_thresh))
    next
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
    plink_command <- paste("module load plink/1.9.0-beta4.4;",
                           "plink --bfile", tools::file_path_sans_ext(chromosome_file),
                           "--clump", temp_write_path,
                           "--clump-p1 0.5",
                           "--clump-r2 0.1",
                           "--clump-kb 2000",
                           "--out", output_prefix, sep=" ")
    system(plink_command)
  }
  
  chromosome_files <- list.files(path = config$chromosome_dir, pattern = "\\.bed$", full.names = TRUE)
  if (length(chromosome_files) == 0) {
    stop("No PLINK reference files were found in the chromosome directory.")
  }
  
  for (file in chromosome_files) {
    run_plink_clumping(file, temp_write_path)
  }
  
  message("Clumping finished for p_thresh: ", p_thresh)
  
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
  print(head(clumps))
  
  iv_df <- significant_snps %>% filter(significant_snps$rsid %in% clumps$SNP) %>%
    select(rsid, BETA, SE, beta_outcome, se_outcome)
  
  colnames(iv_df) <- c("SNPID", "beta_exposure", "se_exposure", "beta_outcome", "se_outcome")
  print(head(iv_df))
  iv_df <- as.data.frame(iv_df)
  
  tempLabel <- protein_full_name
  tempSNPVec <- clumps$SNP
  tempIdx <- which(iv_df$SNPID %in% tempSNPVec)
  valid_idx <- tempIdx[iv_df[tempIdx, "se_outcome"] > 0]
  
  if (length(valid_idx) == 0) {
    warning("No valid SNPs for MR analysis at p-value threshold:", p_thresh)
    next
  }
  
  # Create extended input data for saving
  tempInput_df <- significant_snps %>%
    filter(rsid %in% iv_df$SNPID) %>%
    select(
      SNPID = rsid,
      chr,
      pos, 
      effect_allele = effectAllele,
      other_allele = otherAllele,
      decode_A1_freq,
      beta_exposure = BETA,
      se_exposure = SE,
      pvalue_exposure = PVALUE,
      beta_outcome,
      se_outcome,
      pvalue_outcome
    )

  # Save the extended input data
  output_file_path <- file.path(config$input_files_dir, paste0(protein_full_name, ".txt"))
  write.table(tempInput_df, file = output_file_path, sep = "|", row.names = FALSE, quote = FALSE)
  
  # Prepare MR input
  tempInput <- mr_input(
    bx = iv_df[valid_idx, "beta_exposure"],
    bxse = iv_df[valid_idx, "se_exposure"],
    by = iv_df[valid_idx, "beta_outcome"],
    byse = iv_df[valid_idx, "se_outcome"],
    exposure = protein_full_name,
    outcome = "Breast Cancer",
    snps = iv_df[valid_idx, "SNPID"]
  )
  
  # IVW analysis (only for genome-wide significant threshold)
  tempOutput_ivw <- if (p_thresh == 5e-08) tryCatch({
    mr_ivw(tempInput)
  }, error = function(e) NA)
  
  # Median analysis
  tempOutput_median <- if (p_thresh == 5e-08) tryCatch({
    mr_median(tempInput)
  }, error = function(e) NA)

  # Egger analysis
  tempOutput_egger <- if (p_thresh == 5e-08) tryCatch({
    mr_egger(tempInput)
  }, error = function(e) NA)
  
  # PIVW analysis
  tempOutput_pivw <- tryCatch({
    mr_pivw(Bx = iv_df[valid_idx, "beta_exposure"], 
            Bxse = iv_df[valid_idx, "se_exposure"], 
            By = iv_df[valid_idx, "beta_outcome"], 
            Byse = iv_df[valid_idx, "se_outcome"], 
            n.boot=10000)
  }, error = function(e) NA)
  
  # RAPS analysis
  tempOutput_raps <- tryCatch({
    mr.raps(data = data.frame(
      beta.exposure = iv_df[valid_idx, "beta_exposure"],
      beta.outcome = iv_df[valid_idx, "beta_outcome"],
      se.exposure = iv_df[valid_idx, "se_exposure"],
      se.outcome = iv_df[valid_idx, "se_outcome"]))
  }, error = function(e) NA)
  
  if (class(tempOutput_raps) == "mr.raps") {
    CI_bounds <- c(tempOutput_raps$beta.hat - 1.96 * tempOutput_raps$beta.se, 
                   tempOutput_raps$beta.hat + 1.96 * tempOutput_raps$beta.se)
    z_score <- tempOutput_raps$beta.hat/tempOutput_raps$beta.se
    p_value <- 2 * pnorm(abs(z_score), lower.tail = FALSE) 
    raps_estimate <- tempOutput_raps$beta.hat
    raps_se <- tempOutput_raps$beta.se
  } else {
    CI_bounds <- c(NA, NA)
    p_value <- NA
    raps_estimate <- NA
    raps_se <- NA
  }
  
  results_list[[as.character(p_thresh)]] <- data.frame(
    Protein = protein_full_name,
    Protein_Symbol = protein_symbol,
    SNP_Count = nrow(significant_snps),
    N = length(valid_idx),
    ivs = paste(iv_df[valid_idx,]$SNPID, collapse = ", "),
    Estimate_ivw = if (!is.null(tempOutput_ivw) && !is.na(tempOutput_ivw)) tempOutput_ivw@Estimate else NA,
    SE_ivw = if (!is.null(tempOutput_ivw) && !is.na(tempOutput_ivw)) tempOutput_ivw@StdError else NA,
    CILower_ivw = if (!is.null(tempOutput_ivw) && !is.na(tempOutput_ivw)) tempOutput_ivw@CILower else NA,
    CIUpper_ivw = if (!is.null(tempOutput_ivw) && !is.na(tempOutput_ivw)) tempOutput_ivw@CIUpper else NA,
    P_Value_ivw = if (!is.null(tempOutput_ivw) && !is.na(tempOutput_ivw)) tempOutput_ivw@Pvalue else NA,
    Estimate_median = if (!is.null(tempOutput_median) && !is.na(tempOutput_median)) tempOutput_median@Estimate else NA,
    SE_median = if (!is.null(tempOutput_median) && !is.na(tempOutput_median)) tempOutput_median@StdError else NA,
    CILower_median = if (!is.null(tempOutput_median) && !is.na(tempOutput_median)) tempOutput_median@CILower else NA,
    CIUpper_median = if (!is.null(tempOutput_median) && !is.na(tempOutput_median)) tempOutput_median@CIUpper else NA,
    P_Value_median = if (!is.null(tempOutput_median) && !is.na(tempOutput_median)) tempOutput_median@Pvalue else NA,
    Estimate_egger = if (!is.null(tempOutput_egger) && !is.na(tempOutput_egger)) tempOutput_egger@Estimate else NA,
    SE_egger = if (!is.null(tempOutput_egger) && !is.na(tempOutput_egger)) tempOutput_egger@StdError.Est else NA,
    CILower_egger = if (!is.null(tempOutput_egger) && !is.na(tempOutput_egger)) tempOutput_egger@CILower.Est else NA,
    CIUpper_egger = if (!is.null(tempOutput_egger) && !is.na(tempOutput_egger)) tempOutput_egger@CIUpper.Est else NA,
    P_Value_egger = if (!is.null(tempOutput_egger) && !is.na(tempOutput_egger)) tempOutput_egger@Pvalue.Est else NA,
    Intercept_P_Value_egger = if (!is.null(tempOutput_egger) && !is.na(tempOutput_egger)) tempOutput_egger@Pvalue.Int else NA,
    Estimate_pivw = if (!is.null(tempOutput_pivw) && !is.na(tempOutput_pivw)) tempOutput_pivw@Estimate else NA,
    SE_pivw = if (!is.null(tempOutput_pivw) && !is.na(tempOutput_pivw)) tempOutput_pivw@StdError else NA,
    CILower_pivw = if (!is.null(tempOutput_pivw) && !is.na(tempOutput_pivw)) tempOutput_pivw@CILower else NA,
    CIUpper_pivw = if (!is.null(tempOutput_pivw) && !is.na(tempOutput_pivw)) tempOutput_pivw@CIUpper else NA,
    P_Value_pivw = if (!is.null(tempOutput_pivw) && !is.na(tempOutput_pivw)) tempOutput_pivw@Pvalue else NA,
    Estimate_raps = raps_estimate,
    SE_raps = raps_se,
    CILower_raps = CI_bounds[1],
    CIUpper_raps = CI_bounds[2],
    P_Value_raps = p_value
  )
}

str(results_list)

if (length(results_list) == 0) {
  combined_results <- data.frame(
    Protein = protein_full_name,
    Protein_Symbol = protein_symbol
  )
} else {
  combined_results <- do.call(cbind, results_list)
}

# Save results to breast cancer decode directory
output_file_path <- file.path(config$results_dir, paste0("result_", job_index, ".txt"))
write.table(combined_results, file = output_file_path, sep="|", row.names = FALSE, quote = FALSE)

# Cleanup
unlink(temp_dir, recursive = TRUE)
dir.create(temp_dir)
rm(list = ls())
