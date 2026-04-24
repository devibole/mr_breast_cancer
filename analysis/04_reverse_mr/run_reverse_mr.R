# Standalone reverse MR workflow with breast cancer liability as the exposure.
# Hardcoded HPC paths are preserved for provenance and must be replaced for local reruns.

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

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(vroom)
  library(MendelianRandomization)
})

# ---- Paths (GRCh37 everywhere) ----
config <- list(
  tar_files_dir = "/data/BB_Bioinformatics/ProjectData/UKB_protein_sumstat/UKB-PPP_pGWAS_summary_statistics/European_discovery",
  chromosome_dir = "/data/BB_Bioinformatics/ProjectData/1000G_unrelated/1000G_full_data/GRCh37/EUR",
  out_dir = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/reverse_mr/results",
  input_dir = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/reverse_mr/input_files",
  all_rsids_rdata = "/data/BB_Bioinformatics/DG/MR_bc/all_rsids.RData",
  bc_gwas_path = "/data/BB_Bioinformatics/ProjectData/breast_cancer_sum_data/BCAC_EUR/bc_summary_gwas.txt",
  plink_binary = "/usr/local/apps/plink/1.9.0-beta4.4/plink"
)

dir.create(config$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$input_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Data loads ----
load(config$all_rsids_rdata)   # expects ID, rsid, POS19
exposure_data <- fread(config$bc_gwas_path) %>% as.data.table()

# ---- Helpers ----
standardize_alleles <- function(allele1, allele2) {
  idx <- allele1 > allele2
  allele1_std <- ifelse(idx, allele2, allele1)
  allele2_std <- ifelse(idx, allele1, allele2)
  list(allele1 = allele1_std, allele2 = allele2_std)
}

run_plink_clumping <- function(chromosome_file, temp_write_path, temp_dir) {
  output_prefix <- paste0(temp_dir, "/clump_out_", tools::file_path_sans_ext(basename(chromosome_file)))
  cmd <- paste(
    config$plink_binary,
    "--bfile", tools::file_path_sans_ext(chromosome_file),
    "--clump", temp_write_path,
    "--clump-p1 0.5",
    "--clump-r2 0.01",
    "--clump-kb 2000",
    "--out", output_prefix
  )
  system(cmd)
}

combine_clumped_files <- function(clumped_dir, output_file) {
  command <- sprintf("
    > '%s'
    CLUMPED_FILES=(%s/clump_out_chr*.clumped)
    if [ ${#CLUMPED_FILES[@]} -eq 0 ]; then
      echo 'No .clumped files found in %s'
      exit 1
    fi
    awk 'NF > 0' \"${CLUMPED_FILES[0]}\" > '%s'
    for ((i=1; i<${#CLUMPED_FILES[@]}; i++)); do
      awk 'NR > 1 && NF > 0' \"${CLUMPED_FILES[i]}\" >> '%s'
    done
  ", output_file, clumped_dir, clumped_dir, output_file, output_file)
  system(command, intern = FALSE)
}

# ---- Pick protein tar (outcome) ----
tar_files <- list.files(config$tar_files_dir, pattern = "\\.tar$", full.names = TRUE)
if (length(tar_files) == 0) {
  stop("No tar files were found in the exposure summary-statistics directory.")
}
if (job_index > length(tar_files)) {
  stop(sprintf("Requested tar-file index %s, but only %s tar files are available.", job_index, length(tar_files)))
}

tar_file_to_process <- tar_files[job_index]
protein_name <- sub(".*/([^/]+)\\.tar$", "\\1", tar_file_to_process)

untar(tar_file_to_process, exdir = temp_dir)
data_files <- list.files(temp_dir, pattern = "\\.gz$", full.names = TRUE, recursive = TRUE)
if (length(data_files) == 0) stop("No files found in the tar archive.")
data_files <- data_files[!grepl("chrX", data_files)]

outcome_raw <- rbindlist(lapply(data_files, fread))
setDT(outcome_raw)

# Deduplicate outcome by ID (keep largest N)
setkey(outcome_raw, ID, N)
outcome_raw <- outcome_raw[!duplicated(ID, fromLast = TRUE)]

# Standardize alleles & add rsids/pos37
outcome_raw[, c("allele1", "allele2") := standardize_alleles(ALLELE0, ALLELE1)]
outcome_data <- merge(outcome_raw, all_rsids[, .(ID, rsid, POS19)], by = "ID")
outcome_data[, `:=`(chr = as.numeric(CHROM), pos = as.numeric(POS19))]

# Exposure (BCAC) — standardize & add chr/pos to match
exposure_data[, c("allele1", "allele2") := standardize_alleles(Effect.Meta, Baseline.Meta)]
exposure_data <- exposure_data %>%
  transmute(chr = as.numeric(chr.Onco), pos = as.numeric(Position.Onco),
            EAFcontrols.Onco, Effect.Meta, Baseline.Meta, Beta.meta, sdE.meta, p.meta,
            allele1, allele2) %>% as.data.table()

# Merge on chr/pos/alleles
merged <- merge(exposure_data, outcome_data, by = c("chr", "pos", "allele1", "allele2"))

# Align exposure (BCAC) alleles to outcome ALLELE1
merged <- merged %>%
  mutate(
    swap = (Effect.Meta != ALLELE1),
    ALLELE1 = ifelse(swap, Baseline.Meta, Effect.Meta),
    ALLELE0 = ifelse(swap, Effect.Meta, Baseline.Meta),
    Beta.meta    = ifelse(swap, -Beta.meta, Beta.meta),
    A1FREQ  = ifelse(swap, 1 - EAFcontrols.Onco, EAFcontrols.Onco)
)

# Filters
merged <- merged %>% filter(A1FREQ >= 0.01, A1FREQ <= 0.99, INFO > 0.3)

# Significant exposure SNPs at 5e-8
significant <- merged %>% filter(p.meta <= 5e-8)

# Prepare clumping input
tempWrite <- significant %>% transmute(
  SNP = rsid, CHR = chr, BP = pos, P = p.meta, A1 = ALLELE1, BETA = Beta.meta
)

# Early out if none
out_file <- file.path(config$out_dir, paste0("result_", job_index, ".csv"))
if (nrow(tempWrite) == 0) {
  write.csv(data.frame(
    Protein = protein_name, 
    SNP_Count = 0,
    N = 0,
    ivs = NA,
    Estimate_ivw = NA, 
    SE_ivw = NA,
    CILower_ivw = NA, 
    CIUpper_ivw = NA, 
    P_Value_ivw = NA
  ), out_file, row.names = FALSE)
  unlink(temp_dir, recursive = TRUE)
  quit(save="no")
}

# Clumping
temp_in <- file.path(temp_dir, "clump_in.txt")
write.table(tempWrite, file=temp_in, row.names=FALSE, quote=FALSE)

chrom_files <- list.files(path = config$chromosome_dir, pattern = "\\.bed$", full.names = TRUE)
if (length(chrom_files) == 0) {
  stop("No PLINK reference files were found in the chromosome directory.")
}

for (cf in chrom_files) run_plink_clumping(cf, temp_in, temp_dir)

combined_path <- file.path(temp_dir, "combined_clumps.txt")
combine_clumped_files(temp_dir, combined_path)

clumps <- fread(combined_path, header=TRUE)

# Build IVs from clumped SNPs
iv <- significant %>%
  dplyr::filter(rsid %in% clumps$SNP) %>%
  dplyr::select(
    SNPID        = rsid,
    beta_exposure = Beta.meta,   # <- fix: was BETA.meta
    se_exposure   = sdE.meta,
    beta_outcome  = BETA,
    se_outcome    = SE
  )

valid_idx <- which(iv$se_outcome > 0)
if (length(valid_idx) == 0) {
  write.csv(data.frame(
    Protein = protein_name, 
    SNP_Count = nrow(significant),
    N = 0,
    ivs = NA,
    Estimate_ivw = NA, 
    SE_ivw = NA,
    CILower_ivw = NA, 
    CIUpper_ivw = NA, 
    P_Value_ivw = NA
  ), out_file, row.names = FALSE)
  unlink(temp_dir, recursive = TRUE)
  quit(save="no")
}

significant$PVALUE <- 10^(-significant$LOG10P)

tempInput_df <- significant %>%
  dplyr::filter(rsid %in% iv$SNPID) %>%  
  dplyr::select(
    SNPID = rsid,
    chr,
    pos, 
    effect_allele = ALLELE1,
    other_allele = ALLELE0,
    beta_exposure = Beta.meta,
    se_exposure = sdE.meta,
    freq_gwas = A1FREQ,
    pvalue_exposure = p.meta,
    beta_outcome = BETA,
    se_outcome = SE,
    pvalue_outcome = PVALUE
  )

# Save input file for future reference
input_file_path <- file.path(config$input_dir, paste0(protein_name, ".txt"))
write.table(tempInput_df, file = input_file_path, sep = "|", row.names = FALSE, quote = FALSE)

# Format input for MR packages
mr_data <- mr_input(
  bx = iv$beta_exposure[valid_idx],
  bxse = iv$se_exposure[valid_idx],
  by = iv$beta_outcome[valid_idx],
  byse = iv$se_outcome[valid_idx],
  exposure = "Breast Cancer",
  outcome = protein_name,
  snps = iv$SNPID[valid_idx]
)

# MR IVW analysis
ivw_result <- tryCatch(
  mr_ivw(mr_data),
  error = function(e) NULL
)

if (!is.null(ivw_result)) {
  res <- data.frame(
    Protein = protein_name,
    SNP_Count = nrow(significant),
    N = length(valid_idx),
    ivs = paste(iv$SNPID[valid_idx], collapse = ", "),
    Estimate_ivw = ivw_result@Estimate,
    SE_ivw = ivw_result@StdError,
    CILower_ivw = ivw_result@CILower,
    CIUpper_ivw = ivw_result@CIUpper,
    P_Value_ivw = ivw_result@Pvalue
  )
} else {
  res <- data.frame(
    Protein = protein_name, 
    SNP_Count = nrow(significant),
    N = length(valid_idx),
    ivs = paste(iv$SNPID[valid_idx], collapse = ", "),
    Estimate_ivw = NA, 
    SE_ivw = NA,
    CILower_ivw = NA, 
    CIUpper_ivw = NA, 
    P_Value_ivw = NA
  )
}

write.csv(res, out_file, row.names = FALSE)
unlink(temp_dir, recursive = TRUE)
