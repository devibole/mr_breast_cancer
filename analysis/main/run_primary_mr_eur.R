# Purpose: Ancestry-specific forward MR for European breast cancer GWAS using UKB-PPP EUR cis-pQTL instruments.

# Inputs:
#   - UKB-PPP EUR protein summary-statistic tar archives
#   - EUR breast cancer GWAS summary statistics
#   - rsID lookup RData containing `all_rsids`
#   - Olink annotation TSV
# Outputs:
#   - one MR input file per protein
#   - one MR result file per protein index

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript run_primary_mr_eur.R <protein_index>")
}
i <- as.numeric(args[1])

library(tidyverse)
library(data.table)
library(vroom)
library(MendelianRandomization)
library(mr.pivw)
library(mr.raps)
library(biomaRt)

create_temp_run_dir <- function(run_index, scratch_root = NULL) {
  slurm_id <- Sys.getenv("SLURM_JOB_ID", unset = "local")
  root_dir <- scratch_root
  if (is.null(root_dir) || identical(root_dir, "")) {
    tmp_env <- Sys.getenv("TMPDIR", unset = "")
    root_dir <- if (nzchar(tmp_env)) tmp_env else tempdir()
  }

  run_dir <- file.path(root_dir, paste0("mr_bc_", slurm_id, "_", run_index))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  run_dir
}

load_rsid_lookup <- function(rdata_path) {
  lookup_env <- new.env(parent = emptyenv())
  load(rdata_path, envir = lookup_env)
  if (!exists("all_rsids", envir = lookup_env, inherits = FALSE)) {
    stop("Expected object `all_rsids` was not found in: ", rdata_path)
  }
  data.table::as.data.table(get("all_rsids", envir = lookup_env))
}

standardize_alleles <- function(allele1, allele2) {
  idx <- allele1 > allele2
  list(
    allele1 = ifelse(idx, allele2, allele1),
    allele2 = ifelse(idx, allele1, allele2)
  )
}

convert_protein_name_to_oid <- function(protein_name) {
  sub(".*_(OID\\d+).*", "\\1", protein_name)
}

get_cis_regions_from_biomart <- function(protein_symbol, upstream = 1e6, downstream = 1e6) {
  ensembl <- biomaRt::useMart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    host = "https://grch37.ensembl.org"
  )

  genes <- biomaRt::getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"),
    filters = "hgnc_symbol",
    values = protein_symbol,
    mart = ensembl
  )

  if (nrow(genes) == 0) {
    return(genes)
  }

  genes$cis_start <- ifelse(genes$strand == 1, genes$start_position - upstream, genes$end_position - downstream)
  genes$cis_end <- ifelse(genes$strand == 1, genes$end_position + downstream, genes$start_position + upstream)
  genes$cis_start <- pmax(genes$cis_start, 1)
  genes
}

get_cis_region_table <- function(protein_name, olink_map_path, upstream = 1e6, downstream = 1e6) {
  olink_data <- data.table::fread(olink_map_path)
  protein_oid <- convert_protein_name_to_oid(protein_name)
  cis_regions <- olink_data[grepl(protein_oid, UKBPPP_ProteinID), ]

  if (nrow(cis_regions) > 0) {
    return(
      mutate(
        cis_regions,
        cis_start = ifelse(Strand == 1, gene_start - upstream, gene_end - downstream),
        cis_end = ifelse(Strand == 1, gene_end + downstream, gene_start + upstream),
        cis_start = pmax(cis_start, 1)
      )
    )
  }

  protein_symbol <- strsplit(protein_name, "_", fixed = TRUE)[[1]][1]
  biomart_regions <- get_cis_regions_from_biomart(protein_symbol, upstream = upstream, downstream = downstream)
  if (nrow(biomart_regions) == 0) {
    stop("No cis-region could be identified for protein: ", protein_name)
  }

  mutate(
    rename(biomart_regions, chr = chromosome_name),
    chr = as.numeric(chr),
    cis_start = ifelse(strand == 1, start_position - upstream, end_position - downstream),
    cis_end = ifelse(strand == 1, end_position + downstream, start_position + upstream),
    cis_start = pmax(cis_start, 1)
  )
}

read_exposure_from_tar <- function(tar_file, temp_dir, rsid_lookup) {
  utils::untar(tar_file, exdir = temp_dir)
  data_files <- list.files(temp_dir, pattern = "\\.gz$", full.names = TRUE, recursive = TRUE)
  if (length(data_files) == 0) {
    stop("No gzipped summary-statistic files were found in archive: ", tar_file)
  }

  data_files <- data_files[!grepl("chrX", data_files)]
  exposure_list <- lapply(data_files, data.table::fread)
  exposure_data <- data.table::rbindlist(exposure_list)
  exposure_data <- data.table::setDT(exposure_data)
  data.table::setkey(exposure_data, ID, N)
  exposure_data <- exposure_data[!duplicated(ID, fromLast = TRUE)]
  exposure_data[, c("allele1", "allele2") := standardize_alleles(ALLELE0, ALLELE1)]
  exposure_data <- merge(exposure_data, rsid_lookup[, c("ID", "rsid", "POS19")], by = "ID")
  exposure_data[, chr := as.numeric(CHROM)]
  exposure_data[, pos := as.numeric(POS19)]
  exposure_data
}

filter_to_cis_region <- function(merged_data, cis_regions) {
  pos_col <- if ("GENPOS" %in% colnames(merged_data)) "GENPOS" else "pos"
  merged_data <- data.table::setDT(merged_data)
  keep_idx <- apply(cis_regions, 1, function(cis_row) {
    merged_data[[pos_col]] >= as.numeric(cis_row["cis_start"]) &
      merged_data[[pos_col]] <= as.numeric(cis_row["cis_end"]) &
      merged_data$chr == as.numeric(cis_row["chr"])
  })
  merged_data[rowSums(keep_idx) > 0]
}

run_plink_clumping <- function(chromosome_file, temp_write_path, temp_dir, plink_binary) {
  chromosome_file_name <- tools::file_path_sans_ext(basename(chromosome_file))
  output_prefix <- file.path(temp_dir, paste0("clump_out_", chromosome_file_name))
  plink_command <- paste(
    plink_binary,
    "--bfile", tools::file_path_sans_ext(chromosome_file),
    "--clump", temp_write_path,
    "--clump-p1 0.5",
    "--clump-r2 0.01",
    "--clump-kb 2000",
    "--out", output_prefix
  )
  system(plink_command)
}

combine_clumped_files <- function(clumped_dir, output_file) {
  clumped_files <- list.files(clumped_dir, pattern = "^clump_out_.*\\.clumped$", full.names = TRUE)
  if (length(clumped_files) == 0) {
    stop("No .clumped files found in ", clumped_dir)
  }
  clumped_tables <- lapply(clumped_files, data.table::fread, header = TRUE)
  combined <- data.table::rbindlist(clumped_tables, fill = TRUE)
  data.table::fwrite(combined, output_file, sep = "\t")
  combined
}

safe_mr_result <- function(expr) {
  tryCatch(expr, error = function(e) NA)
}

build_mr_results_row <- function(iv_df, protein_label, snp_count, p_thresh) {
  mr_input_obj <- MendelianRandomization::mr_input(
    bx = iv_df$beta_exposure,
    bxse = iv_df$se_exposure,
    by = iv_df$beta_outcome,
    byse = iv_df$se_outcome,
    exposure = protein_label,
    outcome = "Breast Cancer",
    snps = iv_df$SNPID
  )

  ivw_fit <- if (identical(p_thresh, 5e-08)) safe_mr_result(MendelianRandomization::mr_ivw(mr_input_obj)) else NA
  median_fit <- if (identical(p_thresh, 5e-08)) safe_mr_result(MendelianRandomization::mr_median(mr_input_obj)) else NA
  egger_fit <- if (identical(p_thresh, 5e-08)) safe_mr_result(MendelianRandomization::mr_egger(mr_input_obj)) else NA
  pivw_fit <- safe_mr_result(
    mr.pivw::mr_pivw(
      Bx = iv_df$beta_exposure,
      Bxse = iv_df$se_exposure,
      By = iv_df$beta_outcome,
      Byse = iv_df$se_outcome,
      n.boot = 10000
    )
  )
  raps_fit <- safe_mr_result(
    mr.raps::mr.raps(
      data = data.frame(
        beta.exposure = iv_df$beta_exposure,
        beta.outcome = iv_df$beta_outcome,
        se.exposure = iv_df$se_exposure,
        se.outcome = iv_df$se_outcome
      )
    )
  )

  if (inherits(raps_fit, "mr.raps")) {
    raps_lower <- raps_fit$beta.hat - 1.96 * raps_fit$beta.se
    raps_upper <- raps_fit$beta.hat + 1.96 * raps_fit$beta.se
    raps_z <- raps_fit$beta.hat / raps_fit$beta.se
    raps_p <- 2 * stats::pnorm(abs(raps_z), lower.tail = FALSE)
  } else {
    raps_lower <- NA
    raps_upper <- NA
    raps_p <- NA
  }

  has_ivw <- isS4(ivw_fit)
  has_median <- isS4(median_fit)
  has_egger <- isS4(egger_fit)
  has_pivw <- isS4(pivw_fit)

  data.frame(
    Protein = protein_label,
    SNP_Count = snp_count,
    N = nrow(iv_df),
    ivs = paste(iv_df$SNPID, collapse = ", "),
    Estimate_ivw = if (has_ivw) ivw_fit@Estimate else NA,
    SE_ivw = if (has_ivw) ivw_fit@StdError else NA,
    CILower_ivw = if (has_ivw) ivw_fit@CILower else NA,
    CIUpper_ivw = if (has_ivw) ivw_fit@CIUpper else NA,
    P_Value_ivw = if (has_ivw) ivw_fit@Pvalue else NA,
    Estimate_median = if (has_median) median_fit@Estimate else NA,
    SE_median = if (has_median) median_fit@StdError else NA,
    CILower_median = if (has_median) median_fit@CILower else NA,
    CIUpper_median = if (has_median) median_fit@CIUpper else NA,
    P_Value_median = if (has_median) median_fit@Pvalue else NA,
    Estimate_egger = if (has_egger) egger_fit@Estimate else NA,
    SE_egger = if (has_egger) egger_fit@StdError.Est else NA,
    CILower_egger = if (has_egger) egger_fit@CILower.Est else NA,
    CIUpper_egger = if (has_egger) egger_fit@CIUpper.Est else NA,
    P_Value_egger = if (has_egger) egger_fit@Pvalue.Est else NA,
    Intercept_P_Value_egger = if (has_egger) egger_fit@Pvalue.Int else NA,
    Estimate_pivw = if (has_pivw) pivw_fit@Estimate else NA,
    SE_pivw = if (has_pivw) pivw_fit@StdError else NA,
    CILower_pivw = if (has_pivw) pivw_fit@CILower else NA,
    CIUpper_pivw = if (has_pivw) pivw_fit@CIUpper else NA,
    P_Value_pivw = if (has_pivw) pivw_fit@Pvalue else NA,
    Estimate_raps = if (inherits(raps_fit, "mr.raps")) raps_fit$beta.hat else NA,
    SE_raps = if (inherits(raps_fit, "mr.raps")) raps_fit$beta.se else NA,
    CILower_raps = raps_lower,
    CIUpper_raps = raps_upper,
    P_Value_raps = raps_p
  )
}

write_pipe_delim <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.table(x, file = file, sep = "|", row.names = FALSE, quote = FALSE)
}

config <- list(
  protein_archive_dir = "/data/BB_Bioinformatics/ProjectData/UKB_protein_sumstat/UKB-PPP_pGWAS_summary_statistics/European_discovery",
  olink_map_path = "/data/BB_Bioinformatics/ProjectData/UKB_protein_sumstat/Metadata/Protein_annotation/olink_protein_map_3k_v1.tsv",
  rsid_lookup_rdata = "/data/BB_Bioinformatics/DG/MR_bc/all_rsids.RData",
  outcome_path = "/data/BB_Bioinformatics/ProjectData/breast_cancer_sum_data/BCAC_EUR/bc_summary_gwas.txt",
  ld_reference_dir = "/data/BB_Bioinformatics/ProjectData/1000G_unrelated/1000G_full_data/GRCh37/EUR",
  plink_binary = "/usr/local/apps/plink/1.9.0-beta4.4/plink",
  input_export_dir = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/input_files/eur37",
  results_export_dir = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/results/eur37",
  scratch_root = "/lscratch",
  p_thresholds = c(5e-08)
)

temp_dir <- create_temp_run_dir(i, scratch_root = config$scratch_root)
on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

rsid_lookup <- load_rsid_lookup(config$rsid_lookup_rdata)
outcome_data <- vroom::vroom(config$outcome_path)
outcome_data <- data.table::setDT(outcome_data)
outcome_data[, c("allele1", "allele2") := standardize_alleles(Effect.Meta, Baseline.Meta)]
outcome_data[, c("chr", "pos") := list(chr.Onco, Position.Onco)]

tar_files <- list.files(config$protein_archive_dir, pattern = "\\.tar$", full.names = TRUE)
tar_file_to_process <- tar_files[i]
protein_label <- sub(".*/([^/]+)\\.tar$", "\\1", tar_file_to_process)

exposure_data <- read_exposure_from_tar(tar_file_to_process, temp_dir, rsid_lookup)

merged_data <- merge(
  exposure_data,
  outcome_data,
  by.x = c("chr", "pos", "allele1", "allele2"),
  by.y = c("chr", "pos", "allele1", "allele2")
)
merged_data <- data.table::setDT(merged_data)
merged_data <- merged_data[!duplicated(rsid)]

cis_regions <- get_cis_region_table(protein_label, config$olink_map_path)
merged_data <- filter_to_cis_region(merged_data, cis_regions)

merged_data <- merged_data %>%
  mutate(
    swap = ALLELE1 != Effect.Meta,
    ALLELE1_new = ifelse(swap, ALLELE0, ALLELE1),
    ALLELE0_new = ifelse(swap, ALLELE1, ALLELE0),
    BETA = ifelse(swap, -BETA, BETA),
    A1FREQ = ifelse(swap, 1 - A1FREQ, A1FREQ)
  ) %>%
  select(-ALLELE1, -ALLELE0) %>%
  rename(ALLELE1 = ALLELE1_new, ALLELE0 = ALLELE0_new)

merged_data$PVALUE <- 10^(-merged_data$LOG10P)
merged_data <- merged_data %>%
  filter(A1FREQ >= 0.01, A1FREQ <= 0.99, INFO > 0.3)

results_list <- list()

for (p_thresh in config$p_thresholds) {
  significant_snps <- merged_data[merged_data$PVALUE <= p_thresh, ]
  if (nrow(significant_snps) == 0) {
    warning("No SNPs passed preprocessing for threshold ", p_thresh, " in ", protein_label)
    next
  }

  clump_input <- data.frame(
    SNP = significant_snps$rsid,
    CHR = significant_snps$chr,
    BP = significant_snps$pos,
    P = significant_snps$PVALUE,
    A1 = significant_snps$ALLELE1,
    BETA = significant_snps$BETA
  )

  clump_input_path <- file.path(temp_dir, "clump_in.txt")
  write.table(clump_input, file = clump_input_path, row.names = FALSE, quote = FALSE)

  chromosome_files <- list.files(config$ld_reference_dir, pattern = "\\.bed$", full.names = TRUE)
  for (chromosome_file in chromosome_files) {
    run_plink_clumping(chromosome_file, clump_input_path, temp_dir, config$plink_binary)
  }

  clumps <- combine_clumped_files(temp_dir, file.path(temp_dir, "combined_clumps.txt"))

  iv_df <- significant_snps %>%
    filter(rsid %in% clumps$SNP) %>%
    transmute(
      SNPID = rsid,
      beta_exposure = BETA,
      se_exposure = SE,
      beta_outcome = Beta.meta,
      se_outcome = sdE.meta
    ) %>%
    as.data.frame()

  iv_df <- iv_df[iv_df$se_outcome > 0, , drop = FALSE]
  if (nrow(iv_df) == 0) {
    warning("No valid IVs remained after clumping for ", protein_label)
    next
  }

  mr_input_export <- significant_snps %>%
    filter(rsid %in% iv_df$SNPID) %>%
    transmute(
      SNPID = rsid,
      chr,
      pos,
      effect_allele = ALLELE1,
      other_allele = ALLELE0,
      beta_exposure = BETA,
      se_exposure = SE,
      freq_exposure = A1FREQ,
      pvalue_exposure = PVALUE,
      beta_outcome = Beta.meta,
      se_outcome = sdE.meta,
      pvalue_outcome = p.meta
    )

  write_pipe_delim(
    mr_input_export,
    file.path(config$input_export_dir, paste0(protein_label, ".txt"))
  )

  results_list[[as.character(p_thresh)]] <- build_mr_results_row(
    iv_df = iv_df,
    protein_label = protein_label,
    snp_count = nrow(significant_snps),
    p_thresh = p_thresh
  )
}

if (length(results_list) == 0) {
  stop("No MR results were generated for protein index ", i, " (", protein_label, ").")
}

combined_results <- bind_cols(results_list)
write_pipe_delim(
  combined_results,
  file.path(config$results_export_dir, paste0("result_", i, ".txt"))
)
