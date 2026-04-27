suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(purrr)
  library(fs)
  library(writexl)
})

config <- list(
  base_dir = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/coloc/results37",
  combined_all_results_path = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/coloc/results37/combined_all_results.xlsx",
  combined_high_sig_path = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/coloc/results37/combined_highly_significant_snps.xlsx"
)

# Get list of all gene folders
gene_folders <- list.dirs(config$base_dir, recursive = FALSE, full.names = TRUE)

# Function to extract gene name from folder path
get_gene_name <- function(folder_path) {
  basename(folder_path)
}

# Function to find Excel file in each folder
find_excel_file <- function(folder_path) {
  gene_name <- get_gene_name(folder_path)
  expected_file <- file.path(folder_path, paste0(gene_name, "_coloc_results.xlsx"))
  
  if (file.exists(expected_file)) {
    return(expected_file)
  } else {
    # If exact name doesn't exist, look for any Excel file
    excel_files <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)
    if (length(excel_files) > 0) {
      return(excel_files[1])
    } else {
      return(NA)
    }
  }
}

# Get all Excel file paths
excel_files <- map_chr(gene_folders, find_excel_file)
excel_files <- excel_files[!is.na(excel_files)]

# Expected column names for validation
expected_cols_all_results <- c("protein", "protein_symbol", "analyzed_snp", "region", "n_snps", "PPH4", "PPH3")
expected_cols_highly_sig <- c("protein", "protein_symbol", "query_snp", "region", "H4", "rsid", 
                              "protein_beta", "protein_pval", "bc_beta", "bc_pval", "SNP_PP_H4", 
                              "chr", "pos", "joint_signal", "in_credible_set", "query_snp.1")

# Function to validate columns
validate_columns <- function(data, expected_cols, sheet_name, file_path) {
  actual_cols <- colnames(data)
  missing_cols <- setdiff(expected_cols, actual_cols)
  extra_cols <- setdiff(actual_cols, expected_cols)
  
  if (length(missing_cols) > 0) {
    warning(paste("Missing expected columns in", sheet_name, "from", basename(file_path), ":", paste(missing_cols, collapse = ", ")))
  }
  
  if (length(extra_cols) > 0) {
    cat("Extra columns found in", sheet_name, "from", basename(file_path), ":", paste(extra_cols, collapse = ", "), "\n")
  }
  
  return(data)
}

# Function to standardize data types for "Highly Significant SNPs"
standardize_highly_sig_types <- function(data) {
  # Convert specific columns to consistent types
  if ("joint_signal" %in% colnames(data)) {
    data$joint_signal <- as.character(data$joint_signal)
  }
  if ("in_credible_set" %in% colnames(data)) {
    data$in_credible_set <- as.character(data$in_credible_set)
  }
  if ("chr" %in% colnames(data)) {
    data$chr <- as.character(data$chr)
  }
  if ("pos" %in% colnames(data)) {
    data$pos <- as.numeric(data$pos)
  }
  if ("H4" %in% colnames(data)) {
    data$H4 <- as.numeric(data$H4)
  }
  if ("protein_beta" %in% colnames(data)) {
    data$protein_beta <- as.numeric(data$protein_beta)
  }
  if ("protein_pval" %in% colnames(data)) {
    data$protein_pval <- as.numeric(data$protein_pval)
  }
  if ("bc_beta" %in% colnames(data)) {
    data$bc_beta <- as.numeric(data$bc_beta)
  }
  if ("bc_pval" %in% colnames(data)) {
    data$bc_pval <- as.numeric(data$bc_pval)
  }
  if ("SNP_PP_H4" %in% colnames(data)) {
    data$SNP_PP_H4 <- as.numeric(data$SNP_PP_H4)
  }
  
  return(data)
}

# Function to standardize data types for "All Results"
standardize_all_results_types <- function(data) {
  if ("n_snps" %in% colnames(data)) {
    data$n_snps <- as.numeric(data$n_snps)
  }
  if ("PPH4" %in% colnames(data)) {
    data$PPH4 <- as.numeric(data$PPH4)
  }
  if ("PPH3" %in% colnames(data)) {
    data$PPH3 <- as.numeric(data$PPH3)
  }
  
  return(data)
}

# Function to read and add source information
read_sheet_with_source <- function(file_path, sheet_name) {
  gene_name <- get_gene_name(dirname(file_path))
  
  tryCatch({
    data <- read_excel(file_path, sheet = sheet_name)
    
    # Validate columns based on sheet type
    if (sheet_name == "All Results") {
      data <- validate_columns(data, expected_cols_all_results, sheet_name, file_path)
      data <- standardize_all_results_types(data)
    } else if (sheet_name == "Highly Significant SNPs") {
      data <- validate_columns(data, expected_cols_highly_sig, sheet_name, file_path)
      data <- standardize_highly_sig_types(data)
    }
    
    # Add source information
    data$Gene <- gene_name
    data$Source_File <- basename(file_path)
    return(data)
  }, error = function(e) {
    warning(paste("Could not read", sheet_name, "from", file_path, ":", e$message))
    return(NULL)
  })
}

# Read and combine "All Results" sheets
cat("Reading 'All Results' sheets...\n")
all_results_list <- map(excel_files, ~read_sheet_with_source(.x, "All Results"))
all_results_list <- compact(all_results_list)  # Remove NULL entries

if (length(all_results_list) > 0) {
  all_results_combined <- bind_rows(all_results_list)
  cat("Combined", length(all_results_list), "files for 'All Results'\n")
  cat("Total rows in combined 'All Results':", nrow(all_results_combined), "\n")
} else {
  warning("No 'All Results' sheets were successfully read")
  all_results_combined <- NULL
}

# Read and combine "Highly Significant SNPs" sheets
cat("Reading 'Highly Significant SNPs' sheets...\n")
highly_sig_list <- map(excel_files, ~read_sheet_with_source(.x, "Highly Significant SNPs"))
highly_sig_list <- compact(highly_sig_list)  # Remove NULL entries

if (length(highly_sig_list) > 0) {
  highly_sig_combined <- bind_rows(highly_sig_list)
  cat("Combined", length(highly_sig_list), "files for 'Highly Significant SNPs'\n")
  cat("Total rows in combined 'Highly Significant SNPs':", nrow(highly_sig_combined), "\n")
} else {
  warning("No 'Highly Significant SNPs' sheets were successfully read")
  highly_sig_combined <- NULL
}

# Display summary
cat("\n=== SUMMARY ===\n")
cat("Files processed:", length(excel_files), "\n")

if (!is.null(all_results_combined)) {
  cat("'All Results' - Rows:", nrow(all_results_combined), "Columns:", ncol(all_results_combined), "\n")
  cat("'All Results' - Genes included:", paste(unique(all_results_combined$Gene), collapse = ", "), "\n")
}

if (!is.null(highly_sig_combined)) {
  cat("'Highly Significant SNPs' - Rows:", nrow(highly_sig_combined), "Columns:", ncol(highly_sig_combined), "\n")
  cat("'Highly Significant SNPs' - Genes included:", paste(unique(highly_sig_combined$Gene), collapse = ", "), "\n")
}

# Optional: Save combined results to new Excel files
# Uncomment the lines below if you want to save the results

if (!is.null(all_results_combined)) {
  write_xlsx(all_results_combined, config$combined_all_results_path)
  cat("Saved 'All Results' to combined_all_results.xlsx\n")
}

if (!is.null(highly_sig_combined)) {
  write_xlsx(highly_sig_combined, config$combined_high_sig_path)
  cat("Saved 'Highly Significant SNPs' to combined_highly_significant_snps.xlsx\n")
}

# Display first few rows of each dataset
if (!is.null(all_results_combined)) {
  cat("\n=== First 5 rows of 'All Results' ===\n")
  print(head(all_results_combined, 5))
}

if (!is.null(highly_sig_combined)) {
  cat("\n=== First 5 rows of 'Highly Significant SNPs' ===\n")
  print(head(highly_sig_combined, 5))
}
