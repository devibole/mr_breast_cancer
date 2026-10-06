#!/usr/bin/env Rscript

# Purpose: combine per-protein MR result files into one ancestry-level table.
# This is analysis-output collation required before cross-ancestry meta-analysis;
# it does not create or format any manuscript or supplementary table.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript compile_primary_results.R ",
    "<per_protein_results_dir> <output_csv> [expected_file_count]"
  )
}

results_dir <- args[[1]]
output_path <- args[[2]]
expected_file_count <- if (length(args) >= 3L) as.integer(args[[3]]) else NA_integer_

result_files <- list.files(
  results_dir,
  pattern = "^result_[0-9]+\\.txt$",
  full.names = TRUE
)

if (length(result_files) == 0L) {
  stop("No result_<index>.txt files were found in: ", results_dir)
}

result_index <- as.integer(sub("^result_([0-9]+)\\.txt$", "\\1", basename(result_files)))
result_files <- result_files[order(result_index)]
result_index <- sort(result_index)

if (anyDuplicated(result_index)) {
  stop("Duplicate result indexes were found in: ", results_dir)
}
if (!is.na(expected_file_count) && length(result_files) != expected_file_count) {
  stop(
    "Expected ", expected_file_count, " result files but found ",
    length(result_files), "."
  )
}
if (any(file.info(result_files)$size <= 0)) {
  stop("One or more per-protein result files are empty.")
}

read_result <- function(path) {
  result <- fread(path, sep = "|", check.names = FALSE)
  result[, result_index := as.integer(sub(
    "^result_([0-9]+)\\.txt$", "\\1", basename(path)
  ))]
  result
}

compiled <- rbindlist(lapply(result_files, read_result), fill = TRUE)
protein_columns <- grep("(^|\\.)Protein$", names(compiled), value = TRUE)
if (length(protein_columns) == 0L) {
  stop("No protein identifier column was found in the result files.")
}

protein_values <- lapply(protein_columns, function(column) as.character(compiled[[column]]))
compiled[, Protein := Reduce(dplyr::coalesce, protein_values)]

if (anyNA(compiled$Protein) || any(compiled$Protein == "")) {
  stop("At least one result row is missing a protein identifier.")
}
if (anyDuplicated(compiled$Protein)) {
  duplicates <- unique(compiled$Protein[duplicated(compiled$Protein)])
  stop("Duplicate proteins were found: ", paste(duplicates, collapse = ", "))
}

compiled[, (protein_columns) := NULL]
setcolorder(compiled, c("Protein", setdiff(names(compiled), "Protein")))
setorder(compiled, result_index)
compiled[, result_index := NULL]

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
fwrite(compiled, output_path)
message("Compiled ", nrow(compiled), " protein results to: ", output_path)
