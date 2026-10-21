# Release runner for the BCAC five-subtype MR heterogeneity approximation.
# Uses base R only. No example or study results are supplied or fabricated.
# CLI: Rscript run_analysis.R retained_inputs.csv original_results.csv new_output BCAC_only
#
# Inputs: two existing CSVs in the supplied template schemas; a NEW output
# directory; explicit tnbc_source = "BCAC_only"; code_dir containing this file
# and both core R scripts. Requires exactly 12 proteins, all five BCAC subtypes,
# and this project's confirmed common retained SNP set within each protein.
# Original MR estimates and default marginal SEs must reproduce before output.
# Return (invisibly): list(results, diagnostics, fits, covariance_comparison),
# with analysis roles added to results and the three-model comparison included.
# Writes all nine scenarios plus three clearly separated effective-count
# summaries and their comparison. These labels do not choose a manuscript
# primary result. Both MRE covariance constructions are sensitivity assumptions.

run_release_analysis <- function(input_csv, original_results_csv, output_dir,
                                 tnbc_source, code_dir) {
  code_files <- c("run_analysis.R", "subtype_heterogeneity.R", "shared_control_pipeline.R")
  if (any(!file.exists(file.path(code_dir, code_files)))) {
    stop("The runner and both core R scripts must be present in code_dir.")
  }
  if (!identical(tnbc_source, "BCAC_only")) {
    stop("Confirm BCAC alone by passing BCAC_only. BCAC plus CIMBA is unsupported.")
  }
  if (!file.exists(input_csv) || !file.exists(original_results_csv)) {
    stop("Both actual retained inputs and full-precision original MR results are required.")
  }
  if (file.exists(output_dir) || dir.exists(output_dir)) {
    stop("Choose a new output directory. Previous results will not be overwritten.")
  }
  d <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
  old <- read.csv(original_results_csv, stringsAsFactors = FALSE, check.names = FALSE)
  required_input <- c("protein", "snp", "subtype", "bx", "bxse", "by", "byse",
                      "effect_allele", "other_allele")
  required_old <- c("protein", "subtype", "estimate", "se_default")
  if (anyDuplicated(names(d)) || !all(required_input %in% names(d))) {
    stop("Retained-input columns are missing or duplicated; use the supplied template.")
  }
  if (anyDuplicated(names(old)) || !all(required_old %in% names(old))) {
    stop("Original-results columns are missing or duplicated; use the supplied template.")
  }
  if (nrow(d) == 0L || anyNA(d$protein) || any(trimws(d$protein) == "") ||
      length(unique(d$protein)) != 12L) {
    stop("The release run requires exactly 12 nonempty protein identifiers.")
  }
  if (nrow(old) != 60L || anyNA(old[required_old]) ||
      anyDuplicated(old[c("protein", "subtype")]) ||
      !setequal(unique(old$protein), unique(d$protein))) {
    stop("Original results must contain exactly the same 12 proteins and 60 unique subtype rows.")
  }
  for (key in c("estimate", "se_default")) {
    if (!is.numeric(old[[key]]) || any(!is.finite(old[[key]]))) {
      stop("Original MR ", key, " must contain full-precision finite numbers.")
    }
  }
  if (any(old$se_default <= 0)) stop("Original MR standard errors must be positive.")
  for (key in c("snp", "subtype")) {
    if (anyNA(d[[key]]) || any(trimws(d[[key]]) == "")) {
      stop("Retained-input identifiers must be nonempty: ", key)
    }
  }
  types <- bcac_subtype_counts()$subtype
  for (protein in unique(d$protein)) {
    if (!setequal(d$subtype[d$protein == protein], types) ||
        !setequal(old$subtype[old$protein == protein], types)) {
      stop("Both inputs require all five named BCAC subtypes for protein ", protein, ".")
    }
    retained <- d[d$protein == protein, ]
    reference_snps <- retained$snp[retained$subtype == types[1]]
    same_sets <- vapply(types, function(s) {
      setequal(retained$snp[retained$subtype == s], reference_snps)
    }, logical(1))
    if (!all(same_sets)) {
      stop("This project's release run requires the same retained SNP set across all five subtypes for protein ",
           protein, ". Resolve the source inputs; do not silently intersect the sets.")
    }
  }
  # Optional original n_snps values add a separate retained-instrument audit.
  if ("n_snps" %in% names(old)) {
    if (!is.numeric(old$n_snps) || any(!is.finite(old$n_snps)) ||
        any(old$n_snps < 1) || any(old$n_snps != floor(old$n_snps))) {
      stop("If supplied, original n_snps must be positive integers in all rows.")
    }
    observed <- aggregate(snp ~ protein + subtype, d, length)
    names(observed)[3] <- "n_retained"
    sizes <- merge(old[c("protein", "subtype", "n_snps")], observed,
                   by = c("protein", "subtype"), all = TRUE)
    if (anyNA(sizes) || any(sizes$n_snps != sizes$n_retained)) {
      stop("Original SNP counts do not match retained inputs. Resolve before analysis.")
    }
  }

  # The core validates harmonized data, covariance, and reproduction before
  # creating output. Its fixed tolerance is 1e-8 + 1e-6 * abs(original).
  result <- run_shared_control_pipeline(
    input_csv = input_csv, output_dir = output_dir, tnbc_source = tnbc_source,
    n_proteins = 12L, original_results_csv = original_results_csv)
  scenarios <- result$results
  role <- rep("count sensitivity or independence reference", nrow(scenarios))
  sampling <- scenarios$outcome_correlation == "effective_cases" &
              scenarios$variance_assumption == "sampling_only"
  independent_excess <- scenarios$outcome_correlation == "effective_cases" &
                        scenarios$variance_assumption == "default_independent_excess"
  scaled_correlation <- scenarios$outcome_correlation == "effective_cases" &
                        scenarios$variance_assumption == "default_scaled_correlation"
  independence <- scenarios$outcome_correlation == "independence_reference" &
                  scenarios$variance_assumption == "default_independent_excess"
  role[sampling] <- "shared-control sampling-covariance analysis"
  role[independent_excess] <- "MRE sensitivity: independent excess covariance"
  role[scaled_correlation] <- "MRE sensitivity: scaled total covariance"
  role[independence] <- "ordinary independence-based Q reference"
  scenarios$analysis_role <- role
  write.csv(scenarios, file.path(output_dir, "heterogeneity_scenarios.csv"), row.names = FALSE)
  write.csv(scenarios[sampling, ],
            file.path(output_dir, "shared_control_sampling_results.csv"), row.names = FALSE)
  write.csv(scenarios[independent_excess, ],
            file.path(output_dir, "sensitivity_independent_excess.csv"), row.names = FALSE)
  write.csv(scenarios[scaled_correlation, ],
            file.path(output_dir, "sensitivity_scaled_correlation.csv"), row.names = FALSE)
  write.csv(scenarios[independence, ],
            file.path(output_dir, "ordinary_Q_reference.csv"), row.names = FALSE)

  # Compare all three effective-count covariance choices without selecting one.
  # Suffixes identify the variance assumption for every numerical column.
  comparison_rows <- function(use, suffix) {
    ans <- scenarios[use, c("protein", "Q", "df", "p_value", "threshold", "below_threshold")]
    names(ans)[-1] <- paste0(names(ans)[-1], "_", suffix)
    ans
  }
  paired <- Reduce(function(x, y) merge(x, y, by = "protein"), list(
    comparison_rows(sampling, "sampling_only"),
    comparison_rows(independent_excess, "independent_excess"),
    comparison_rows(scaled_correlation, "scaled_correlation")))
  decisions <- paired[c("below_threshold_sampling_only", "below_threshold_independent_excess",
                        "below_threshold_scaled_correlation")]
  paired$threshold_decision_sensitive <- apply(decisions, 1L, function(z) length(unique(z)) > 1L)
  write.csv(paired, file.path(output_dir, "covariance_sensitivity_comparison.csv"), row.names = FALSE)

  paths <- c(input_csv, original_results_csv, file.path(code_dir, code_files))
  if (any(!file.exists(paths))) stop("Missing code file while recording provenance.")
  checksums <- data.frame(file = normalizePath(paths), md5 = unname(tools::md5sum(paths)))
  write.csv(checksums, file.path(output_dir, "input_and_code_checksums.csv"), row.names = FALSE)
  # Preserve exact submitted inputs alongside file hashes for reproducibility.
  if (!all(file.copy(c(input_csv, original_results_csv),
                     file.path(output_dir, c("retained_inputs_used.csv", "original_results_used.csv")),
                     overwrite = FALSE))) stop("Could not archive both actual input files.")
  release_notes <- c(
    "Release runner completed: all 12 proteins, all five BCAC subtypes, full-precision original MR comparison passed.",
    "Common retained SNP sets across all five subtypes checked separately for every protein.",
    "Shared-control sampling covariance: effective_cases/sampling_only.",
    "MRE sensitivity models: effective_cases/default_independent_excess and effective_cases/default_scaled_correlation.",
    "All nine scenarios are computed automatically; these output labels do not designate a manuscript primary result.",
    "All exposure estimates are conditioned on. bxse does not enter the variance formula.",
    "Independent-excess sensitivity alone assumes no cross-subtype covariance in the additional residual component.",
    "Scaled-correlation sensitivity separately specifies total covariance; it need not be an additive independent-error model.",
    "Effective case counts approximate shared-control correlations; original EM covariance is not recovered.",
    "Neither the MRE assumptions nor the count approximation are validated by this software run.",
    "P-value threshold: 0.05/12; four contrast degrees of freedom; plug-in chi-square approximation.",
    paste("Optional original instrument counts checked:", "n_snps" %in% names(old)),
    paste("UTC completion:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Command:", paste(commandArgs(), collapse = " ")))
  writeLines(release_notes, file.path(output_dir, "release_notes.txt"))
  message("Completed all 12 proteins with common subtype SNP sets. Review covariance_sensitivity_comparison.csv.")
  result$results <- scenarios
  result$covariance_comparison <- paired
  invisible(result)
}

# Permit source() for tests and interactive use; run automatically only via Rscript.
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 4L) {
    stop("Usage: Rscript run_analysis.R retained_inputs.csv original_results.csv NEW_output_directory BCAC_only", call. = FALSE)
  }
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  code_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
  source(file.path(code_dir, "subtype_heterogeneity.R"))
  source(file.path(code_dir, "shared_control_pipeline.R"))
  run_release_analysis(args[1], args[2], args[3], args[4], code_dir)
}
