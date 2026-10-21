# BCAC-only five-subtype shared-control approximation, 13 October 2026.
# Source subtype_heterogeneity.R first. Base R only; no automatic analysis.
# Read ANALYST_GUIDE.html or ANALYST_GUIDE.md before interpreting the results.

# No arguments. Returns five labeled rows of effective/observed CASE counts
# and the common control count. No data are downloaded or inferred.
bcac_subtype_counts <- function() {
  # User-supplied workbook, Sheet1!A5:D9. These are effective CASE counts,
  # not balanced case-control effective total sample sizes.
  data.frame(
    subtype = c("Luminal_A", "Luminal_B_HER2Neg", "Luminal_B",
                "HER2_Enriched", "Triple_Neg"),
    effective_cases = c(45253, 6350, 6427, 2884, 8602),
    observed_cases = c(27481, 6799, 6502, 2793, 7174),
    shared_controls = rep(91477, 5), stringsAsFactors = FALSE)
}

# Simple null approximation with disjoint cases and identical controls.
# Substituting effective cases is a WORKING approximation for the EM GWAS.
# Inputs: positive named case-count vector and one positive shared-control
# count. Returns a named K x K outcome-error correlation matrix.
shared_control_correlation <- function(cases, controls) {
  if (!is.numeric(cases) || length(cases) < 2L ||
      any(!is.finite(cases)) || any(cases <= 0) ||
      is.null(names(cases)) || anyDuplicated(names(cases)) ||
      anyNA(names(cases)) || any(names(cases) == "")) stop("Invalid named case counts.")
  if (!is.numeric(controls) || length(controls) != 1L ||
      !is.finite(controls) || controls <= 0) stop("One positive control count required.")
  loading <- sqrt(cases / (cases + controls))
  R <- outer(loading, loading)
  diag(R) <- 1
  dimnames(R) <- list(names(cases), names(cases))
  R
}

# Reproduce the marginal default mr_ivw simple-weight SEs, independent SNPs.
# bxse is intentionally not used in this first-order calculation.
# Inputs use the same aligned bx, by, byse as ivw_covariance_first_order();
# fit is that function's result. Returns one row per subtype with SNP count,
# model, estimate, first-order SE, residual Q/df/SE, inflation, and default SE.
ivw_default_dispersion <- function(bx, by, byse, fit) {
  available <- !is.na(by)
  residual <- by - outer(bx, fit$theta)
  terms <- (residual / byse)^2
  terms[!available] <- 0
  q <- colSums(terms)
  n <- colSums(available)
  rse <- setNames(rep(NA_real_, length(n)), names(n))
  rse[n > 1L] <- sqrt(q[n > 1L] / (n[n > 1L] - 1L))
  inflation <- setNames(rep(1, length(n)), names(n))
  inflation[n >= 4L] <- pmax(1, rse[n >= 4L])
  if (any(!is.finite(q)) || any(!is.finite(inflation))) stop("Nonfinite IVW dispersion.")
  data.frame(subtype = names(n), n_snps = as.integer(n),
             model = ifelse(n == 1L, "Wald", ifelse(n <= 3L, "fixed", "multiplicative_random")),
             estimate = unname(fit$theta), se_first_order = unname(fit$se_first_order),
             residual_Q = unname(q), residual_df = as.integer(n - 1L),
             residual_SE = unname(rse), se_inflation = unname(inflation),
             se_default = unname(fit$se_first_order * inflation),
             stringsAsFactors = FALSE)
}

# Two optional JOINT covariance assumptions with the same default marginal SEs.
# Inputs: propagated sampling covariance V0 and aligned named inflation vector.
# Returns sampling_only, default_scaled_correlation, default_independent_excess.
# Neither sensitivity assumption follows from the separate mr_ivw fits, shared
# controls, or conditioning on exposure effects. Computing both does not select
# either as a final manuscript analysis.
ivw_working_covariances <- function(V0, inflation) {
  V0 <- .check_cov(V0, "First-order covariance")
  if (!is.numeric(inflation) || !identical(names(inflation), rownames(V0)) ||
      any(!is.finite(inflation)) || any(inflation < 1)) stop("Invalid named SE inflation.")
  # Sensitivity A specifies total residual covariance by scaling V0. It retains
  # the propagated MR correlation, not necessarily the outcome correlation.
  # Its difference from V0 need not be positive semidefinite, so this is not
  # generally an additive model with an independent extra-error component.
  scaled <- V0 * outer(inflation, inflation)
  # Sensitivity B retains V0 and adds a diagonal excess covariance at the MR
  # level. Zero covariance of that excess component across subtypes is a
  # scientific working assumption, not something estimated by this pipeline.
  independent_excess <- V0
  diag(independent_excess) <- diag(V0) * inflation^2
  list(sampling_only = V0, default_scaled_correlation = scaled,
       default_independent_excess = independent_excess)
}

# Input: actual retained/harmonized instruments, one row per protein/SNP/subtype.
# Omitted rows denote SNPs unavailable in that subtype; every subtype must remain.
# Effects in all rows must refer to effect_allele. This checks cross-subtype
# consistency, not the original exposure/outcome harmonization itself.
# Returns aligned bx, by/byse matrices, and the original validated input rows.
# This lower-level helper remains general; release-layer common-set checks
# apply specifically to the current five-subtype project.
.prepare_shared_input <- function(d, subtypes) {
  required <- c("protein", "snp", "subtype", "bx", "bxse", "by", "byse",
                "effect_allele", "other_allele")
  if (!is.data.frame(d) || !all(required %in% names(d)) || nrow(d) == 0L) {
    stop("Input needs nonempty rows and columns: ", paste(required, collapse = ", "))
  }
  for (key in c("protein", "snp", "subtype", "effect_allele", "other_allele")) {
    if (anyNA(d[[key]]) || any(trimws(d[[key]]) == "")) stop("Missing identifier: ", key)
    d[[key]] <- as.character(d[[key]])
  }
  if (length(unique(d$protein)) != 1L) stop("Prepare one protein at a time.")
  if (!setequal(unique(d$subtype), subtypes)) stop("Exactly the five BCAC subtypes are required.")
  if (anyDuplicated(d[c("snp", "subtype")])) stop("Duplicate SNP/subtype rows.")
  for (key in c("bx", "bxse", "by", "byse")) {
    if (!is.numeric(d[[key]]) || any(!is.finite(d[[key]]))) stop("Nonfinite/nonnumeric ", key)
  }
  if (any(d$bxse <= 0) || any(d$byse <= 0)) stop("All supplied SEs must be positive.")
  if (any(d$effect_allele == d$other_allele)) stop("Effect and other allele coincide.")
  for (g in split(d, d$snp)) {
    for (key in c("bx", "bxse", "effect_allele", "other_allele")) {
      if (length(unique(g[[key]])) != 1L) stop("Across-subtype mismatch for ", g$snp[1], ": ", key)
    }
  }
  ids <- sort(unique(d$snp))
  first <- match(ids, d$snp)
  bx <- setNames(d$bx[first], ids)
  by <- byse <- matrix(NA_real_, length(ids), length(subtypes), dimnames = list(ids, subtypes))
  at <- cbind(match(d$snp, ids), match(d$subtype, subtypes))
  by[at] <- d$by
  byse[at] <- d$byse
  list(bx = bx, by = by, byse = byse, input = d)
}

# Returns all assumptions without selecting whichever gives the smallest P.
# n_proteins controls the prespecified family size, even for an incomplete run.
# Input d contains one protein's actual SNP/subtype rows. Returns results
# (nine scenarios), diagnostics, all covariance details, aligned_input,
# counts, and assumptions. It does not read files or choose a result hierarchy.
analyse_shared_control_protein <- function(d, n_proteins = 12L) {
  counts <- bcac_subtype_counts()
  types <- counts$subtype
  x <- .prepare_shared_input(d, types)
  correlations <- list(
    effective_cases = shared_control_correlation(setNames(counts$effective_cases, types), 91477),
    observed_cases = shared_control_correlation(setNames(counts$observed_cases, types), 91477),
    independence_reference = diag(5))
  dimnames(correlations$independence_reference) <- list(types, types)
  results <- list(); details <- list(); diagnostics <- NULL
  for (scenario in names(correlations)) {
    Cy <- outcome_covariance_from_correlation(x$byse, correlations[[scenario]])
    fit <- ivw_covariance_first_order(x$bx, x$by, x$byse, Cy)
    disp <- ivw_default_dispersion(x$bx, x$by, x$byse, fit)
    inflation <- setNames(disp$se_inflation, types)
    Vs <- ivw_working_covariances(fit$V, inflation)
    tests <- lapply(Vs, function(V) subtype_heterogeneity(fit$theta, V, n_proteins))
    for (variance in names(tests)) {
      z <- tests[[variance]]
      results[[length(results) + 1L]] <- data.frame(
        protein = unique(d$protein), outcome_correlation = scenario,
        variance_assumption = variance, Q = z$Q, df = z$df,
        p_value = z$p_value, threshold = z$threshold, below_threshold = z$below_threshold,
        condition_number = z$condition_number, approximate_covariance = TRUE,
        stringsAsFactors = FALSE)
    }
    details[[scenario]] <- list(outcome_correlation = correlations[[scenario]],
      outcome_covariances = Cy, fit = fit, mr_covariances = Vs, tests = tests)
    diagnostics <- disp
  }
  diagnostics$protein <- unique(d$protein)
  list(results = do.call(rbind, results), diagnostics = diagnostics, details = details,
       aligned_input = x$input, counts = counts,
       assumptions = c("BCAC alone, 91477 identical shared controls",
         "Count-based correlations approximate EM outcome-error correlations",
         "First-order simple IVW, fixed exposure effects, negligible residual SNP LD",
         "MRE covariance models are optional sensitivity assumptions, not fitted joint covariance",
         "Chi-square calibration is asymptotic with plug-in dispersion and correlation"))
}

# Run after exporting the ACTUAL retained instruments from the original pipeline.
# New output directory required to preserve any previous results.
# tnbc_source is explicit to prevent accidental BCAC+CIMBA use.
# General file runner: supports partial proteins and differing available SNP
# sets with warnings/checks as documented. Use run_analysis.R for this
# project's strict 12-protein/common-SNP-set release. Returns invisibly a list
# of results, diagnostics, and per-protein fits; writes CSV/RDS/provenance files.
run_shared_control_pipeline <- function(input_csv, output_dir, tnbc_source,
                                        n_proteins = 12L, original_results_csv = NULL) {
  if (!identical(tnbc_source, "BCAC_only")) stop("This pipeline supports BCAC_only outcomes.")
  if (length(n_proteins) != 1L || !is.numeric(n_proteins) || !is.finite(n_proteins) ||
      n_proteins < 1 || n_proteins != as.integer(n_proteins)) stop("Invalid protein family size.")
  if (dir.exists(output_dir) || file.exists(output_dir)) stop("Choose a new output directory.")
  d <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(d) == 0L || !"protein" %in% names(d) || anyNA(d$protein) ||
      any(trimws(d$protein) == "")) stop("Nonempty input with valid protein IDs required.")
  if (length(unique(d$protein)) > n_proteins) stop("More proteins than the declared test family.")
  fits <- lapply(split(d, d$protein), analyse_shared_control_protein, n_proteins = n_proteins)
  results <- do.call(rbind, lapply(fits, `[[`, "results"))
  diagnostics <- do.call(rbind, lapply(fits, `[[`, "diagnostics"))
  rownames(results) <- rownames(diagnostics) <- NULL
  comparison <- NULL
  if (!is.null(original_results_csv)) {
    old <- read.csv(original_results_csv, stringsAsFactors = FALSE)
    needed <- c("protein", "subtype", "estimate", "se_default")
    if (!all(needed %in% names(old)) || anyNA(old[needed]) ||
        anyDuplicated(old[c("protein", "subtype")])) stop("Invalid original MR results file.")
    comparison <- merge(diagnostics[c("protein", "subtype", "estimate", "se_default")],
                        old[needed], by = c("protein", "subtype"), all = TRUE,
                        suffixes = c("_recomputed", "_original"))
    if (anyNA(comparison)) stop("Original and recomputed result keys must match.")
    for (key in c("estimate", "se_default")) {
      new <- comparison[[paste0(key, "_recomputed")]]
      old_value <- comparison[[paste0(key, "_original")]]
      if (!is.numeric(old_value) || any(!is.finite(old_value)) ||
          any(abs(new - old_value) > 1e-8 + 1e-6 * abs(old_value))) {
        stop("Existing MR ", key, " does not reproduce. Resolve inputs/settings before analysis.")
      }
    }
  }
  complete <- length(fits) == n_proteins
  if (!complete) warning("Partial run: the multiplicity denominator remains ", n_proteins, ".")
  if (is.null(comparison)) warning("Original MR estimates/SEs have not yet been independently compared.")
  dir.create(output_dir, recursive = TRUE)
  write.csv(results, file.path(output_dir, "heterogeneity_scenarios.csv"), row.names = FALSE)
  write.csv(diagnostics, file.path(output_dir, "mr_diagnostics.csv"), row.names = FALSE)
  write.csv(bcac_subtype_counts(), file.path(output_dir, "sample_sizes_used.csv"), row.names = FALSE)
  if (!is.null(comparison)) write.csv(comparison, file.path(output_dir, "original_MR_comparison.csv"), row.names = FALSE)
  saveRDS(fits, file.path(output_dir, "full_covariance_results.rds"))
  provenance <- c(paste("Input MD5:", tools::md5sum(input_csv)),
    paste("Run time:", format(Sys.time(), tz = "UTC")),
    paste("All declared proteins present:", complete),
    paste("Original MR estimates and SEs reproduced:", !is.null(comparison)),
    fits[[1]]$assumptions,
    "These scenarios do not bound all possible EM covariance structures.",
    capture.output(sessionInfo()))
  writeLines(provenance, file.path(output_dir, "provenance.txt"))
  invisible(list(results = results, diagnostics = diagnostics, fits = fits))
}
