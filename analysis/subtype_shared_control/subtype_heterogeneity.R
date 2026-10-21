# Covariance-aware comparison of subtype MR estimates.
# Analyst template, 2026-10-13. No study data or fitted covariance included.
# See ANALYST_GUIDE.html or ANALYST_GUIDE.md for assumptions and integration.
# Uses base R only. Source this file; it performs no analysis automatically.

# Validate a finite symmetric positive-semidefinite covariance matrix.
# Inputs: V (square numeric matrix), diagnostic label, relative tolerance.
# Returns the symmetric matrix; stops on invalid values instead of repairing
# a substantively invalid covariance. Positive diagonal entries are required.
.check_cov <- function(V, label, tol = 1e-8) {
  if (!is.matrix(V) || !is.numeric(V) || nrow(V) != ncol(V) ||
      any(!is.finite(V)) || any(diag(V) <= 0)) {
    stop(label, " must be a finite square covariance matrix with positive diagonal.")
  }
  scale <- max(diag(V))
  if (max(abs(V - t(V))) > tol * scale) stop(label, " is not symmetric.")
  V <- (V + t(V)) / 2
  if (min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) < -tol * scale) {
    stop(label, " is not positive semidefinite. Investigate; do not silently repair.")
  }
  V
}

# theta: named vector of K subtype MR log-OR estimates, for ONE protein.
# V: named K x K sampling covariance matrix of these SAME estimates.
# Names must be in exactly the same order. There is no automatic reordering.
# Tests H0: theta_1 = ... = theta_K. Does not test whether all are zero.
# n_proteins is the prespecified multiplicity denominator, not a count of
# significant proteins. Returns Q, df, p_value, threshold, below_threshold,
# contrasts, contrast_covariance, mr_covariance, and condition_number.
subtype_heterogeneity <- function(theta, V, n_proteins = 12L) {
  if (!is.numeric(theta) || length(theta) < 2L || any(!is.finite(theta)) ||
      is.null(names(theta)) || anyDuplicated(names(theta)) ||
      any(names(theta) == "")) stop("theta must be a finite named vector.")
  if (!identical(rownames(V), names(theta)) ||
      !identical(colnames(V), names(theta))) stop("V names/order must match theta.")
  if (length(n_proteins) != 1L || !is.finite(n_proteins) || n_proteins < 1L) {
    stop("n_proteins must be positive.")
  }
  V <- .check_cov(V, "MR covariance")
  K <- length(theta)
  C <- cbind(-1, diag(K - 1L))
  colnames(C) <- names(theta)
  rownames(C) <- paste(names(theta)[-1L], "minus", names(theta)[1L])
  difference <- drop(C %*% theta)
  W <- C %*% V %*% t(C)
  ev <- eigen(W, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-10 * max(ev)) {
    stop("Contrast covariance is singular or ill-conditioned; inspect the model.")
  }
  U <- chol(W)
  whitened <- forwardsolve(t(U), difference)
  Q <- sum(whitened^2)
  df <- K - 1L
  p <- pchisq(Q, df = df, lower.tail = FALSE)
  list(Q = Q, df = df, p_value = p, threshold = 0.05 / n_proteins,
       below_threshold = p < 0.05 / n_proteins,
       contrasts = difference, contrast_covariance = W,
       mr_covariance = V, condition_number = max(ev) / min(ev))
}

# Build covariance for unpenalized, first-order IVW with independent SNPs.
# CONDITIONAL on the estimated exposure effects (NOME approximation).
# Does NOT implement multiplicative random-effects, robust/penalized weights,
# exposure-effect uncertainty, exposure/outcome overlap, or residual SNP LD.
#
# bx: named vector (J SNPs) of aligned SNP-protein effects.
# by, byse: J x K matrices, rownames = SNP IDs, colnames = subtype names.
# Missing instruments may be NA in BOTH by and byse. No silent filtering.
# Cy: named list, one K x K SNP-outcome covariance matrix per SNP.
#     Each matrix uses the same subtype order as by. Its diagonal must equal
#     byse[j, ]^2 for the subtypes available at SNP j. Unused entries may be NA.
#     Cy must describe the FINAL outcome GWAS, after its cohort meta-analysis.
# Returns theta (K-vector), V (K x K sampling covariance), se_first_order,
# coefficients (J x K matrix), n_snps (K-vector), and a model label.
# The general engine permits missing SNP sets; the project release runner
# separately requires identical retained sets across its five subtypes.
ivw_covariance_first_order <- function(bx, by, byse, Cy) {
  if (!is.matrix(by) || !is.matrix(byse) || !is.numeric(by) ||
      !is.numeric(byse) || !identical(dimnames(by), dimnames(byse)) ||
      !identical(dim(by), dim(byse)) || ncol(by) < 2L) {
    stop("by and byse must be matching numeric J x K matrices, K >= 2.")
  }
  ids <- rownames(by)
  subtypes <- colnames(by)
  if (is.null(ids) || is.null(subtypes) || anyDuplicated(ids) ||
      anyDuplicated(subtypes) || any(ids == "") || any(subtypes == "")) {
    stop("Unique nonempty SNP and subtype names are required.")
  }
  if (!is.numeric(bx) || !identical(names(bx), ids) || any(!is.finite(bx))) {
    stop("bx must be finite and its names/order must match the SNP rows.")
  }
  if (!is.list(Cy) || is.null(names(Cy)) || anyDuplicated(names(Cy)) ||
      !all(ids %in% names(Cy))) stop("Cy must include every SNP by name.")
  if (any(is.na(by) != is.na(byse))) stop("Missing by and byse must coincide.")
  available <- !is.na(by)
  if (any(!is.finite(by[available])) || any(!is.finite(byse[available])) ||
      any(byse[available] <= 0)) stop("Available effects/SEs must be finite; SE > 0.")
  if (any(colSums(available) == 0L)) stop("Each subtype needs at least one SNP.")
  J <- nrow(by)
  K <- ncol(by)
  precision <- matrix(0, J, K, dimnames = dimnames(by))
  precision[available] <- 1 / byse[available]^2
  D <- colSums(precision * bx^2)
  if (any(!is.finite(D)) || any(D <= 0)) stop("Nonpositive IVW information.")
  A <- sweep(precision * bx, 2L, D, "/")
  outcome <- by
  outcome[!available] <- 0
  theta <- colSums(A * outcome)
  V <- matrix(0, K, K, dimnames = list(subtypes, subtypes))
  for (j in seq_len(J)) {
    M <- Cy[[ids[j]]]
    if (!is.matrix(M) || !identical(rownames(M), subtypes) ||
        !identical(colnames(M), subtypes)) stop("Cy order mismatch for ", ids[j])
    use <- which(available[j, ])
    if (length(use) == 0L) next
    M <- .check_cov(M[use, use, drop = FALSE], paste("Cy for", ids[j]))
    target <- byse[j, use]^2
    if (max(abs(diag(M) / target - 1)) > 1e-6) {
      stop("Cy diagonal differs from supplied outcome SEs for ", ids[j])
    }
    V[use, use] <- V[use, use, drop = FALSE] +
      outer(A[j, use], A[j, use]) * M
  }
  V <- .check_cov(V, "Propagated MR covariance")
  if (max(abs(diag(V) * D - 1)) > 1e-6) stop("Internal diagonal check failed.")
  list(theta = theta, V = V, se_first_order = sqrt(diag(V)),
       coefficients = A, n_snps = colSums(available),
       model = "Conditional first-order IVW, independent SNPs, fixed-effect variance")
}

# OPTIONAL approximation: construct Cy from a validated outcome sampling-error
# correlation matrix R_outcome. This is NOT genetic correlation, correlation
# of the 12 proteins' estimates, or automatically the MR-estimate correlation.
# It assumes the same outcome-error correlation at each SNP. No covariance is
# estimated here; the analyst must supply and justify R_outcome.
# Returns a SNP-named list of K x K covariance matrices. Missing outcome SEs
# remain missing, for the general engine's explicit availability handling.
outcome_covariance_from_correlation <- function(byse, R_outcome) {
  if (!is.matrix(byse) || !is.numeric(byse) || is.null(rownames(byse)) ||
      is.null(colnames(byse)) || anyDuplicated(rownames(byse)) ||
      anyDuplicated(colnames(byse))) stop("Named, numeric byse matrix required.")
  if (any(!is.finite(byse[!is.na(byse)])) || any(byse[!is.na(byse)] <= 0)) {
    stop("Available outcome SEs must be positive and finite.")
  }
  if (!identical(rownames(R_outcome), colnames(byse)) ||
      !identical(colnames(R_outcome), colnames(byse))) stop("R_outcome order mismatch.")
  R_outcome <- .check_cov(R_outcome, "Outcome error correlation")
  if (max(abs(diag(R_outcome) - 1)) > 1e-8) stop("R_outcome diagonal must equal one.")
  result <- lapply(seq_len(nrow(byse)), function(j) {
    M <- outer(byse[j, ], byse[j, ]) * R_outcome
    dimnames(M) <- dimnames(R_outcome)
    M
  })
  names(result) <- rownames(byse)
  result
}

# Usage after preparing and validating real inputs:
# source("subtype_heterogeneity.R")
# fit <- ivw_covariance_first_order(bx, by, byse, Cy)
# result <- subtype_heterogeneity(fit$theta, fit$V, n_proteins = 12)
# Save result, fit, covariance provenance, input IDs and sessionInfo().
# If a covariance appropriate to a different estimator is available, call
# subtype_heterogeneity(existing_named_log_ORs, matching_MR_covariance) directly.
