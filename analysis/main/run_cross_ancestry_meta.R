# Purpose: Run fixed-effect cross-ancestry meta-analysis across EUR, EAS, and AFR MR results.

library(vroom)
library(tidyverse)
library(meta)
library(ACAT)

config <- list(
  eas_results_path = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/eas.csv",
  afr_results_path = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/afr.csv",
  eur_results_path = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/eur.csv",
  output_path = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/meta.csv"
)

data_EAS <- vroom(config$eas_results_path)
data_AFR <- vroom(config$afr_results_path)
data_EUR <- vroom(config$eur_results_path)

prepare_data <- function(data, suffix) {
  data %>%
    # Rename all columns except Protein
    rename_with(~paste0(., "_", suffix), -Protein)
}

data_EAS <- prepare_data(data_EAS, "EAS")
data_AFR <- prepare_data(data_AFR, "AFR")
data_EUR <- prepare_data(data_EUR, "EUR")

merged_data <- full_join(data_EAS, data_AFR, by = "Protein") %>%
  full_join(data_EUR, by = "Protein")

convert_to_numeric <- function(df, exclude_cols) {
  df %>%
    mutate(across(!all_of(exclude_cols), ~{
      numeric_values <- suppressWarnings(as.numeric(as.character(.)))
      numeric_values
    }))
}

exclude_cols <- c("Protein", 
                  grep("ivs", names(merged_data), value = TRUE),
                  grep("ACAT|Bonferroni|fdr|duplicate", names(merged_data), value = TRUE))
data <- convert_to_numeric(merged_data, exclude_cols)

perform_meta_analysis <- function(data, protein, method, p_threshold) {
  ancestries <- c("EAS", "EUR", "AFR")
  estimates <- c()
  ses <- c()
  study_labels <- c()
  
  for (ancestry in ancestries) {
    estimate_col <- paste0(p_threshold, ".Estimate_", method, "_", ancestry)
    se_col <- paste0(p_threshold, ".SE_", method, "_", ancestry)
    
    if (all(c(estimate_col, se_col) %in% names(data))) {
      estimate <- data[[estimate_col]][data$Protein == protein]
      se <- data[[se_col]][data$Protein == protein]
      
      if (length(estimate) > 0 && length(se) > 0) {
        if (!is.na(estimate) && !is.na(se) && se > 0) {
          estimates <- c(estimates, estimate)
          ses <- c(ses, se)
          study_labels <- c(study_labels, ancestry)
        }
      }
    }
  }
  
  if (length(estimates) >= 2) {
    tryCatch({
      meta_result <- metagen(
        TE = estimates,
        seTE = ses,
        studlab = study_labels,
        comb.fixed = TRUE,
        comb.random = FALSE
      )
      return(meta_result)
    }, error = function(e) {
      return(NULL)
    })
  } else {
    return(NULL)
  }
}

methods <- c("ivw", "median", "pivw", "raps")
p_thresholds <- c("5e-08", "5e-07", "5e-06")

meta_results <- list()
for (protein in unique(data$Protein)) {
  for (method in methods) {
    for (p_threshold in p_thresholds) {
      key <- paste(protein, method, p_threshold, sep = "_")
      result <- perform_meta_analysis(data, protein, method, p_threshold)
      if (!is.null(result)) {
        meta_results[[key]] <- result
      }
    }
  }
}

summarize_meta_results <- function(meta_results) {
  summary_data <- data.frame()
  
  for (key in names(meta_results)) {
    result <- meta_results[[key]]
    
    if (!is.null(result)) {
      parts <- strsplit(key, "_")[[1]]
      method <- parts[length(parts) - 1]
      p_threshold <- parts[length(parts)]
      protein <- paste(parts[1:(length(parts) - 2)], collapse = "_")
      
      summary_row <- data.frame(
        Protein = protein,
        Method = method,
        P_Threshold = p_threshold,
        TE_fixed = result$TE.fixed,
        seTE_fixed = result$seTE.fixed,
        Lower_CI_fixed = result$lower.fixed,
        Upper_CI_fixed = result$upper.fixed,
        P_Value_fixed = result$pval.fixed,
        N_studies = result$k,
        Q = result$Q,
        Q_pval = result$pval.Q,
        I2 = result$I2,
        tau2 = result$tau2
      )
      
      summary_data <- rbind(summary_data, summary_row)
    }
  }
  
  return(summary_data)
}

meta_summary <- summarize_meta_results(meta_results)

wide_summary <- meta_summary %>%
  pivot_wider(
    id_cols = Protein,
    names_from = c(Method, P_Threshold),
    values_from = c(TE_fixed, seTE_fixed, Lower_CI_fixed, Upper_CI_fixed, P_Value_fixed,
                    N_studies, Q, Q_pval, I2, tau2), 
    names_sep = "_"
  )

pvalue_columns <- grep("P_Value_fixed", names(wide_summary), value = TRUE)
genes_pvalues <- wide_summary[, c("Protein", pvalue_columns)]

genes_pvalues$ACAT_pvalue <- apply(genes_pvalues[, pvalue_columns], 1, function(row) {
  pvals <- as.numeric(row)
  pvals <- pvals[!is.na(pvals)]
  if (length(pvals) == 0) {
    return(NA)
  } else {
    return(ACAT::ACAT(pvals))
  }
})

genes_pvalues$Bonferroni_pvalue <- p.adjust(genes_pvalues$ACAT_pvalue, method = "bonferroni")
genes_pvalues$BH_pvalue <- p.adjust(genes_pvalues$ACAT_pvalue, method = "BH")

q_pval_columns <- grep("Q_pval", names(wide_summary), value = TRUE)
genes_q_pvalues <- wide_summary[, c("Protein", q_pval_columns)]

genes_q_pvalues$ACAT_Q_pval <- apply(genes_q_pvalues[, q_pval_columns], 1, function(row) {
  pvals <- as.numeric(row)
  pvals <- pvals[!is.na(pvals)]
  if (length(pvals) == 0) {
    return(NA)
  } else {
    return(ACAT::ACAT(pvals))
  }
})

genes_q_pvalues$Bonferroni_Q_pval <- p.adjust(genes_q_pvalues$ACAT_Q_pval, method = "bonferroni")
genes_q_pvalues$BH_Q_pval <- p.adjust(genes_q_pvalues$ACAT_Q_pval, method = "BH")

combined_data <- wide_summary %>%
  left_join(genes_pvalues[, c("Protein", "ACAT_pvalue", "Bonferroni_pvalue", "BH_pvalue")], by = "Protein") %>%
  left_join(genes_q_pvalues[, c("Protein", "ACAT_Q_pval", "Bonferroni_Q_pval", "BH_Q_pval")], by = "Protein") %>%
  arrange(ACAT_pvalue)

head(combined_data)
names(combined_data)
meta_results_table <- combined_data %>%
  arrange(as.numeric(`P_Value_fixed_ivw_5e-08`))

write.csv(meta_results_table, config$output_path, row.names = FALSE)
