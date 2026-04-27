suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(coloc)
  library(ggplot2)
  library(writexl)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Provide the protein index as the first command-line argument.")
}

protein_index <- suppressWarnings(as.integer(args[1]))
if (is.na(protein_index) || protein_index < 1) {
  stop("The protein index must be a positive integer.")
}

slurm_job_id <- Sys.getenv("SLURM_JOB_ID", unset = "")
if (nzchar(slurm_job_id)) {
  temp_dir <- file.path("/lscratch", slurm_job_id, paste0("temp_", protein_index))
} else {
  temp_dir <- file.path(tempdir(), paste0("temp_", protein_index))
}
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
message("Temporary directory: ", temp_dir)

config <- list(
  all_rsids_rdata = "/data/BB_Bioinformatics/DG/MR_bc/all_rsids.RData",
  outcome_data_path = "/data/BB_Bioinformatics/ProjectData/breast_cancer_sum_data/BCAC_EUR/bc_summary_gwas.txt",
  tar_files_dir = "/data/BB_Bioinformatics/ProjectData/UKB_protein_sumstat/UKB-PPP_pGWAS_summary_statistics/European_discovery",
  results_dir = "/data/BB_Bioinformatics/DG/MR_bc/2025_updated/coloc/results37",
  n_protein_discovery = 34557,
  n_bcac = 197755,
  bcac_case_fraction = 106278 / 197755,
  region_width = 1.5e5
)

# Load required data
load(config$all_rsids_rdata)
outcome_data <- fread(config$outcome_data_path)

setDT(outcome_data)
outcome_data[, c("allele1", "allele2") := .(
    pmin(Effect.Meta, Baseline.Meta),
    pmax(Effect.Meta, Baseline.Meta)
  )]
outcome_data[, c("chr", "pos") := list(chr.Onco, Position.Onco)]

# Set directories
plots_dir <- file.path(config$results_dir, "plots")

# Create directories if they don't exist
if (!dir.exists(config$results_dir)) dir.create(config$results_dir, recursive = TRUE)
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Updated protein data for the 12 target proteins
all_protein_data <- list(
  list(
    protein = "CASP8_OID21414",
    snps = c("rs527566839", "rs533429656", "rs80201352", "rs34841024", "rs138352040", "rs575593318"),
    search_name = "CASP8_Q14790_OID21414_v1_Oncology",
    protein_symbol = "CASP8"
  ),
  list(
    protein = "USP28_OID31084",
    snps = c("rs58682649", "rs2459974", "rs4288784"),
    search_name = "USP28_Q96RU2_OID31084_v1_Neurology_II",
    protein_symbol = "USP28"
  ),
  list(
    protein = "RALB_OID30567",
    snps = c("rs148102406"),
    search_name = "RALB_P11234_OID30567_v1_Inflammation_II",
    protein_symbol = "RALB"
  ),
  list(
    protein = "PARK7_OID21160",
    snps = c("rs3766606", "rs6684255"),
    search_name = "PARK7_Q99497_OID21160_v1_Neurology",
    protein_symbol = "PARK7"
  ),
  list(
    protein = "SCAMP3_OID21332",
    snps = c("rs147340687", "rs548079062", "rs1142287"),
    search_name = "SCAMP3_O14828_OID21332_v1_Oncology",
    protein_symbol = "SCAMP3"
  ),
  list(
    protein = "GCLM_OID30059",
    snps = c("rs376647070", "rs12563495"),
    search_name = "GCLM_P48507_OID30059_v1_Cardiometabolic_II",
    protein_symbol = "GCLM"
  ),
  list(
    protein = "LRRC37A2_OID31440",
    snps = c("rs147691494", "rs141461397", "rs149448680", "rs241043", "rs118084908", 
             "rs34804414", "rs7221124", "rs55938136", "rs116876049", "rs4792886", 
             "rs140161293", "rs7216893", "rs554998846", "rs150324371", "rs273537", 
             "rs147405729", "rs145992408", "rs57834988", "rs73987093", "rs182054669", 
             "rs146330795", "rs62073767", "rs143416212", "rs141711906", "rs113552917", 
             "rs111527969", "rs74969489", "rs59850765", "rs186985777", "rs62073975"),
    search_name = "LRRC37A2_A6NM11_OID31440_v1_Oncology_II",
    protein_symbol = "LRRC37A2"
  ),
  list(
    protein = "ANXA4_OID20123",
    snps = c("rs56398743"),
    search_name = "ANXA4_P09525_OID20123_v1_Cardiometabolic",
    protein_symbol = "ANXA4"
  ),
  list(
    protein = "DNPH1_OID20744",
    snps = c("rs72857647", "rs115736156", "rs139999734", "rs3793008", "rs114326925",
             "rs77321231", "rs111607015", "rs553448379", "rs528833950", 
             "rs833063", "rs183381951"),
    search_name = "DNPH1_O43598_OID20744_v1_Inflammation",
    protein_symbol = "DNPH1"
  ),
  list(
    protein = "LRRC25_OID21333",
    snps = c("rs10424779", "rs537679857", "rs4808787", "rs34666550", "rs3859570", 
             "rs183423809", "rs181004295", "rs148140604", "rs118113459", "rs560350712", 
             "rs145618160", "rs111365980", "rs115551853", "rs11881338", "rs113219486"),
    search_name = "LRRC25_Q8N386_OID21333_v1_Oncology",
    protein_symbol = "LRRC25"
  ),
  list(
    protein = "ADM_OID21467",
    snps = c("rs138669902", "rs9332432", "rs7926490", "rs1349326"),
    search_name = "ADM_P35318_OID21467_v1_Oncology",
    protein_symbol = "ADM"
  ),
  list(
    protein = "RSPO3_OID21415",
    snps = c("rs527548446", "rs74949497", "rs142371601", "rs117937804", "rs72982840", "rs1936800"),
    search_name = "RSPO3_Q9BXY4_OID21415_v1_Oncology",
    protein_symbol = "RSPO3"
  ),
  list(
    protein = "TLR1_OID30532",
    snps = c("rs148892047", "rs200817536", "rs17499355", "rs115354241", "rs10213295", 
             "rs77816728", "rs57978465", "rs141192258", "rs7659947", "rs73234989", 
             "rs190184314", "rs9990980", "rs145622274", "rs200127224", "rs6531685", 
             "rs2130298", "rs200061843", "rs115562334"),
    search_name = "TLR1_Q15399_OID30532_v1_Inflammation_II",
    protein_symbol = "TLR1"
  ),
  list(
    protein = "FES_OID21207",
    snps = c("rs1894401"),
    search_name = "FES_P07332_OID21207_v1_Oncology",
    protein_symbol = "FES"
  ),
  list(
    protein = "TNFRSF9_OID20985",
    snps = c("rs78752371", "rs2493214"),
    search_name = "TNFRSF9_Q07011_OID20985_v1_Neurology",
    protein_symbol = "TNFRSF9"
  ),
  list(
    protein = "DDX58_OID21226",
    snps = c("rs140262942", "rs7867666", "rs79108616", "rs56327602", "rs141032208", 
             "rs148040790", "rs1133071", "rs10971004", "rs557605117", "rs593700", 
             "rs148019136", "rs60131401"),
    search_name = "DDX58_O95786_OID21226_v1_Oncology",
    protein_symbol = "DDX58"
  )
)

# Select protein based on index
if (protein_index < 1 || protein_index > length(all_protein_data)) {
  stop("Invalid protein index. Must be between 1 and ", length(all_protein_data))
}
protein_data <- all_protein_data[[protein_index]]
cat("Processing protein:", protein_data$protein_symbol, "\n")

setup_protein_data <- function(protein) {
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
  }
  
  tar_files <- list.files(config$tar_files_dir, pattern = "\\.tar$", full.names = TRUE)
  index <- grep(protein$search_name, tar_files, fixed = TRUE, value = FALSE)
  
  if (length(index) == 0) {
    stop(paste("No tar file found for protein:", protein$search_name))
  }
  if (length(index) > 1) {
    stop(paste("Multiple tar files matched protein:", protein$search_name))
  }
  
  data_dir <- file.path(temp_dir, protein$search_name)
  
  # Create protein-specific directory if it doesn't exist
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
    tar_file <- tar_files[index]
    message("Extracting tar file: ", basename(tar_file))
    tryCatch({
      untar(tar_file, exdir = data_dir)
    }, error = function(e) {
      stop(paste("Error extracting tar file:", e$message))
    })
  }
  
  data_files <- list.files(data_dir, pattern = "\\.gz$", full.names = TRUE, recursive = TRUE)
  data_files <- data_files[!grepl("chrX", data_files)]
  
  list_of_data_frames <- lapply(data_files, function(f) {
    tryCatch(fread(f), error = function(e) NULL)
  })
  list_of_data_frames <- Filter(Negate(is.null), list_of_data_frames)
  if (length(list_of_data_frames) == 0) {
    stop(paste("No readable summary-statistic files found for protein:", protein$search_name))
  }
  
  combined_data <- rbindlist(list_of_data_frames, fill = TRUE)
  
  combined_data[, c("allele1", "allele2") := .(
    pmin(ALLELE0, ALLELE1),
    pmax(ALLELE0, ALLELE1)
  )]
  
  combined_data[, PVALUE := 10^(-LOG10P)]
  
  # Check for duplicate IDs before merge
  if (any(duplicated(combined_data$ID))) {
    message("Found ", sum(duplicated(combined_data$ID)), " duplicate IDs in combined_data")
    combined_data <- unique(combined_data)
    message("After removing duplicates: ", nrow(combined_data), " rows remain")
  }
  
  # Add check for duplicates in all_rsids
  if (any(duplicated(all_rsids$ID))) {
    message("Found ", sum(duplicated(all_rsids$ID)), " duplicate IDs in all_rsids")
  }
  
  # Merge to get rsids and handle duplicates that may arise
  combined_data <- merge(combined_data, all_rsids[, .(ID, rsid, POS19)], by = "ID", all.x = FALSE)
  
  # Check for duplicate rsids after merge
  if (any(duplicated(combined_data$rsid))) {
    message("Found ", sum(duplicated(combined_data$rsid)), " duplicate rsids after merge")
    
    # Example of duplicates
    dups <- combined_data[rsid %in% combined_data$rsid[duplicated(combined_data$rsid)]][order(rsid)]
    message("Example of duplicates:")
    print(head(dups[, .(rsid, ID, ALLELE0, ALLELE1, PVALUE)]))
    
    # Sort by p-value within each rsid and keep the most significant variant
    setorder(combined_data, rsid, PVALUE)
    combined_data <- combined_data[!duplicated(combined_data$rsid)]
    message("After removing duplicates: ", nrow(combined_data), " rows remain")
  }
  
  combined_data[, `:=`(
    chr = as.numeric(CHROM),
    pos = as.numeric(POS19)
  )]
  
  return(combined_data)
}

process_single_snp <- function(protein, single_snp, combined_data) {
  tryCatch({
    all_rsids[, CHROM := as.numeric(gsub(":.*$", "", ID))]
    snp_position <- all_rsids[rsid == single_snp, .(chr = CHROM, pos = POS19)]
    
    if (nrow(snp_position) == 0) {
      message("SNP ", single_snp, " not found in all_rsids database")
      return(NULL)
    }
    
    region_bounds <- as.data.table(snp_position)[, .(
      start = pos - config$region_width,
      end = pos + config$region_width
    ), by = chr]
    
    narrow_data <- combined_data[chr %in% region_bounds$chr &
                                 pos >= min(region_bounds$start) &
                                 pos <= max(region_bounds$end)]
    
    if (nrow(narrow_data) == 0) {
      message("No data found in the colocalization region for SNP ", single_snp)
      return(NULL)
    }
    
    merged_data <- merge(narrow_data, outcome_data,
                         by = c("chr", "pos", "allele1", "allele2"),
                         all.x = FALSE)
    
    if (nrow(merged_data) == 0) {
      message("No matching data after merging with outcome data for SNP ", single_snp)
      return(NULL)
    }
    
    # Check for duplicate rsids in merged data
    if (any(duplicated(merged_data$rsid))) {
      message("Found duplicate rsids in merged data for SNP ", single_snp)
      dup_rsids <- merged_data$rsid[duplicated(merged_data$rsid)]
      message("Found ", length(unique(dup_rsids)), " rsids with duplicates")
      
      # Create a log of duplicates before removing them
      dup_log <- merged_data[rsid %in% dup_rsids]
      dup_log <- dup_log[order(rsid, PVALUE)]
      
      # Print example of duplicates (up to 3)
      for (rs in unique(dup_rsids)[1:min(3, length(unique(dup_rsids)))]) {
        message("Example duplicate for rsid ", rs, ":")
        print(merged_data[rsid == rs, .(rsid, chr, pos, allele1, allele2, ALLELE0, ALLELE1, BETA, SE, PVALUE)])
      }
      
      # Sort by p-value and keep most significant
      setorder(merged_data, rsid, PVALUE)
      merged_data <- merged_data[!duplicated(merged_data$rsid)]
      
      message("After deduplication: ", nrow(merged_data), " SNPs remain")
    }
    
    merged_data[, `:=`(
      swap = (ALLELE1 != Effect.Meta),
      ALLELE1_new = ifelse(ALLELE1 != Effect.Meta, ALLELE0, ALLELE1),
      ALLELE0_new = ifelse(ALLELE1 != Effect.Meta, ALLELE1, ALLELE0),
      BETA = ifelse(ALLELE1 != Effect.Meta, -BETA, BETA),
      A1FREQ = ifelse(ALLELE1 != Effect.Meta, 1 - A1FREQ, A1FREQ)
    )]
    
    # Remove old allele columns and rename new ones
    merged_data[, c("ALLELE1", "ALLELE0") := NULL]
    setnames(merged_data, c("ALLELE1_new", "ALLELE0_new"), c("ALLELE1", "ALLELE0"))
    
    
    merged_data$PVALUE <- 10^(-merged_data$LOG10P)
    merged_data <- merged_data[A1FREQ >= 0.01 & A1FREQ <= 0.99 & INFO > 0.3]
    
    if (nrow(merged_data) == 0) {
      message("No data remaining after filtering for SNP ", single_snp)
      return(NULL)
    }
    
    dataset1 <- list(
      beta = merged_data$BETA,
      varbeta = merged_data$SE^2,
      type = "quant",
      snp = merged_data$rsid,
      N = config$n_protein_discovery,
      MAF = merged_data$A1FREQ,
      sdY = 1
    )
    
    dataset2 <- list(
      beta = merged_data$Beta.meta,
      varbeta = merged_data$sdE.meta^2,
      type = "cc",
      snp = merged_data$rsid,
      N = config$n_bcac,
      s = config$bcac_case_fraction,
      sdY = 1
    )
    
    coloc_results <- coloc.abf(dataset1, dataset2)
    
    # Extract all SNP information including posterior probabilities
    snp_stats <- data.frame(
      rsid = merged_data$rsid,
      protein_beta = merged_data$BETA,
      protein_pval = merged_data$PVALUE,
      bc_beta = merged_data$Beta.meta,
      bc_pval = 2*pnorm(-abs(merged_data$Beta.meta/merged_data$sdE.meta)),
      SNP_PP_H4 = coloc_results$results$SNP.PP.H4,
      chr = merged_data$chr,
      pos = merged_data$pos
    )
    
    snp_stats$joint_signal <- -log10(snp_stats$protein_pval) * -log10(snp_stats$bc_pval)
    snp_stats <- snp_stats[order(-snp_stats$SNP_PP_H4),]
    
    region_info <- sprintf("chr%s:%s-%s",
                         region_bounds$chr,
                         format(region_bounds$start, scientific=FALSE),
                         format(region_bounds$end, scientific=FALSE))
    
    list(
      results = coloc_results,
      protein = protein$protein,
      protein_symbol = protein$protein_symbol,
      analyzed_snp = single_snp,
      region = region_info,
      n_snps = nrow(merged_data),
      snp_stats = snp_stats, 
      iv_snps = protein$snps
    )
  }, error = function(e) {
    message("Error processing SNP ", single_snp, " for protein ", protein$protein, ": ", e$message)
    return(NULL)
  })
}

process_significant_results <- function(all_results) {
  significant_results <- list()
  
  for (result_name in names(all_results)) {
    result <- all_results[[result_name]]
    
    if (!is.null(result) && result$results$summary[6] >= 0.8) {  # H4 >= 0.8
      tryCatch({
        # Get SNPs with PP.H4 > 0.8
        highly_significant_snps <- subset(result$snp_stats, SNP_PP_H4 >= 0.8)
        
        # Store whether we have highly significant SNPs
        has_high_sig_snps <- nrow(highly_significant_snps) > 0
        if (!has_high_sig_snps) {
          message("No highly significant SNPs found for ", result_name)
        }
        
        # Calculate 95% credible set
        ordered_idx <- order(result$snp_stats$SNP_PP_H4, decreasing = TRUE)
        cumsum_pp <- cumsum(result$snp_stats$SNP_PP_H4[ordered_idx])
        
        if (length(cumsum_pp) > 0 && any(cumsum_pp > 0.95)) {
          credible_set_idx <- ordered_idx[1:which(cumsum_pp > 0.95)[1]]
          credible_set_snps <- result$snp_stats[credible_set_idx, ]
          credible_set_snps$in_credible_set <- TRUE
          
          if (has_high_sig_snps) {
            highly_significant_snps$in_credible_set <-
              highly_significant_snps$rsid %in% credible_set_snps$rsid
          }
        } else {
          message("Couldn't reach 95% credible set for ", result_name, ", using all SNPs")
          credible_set_snps <- result$snp_stats
          credible_set_snps$in_credible_set <- TRUE
          if (has_high_sig_snps) {
            highly_significant_snps$in_credible_set <- TRUE
          }
        }
        
        # Add query SNP information columns
        if (has_high_sig_snps) {
          highly_significant_snps$query_snp <- result$analyzed_snp
        }
        credible_set_snps$query_snp <- result$analyzed_snp
        
        # ---- New plotting logic ----
        # Identify the top SNP by SNP_PP_H4
        top_h4_snp <- result$snp_stats$rsid[which.max(result$snp_stats$SNP_PP_H4)]
        
        plot_data <- data.frame(
          SNP = result$snp_stats$rsid,
          protein_p = -log10(result$snp_stats$protein_pval),
          bc_p = -log10(result$snp_stats$bc_pval),
          SNP_PP_H4 = result$snp_stats$SNP_PP_H4
        )
        plot_data$is_query   <- plot_data$SNP == result$analyzed_snp
        plot_data$is_top_h4  <- plot_data$SNP == top_h4_snp
        plot_data$is_iv      <- plot_data$SNP %in% result$iv_snps  # requires change (1) you already added
        
        # Build plot with layer order: background -> IV circles -> diamonds on top
        coloc_plot <- ggplot(plot_data, aes(x = protein_p, y = bc_p)) +
          # Background: all SNPs in light gray
          geom_point(color = "gray", alpha = 0.5) +
          
          # IV SNPs as blue circles (exclude those that will be diamonds)
          geom_point(
            data = subset(plot_data, is_iv & !is_query & !is_top_h4),
            shape = 16, size = 2, color = "blue"
          ) +
          
          # Query SNP as blue diamond (but not if it's also top H4)
          geom_point(
            data = subset(plot_data, is_query & !is_top_h4),
            shape = 18, size = 4, color = "blue"
          ) +
          
          # Top SNP by SNP_PP_H4 as red diamond (takes precedence)
          geom_point(
            data = subset(plot_data, is_top_h4),
            shape = 18, size = 4, color = "red"
          ) +
          
          # Optional labels: label the two diamonds (query + top H4)
          geom_text(
            data = subset(plot_data, is_query | is_top_h4),
            aes(label = SNP),
            vjust = -0.8, size = 3
          ) +
          
          labs(
            title = paste0("Colocalization Plot for ", result$protein_symbol),
            subtitle = paste0(
              "Query SNP: ", result$analyzed_snp,
              "   |   Top H4 SNP: ", top_h4_snp,
              "\nRegion: ", result$region,
              "   |   PP.H4 = ", round(result$results$summary[6], 3)
            ),
            x = "-log10(p-value) for Association with Plasma Protein Level",
            y = "-log10(p-value) for Association with Breast Cancer"
          ) +
          theme_bw() +
          
          # Simple legend explanation (annotation text)
          annotate(
            "text",
            x = max(plot_data$protein_p, na.rm = TRUE) * 0.98,
            y = max(plot_data$bc_p, na.rm = TRUE) * 0.98,
            hjust = 1,
            size = 3.2,
            label = paste(
              "Blue circle: IV SNPs\n",
              "Blue diamond: Query SNP\n",
              "Red diamond: Top SNP by SNP_PP_H4"
            )
          )
        # ---- end new plotting logic ----
        
        # Store results
        significant_results[[result_name]] <- list(
          protein = result$protein,
          protein_symbol = result$protein_symbol,
          region = result$region,
          H4 = result$results$summary[6],
          H3 = result$results$summary[5],
          highly_significant_snps = if (has_high_sig_snps) highly_significant_snps else NULL,
          credible_set_snps = credible_set_snps,
          plot = coloc_plot,
          analyzed_snp = result$analyzed_snp
        )
      }, error = function(e) {
        message("Error processing significant result for ", result_name, ": ", e$message)
        print(traceback())  # Show detailed stack trace
      })
    }
  }
  
  return(significant_results)
}

create_summary_tables <- function(all_results, significant_results) {
  # Full results table
  full_results <- do.call(rbind, lapply(names(all_results), function(name) {
    result <- all_results[[name]]
    if(is.null(result)) return(NULL)
    
    tryCatch({
      data.frame(
        protein = result$protein,
        protein_symbol = result$protein_symbol,
        analyzed_snp = result$analyzed_snp,
        region = result$region,
        n_snps = result$n_snps,
        PPH4 = result$results$summary[6],
        PPH3 = result$results$summary[5],
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message("Error creating row for full results table for ", name, ": ", e$message)
      return(NULL)
    })
  }))
  
  # Filter out NULLs
  full_results <- full_results[!sapply(full_results, is.null)]
  
  # Create tables for highly significant SNPs
  high_sig_snps <- do.call(rbind, lapply(names(significant_results), function(name) {
    tryCatch({
      res <- significant_results[[name]]
      if(is.null(res) || is.null(res$highly_significant_snps)) return(NULL)
      
      # Add 'query_snp' column to make it easier to track
      cbind(
        protein = res$protein,
        protein_symbol = res$protein_symbol,
        query_snp = res$analyzed_snp,
        region = res$region,
        H4 = res$H4,
        res$highly_significant_snps
      )
    }, error = function(e) {
      message("Error creating highly significant SNPs for ", name, ": ", e$message)
      return(NULL)
    })
  }))
  
  # Filter out NULLs from high_sig_snps
  high_sig_snps <- high_sig_snps[!sapply(high_sig_snps, is.null)]
  
  # Create table for credible sets
  credible_set_snps <- do.call(rbind, lapply(names(significant_results), function(name) {
    tryCatch({
      res <- significant_results[[name]]
      if(is.null(res) || is.null(res$credible_set_snps) || nrow(res$credible_set_snps) == 0) {
        message("No credible set SNPs available for ", name)
        return(NULL)
      }
      
      # Add 'query_snp' column to make it easier to track
      cbind(
        protein = res$protein,
        protein_symbol = res$protein_symbol,
        query_snp = res$analyzed_snp,
        region = res$region,
        H4 = res$H4,
        res$credible_set_snps
      )
    }, error = function(e) {
      message("Error creating credible set SNPs for ", name, ": ", e$message)
      print(traceback())  # Show detailed stack trace
      return(NULL)
    })
  }))
  
  # Filter out NULLs from credible_set_snps
  credible_set_snps <- credible_set_snps[!sapply(credible_set_snps, is.null)]
  
  return(list(
    full_results = if(!is.null(full_results) && length(full_results) > 0 && nrow(full_results) > 0) 
                    full_results[order(-full_results$PPH4),] else data.frame(),
    highly_significant_snps = if(!is.null(high_sig_snps) && length(high_sig_snps) > 0 && nrow(high_sig_snps) > 0)
                               high_sig_snps[order(-high_sig_snps$SNP_PP_H4),] else data.frame(),
    credible_set_snps = if(!is.null(credible_set_snps) && length(credible_set_snps) > 0 && nrow(credible_set_snps) > 0) 
                         credible_set_snps[order(-credible_set_snps$SNP_PP_H4),] else data.frame()
  ))
}

# Main execution
all_results <- list()
cat("Setting up protein data for", protein_data$protein_symbol, "\n")
combined_data <- setup_protein_data(protein_data)

for(snp in protein_data$snps) {
  cat("Processing SNP:", snp, "\n")
  result <- process_single_snp(protein_data, snp, combined_data)
  if(!is.null(result)) {
    all_results[[paste(protein_data$protein, snp, sep="_")]] <- result
  }
}

# Process results and create plots
cat("Processing significant results\n")
significant_results <- process_significant_results(all_results)

# Create summary tables
cat("Creating summary tables\n")
results_tables <- create_summary_tables(all_results, significant_results)

# Create a protein-specific output directory
protein_results_dir <- file.path(config$results_dir, protein_data$protein_symbol)
if (!dir.exists(protein_results_dir)) dir.create(protein_results_dir, recursive = TRUE)

# Save Excel file with multiple sheets
excel_file <- file.path(protein_results_dir, paste0(protein_data$protein_symbol, "_coloc_results.xlsx"))
sheets_list <- list(
  "All Results" = results_tables$full_results
)

if(!is.null(results_tables$highly_significant_snps) && nrow(results_tables$highly_significant_snps) > 0) {
  sheets_list[["Highly Significant SNPs"]] <- results_tables$highly_significant_snps
}

if(!is.null(results_tables$credible_set_snps) && nrow(results_tables$credible_set_snps) > 0) {
  sheets_list[["95% Credible Set"]] <- results_tables$credible_set_snps
  cat("Added 95% Credible Set sheet with", nrow(results_tables$credible_set_snps), "SNPs\n")
} else {
  cat("WARNING: No credible set SNPs available to create sheet\n")
}

tryCatch({
  write_xlsx(sheets_list, excel_file)
  cat("Successfully created Excel file:", excel_file, "\n")
}, error = function(e) {
  cat("Error creating Excel file:", e$message, "\n")
  print(traceback())  # Show detailed stack trace
})

# Save plots for significant results - Save to plots_dir instead of protein_results_dir
for (result_name in names(significant_results)) {
  # Create clean filename
  clean_name <- gsub("[^[:alnum:]]", "_", result_name)
  # Save to plots directory instead of protein directory
  pdf_path <- file.path(plots_dir, paste0(protein_data$protein_symbol, "_", clean_name, ".pdf"))
  
  tryCatch({
    pdf(pdf_path, width=10, height=8)  # Larger size for better visibility
    print(significant_results[[result_name]]$plot)
    dev.off()
    cat("Successfully created plot:", pdf_path, "\n")
  }, error = function(e) {
    cat("Error creating plot:", pdf_path, "\n")
    cat("Error message:", e$message, "\n")
  })
}

# Print summary
cat("\nAnalysis complete for", protein_data$protein_symbol, "!\n")
cat("Number of total regions analyzed:", nrow(results_tables$full_results), "\n")
cat("Number of significant regions (H4 >= 0.8):", length(significant_results), "\n")
cat("Number of highly significant SNPs (PP.H4 > 0.8):", 
    ifelse(is.null(results_tables$highly_significant_snps) || nrow(results_tables$highly_significant_snps) == 0, 
           0, nrow(results_tables$highly_significant_snps)), "\n")
cat("Number of SNPs in 95% credible sets:", 
    ifelse(is.null(results_tables$credible_set_snps) || nrow(results_tables$credible_set_snps) == 0, 
           0, nrow(results_tables$credible_set_snps)), "\n")
cat("\nResults saved in:", protein_results_dir, "\n")
cat("Plots saved in:", plots_dir, "\n")
