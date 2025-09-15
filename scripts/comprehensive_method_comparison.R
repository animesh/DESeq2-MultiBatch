#!/usr/bin/env Rscript
# Comprehensive comparison: Full DESeq2 vs LM-only on all factors

suppressMessages({
  library(data.table)
  library(dplyr)
})

cat("=== COMPREHENSIVE COMPARISON: Full DESeq2 vs LM-Only ===\n\n")

# Read all DESeq2 results
deseq2_files <- list.files("output/proteomics_enhanced_analysis/corrected_data_results/", 
                          pattern = "corrected_.*\\.csv", full.names = TRUE)

deseq2_results <- list()
for (file in deseq2_files) {
  factor_name <- gsub(".*corrected_(.*)_(categorical|numerical)\\.csv", "\\1", basename(file))
  deseq2_results[[factor_name]] <- fread(file, data.table = FALSE)
}

cat("DESEQ2 RESULTS (Current Pipeline)\n")
cat("==================================\n")
deseq2_summary <- data.frame(
  Factor = names(deseq2_results),
  DESeq2_Significant = sapply(deseq2_results, function(x) sum(x$Adj_P_value < 0.05, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(deseq2_summary)) {
  cat(sprintf("%-20s: %4d significant proteins\n", 
              deseq2_summary$Factor[i], deseq2_summary$DESeq2_Significant[i]))
}

# Load data for LM analysis
corrected_data <- fread("output/proteomics_enhanced_analysis/log2_intensity_corrected.csv", data.table = FALSE)
rownames(corrected_data) <- corrected_data$V1
corrected_data$V1 <- NULL

annotation <- fread("output/proteomics_enhanced_analysis/sample_annotation.csv", data.table = FALSE)

cat("\nLM-ONLY ANALYSIS (Same Corrected Data)\n")
cat("======================================\n")

# Generic LM analysis function
run_lm_analysis <- function(intensity_matrix, annotation, variable_col, is_numerical = FALSE) {
  var_values <- annotation[[variable_col]]
  valid_samples <- !is.na(var_values)
  
  if (sum(valid_samples) < 10) return(0)
  
  intensity_clean <- intensity_matrix[, valid_samples]
  var_clean <- var_values[valid_samples]
  
  # Create design matrix
  if (is_numerical) {
    design_data <- data.frame(
      Variable = as.numeric(var_clean),
      Batch = factor(annotation$Batch[valid_samples])
    )
  } else {
    design_data <- data.frame(
      Variable = factor(var_clean),
      Batch = factor(annotation$Batch[valid_samples])
    )
    if (length(levels(design_data$Variable)) < 2) return(0)
  }
  
  design_matrix <- model.matrix(~ Variable + Batch, data = design_data)
  
  significant_count <- 0
  p_values <- numeric(nrow(intensity_clean))
  
  for (i in 1:nrow(intensity_clean)) {
    protein_data <- as.numeric(intensity_clean[i, ])
    valid_protein <- !is.na(protein_data)
    if (sum(valid_protein) < 10) next
    
    y <- protein_data[valid_protein]
    X <- design_matrix[valid_protein, ]
    
    tryCatch({
      fit <- lm(y ~ X - 1)
      var_coef_idx <- grep("Variable", colnames(X))
      if (length(var_coef_idx) > 0) {
        p_values[i] <- summary(fit)$coefficients[var_coef_idx[1], 4]
      }
    }, error = function(e) {})
  }
  
  # Adjust p-values
  valid_pvals <- !is.na(p_values) & p_values > 0
  if (sum(valid_pvals) > 0) {
    adj_p_values <- rep(1, length(p_values))
    adj_p_values[valid_pvals] <- p.adjust(p_values[valid_pvals], method = "BH")
    significant_count <- sum(adj_p_values < 0.05, na.rm = TRUE)
  }
  
  return(significant_count)
}

# Run LM analysis for each factor
lm_results <- list(
  Age = run_lm_analysis(corrected_data, annotation, "Age", TRUE),
  Batch = run_lm_analysis(corrected_data, annotation, "Batch", FALSE),
  Days_Until_Relapse = run_lm_analysis(corrected_data, annotation, "Days until relaps", TRUE),
  Location = run_lm_analysis(corrected_data, annotation, "Location", FALSE),
  Relapse_Risk = run_lm_analysis(corrected_data, annotation, "Relapse_Risk", FALSE),
  Sex = run_lm_analysis(corrected_data, annotation, "Gender", FALSE),
  Subject = run_lm_analysis(corrected_data, annotation, "Subject", FALSE),
  Treatment = run_lm_analysis(corrected_data, annotation, "Treatment", FALSE)
)

for (factor in names(lm_results)) {
  cat(sprintf("%-20s: %4d significant proteins\n", factor, lm_results[[factor]]))
}

cat("\nDIRECT COMPARISON\n")
cat("=================\n")

# Create comparison table
comparison_data <- data.frame(
  Factor = names(lm_results),
  LM_Significant = unlist(lm_results),
  stringsAsFactors = FALSE
)

# Match with DESeq2 results
comparison_data$DESeq2_Significant <- sapply(comparison_data$Factor, function(f) {
  idx <- which(deseq2_summary$Factor == f)
  if (length(idx) > 0) deseq2_summary$DESeq2_Significant[idx] else 0
})

comparison_data$Difference <- comparison_data$DESeq2_Significant - comparison_data$LM_Significant
comparison_data$Percent_Change <- ifelse(comparison_data$LM_Significant > 0,
                                        round((comparison_data$Difference / comparison_data$LM_Significant) * 100, 1),
                                        ifelse(comparison_data$DESeq2_Significant > 0, Inf, 0))

# Sort by difference
comparison_data <- comparison_data[order(-comparison_data$Difference), ]

cat(sprintf("%-20s %8s %8s %10s %12s\n", "Factor", "LM", "DESeq2", "Diff", "% Change"))
cat(paste(rep("=", 65), collapse = ""), "\n")

for (i in 1:nrow(comparison_data)) {
  pct_str <- if (is.infinite(comparison_data$Percent_Change[i])) {
    "NEW"
  } else if (comparison_data$Percent_Change[i] == 0) {
    "0%"
  } else {
    paste0(ifelse(comparison_data$Percent_Change[i] > 0, "+", ""), 
           comparison_data$Percent_Change[i], "%")
  }
  
  cat(sprintf("%-20s %8d %8d %10d %12s\n",
              comparison_data$Factor[i],
              comparison_data$LM_Significant[i],
              comparison_data$DESeq2_Significant[i],
              comparison_data$Difference[i],
              pct_str))
}

# Summary statistics
total_lm <- sum(comparison_data$LM_Significant)
total_deseq2 <- sum(comparison_data$DESeq2_Significant)
total_diff <- total_deseq2 - total_lm
total_pct <- round((total_diff / total_lm) * 100, 1)

cat("\nSUMMARY\n")
cat("=======\n")
cat("Total LM significant proteins:", total_lm, "\n")
cat("Total DESeq2 significant proteins:", total_deseq2, "\n")
cat("Overall difference:", total_diff, "proteins\n")
cat("Overall improvement:", total_pct, "%\n\n")

# Identify biggest winners and losers
winners <- comparison_data[comparison_data$Difference > 0, ]
losers <- comparison_data[comparison_data$Difference < 0, ]

if (nrow(winners) > 0) {
  cat("BIGGEST IMPROVEMENTS (DESeq2 > LM):\n")
  for (i in 1:min(3, nrow(winners))) {
    cat(sprintf("- %s: +%d proteins (%s)\n", 
                winners$Factor[i], 
                winners$Difference[i],
                ifelse(is.infinite(winners$Percent_Change[i]), "NEW", 
                       paste0("+", winners$Percent_Change[i], "%"))))
  }
  cat("\n")
}

if (nrow(losers) > 0) {
  cat("AREAS WHERE LM PERFORMED BETTER:\n")
  for (i in 1:min(3, nrow(losers))) {
    cat(sprintf("- %s: %d proteins (%+.1f%%)\n", 
                losers$Factor[i], 
                losers$Difference[i],
                losers$Percent_Change[i]))
  }
  cat("\n")
}

cat("INTERPRETATION\n")
cat("==============\n")
if (total_diff > 0) {
  cat("✓ DESeq2 provides", total_diff, "additional discoveries (", total_pct, "% improvement)\n")
  cat("✓ Better statistical power for most biological factors\n")
  cat("✓ More appropriate statistical model for count-like proteomics data\n")
} else {
  cat("• Methods perform similarly on this dataset\n")
  cat("• Choice depends on computational resources and interpretation needs\n")
}

cat("\nRECOMMENDation:\n")
if (total_pct > 5) {
  cat("Strong recommendation for DESeq2 due to significant improvement in discovery power\n")
} else if (total_pct > 0) {
  cat("Mild preference for DESeq2 due to better statistical framework\n")
} else {
  cat("Either method acceptable, choose based on computational constraints\n")
}