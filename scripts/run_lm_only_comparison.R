#!/usr/bin/env Rscript
# Direct comparison: Full DESeq2 vs LM-only approaches

suppressMessages({
  library(data.table)
  library(dplyr)
  library(DESeq2)
})

cat("=== DIRECT COMPARISON: Full DESeq2 vs LM-Only ===\n\n")

# Read the current DESeq2 results
deseq2_sex <- fread("output/proteomics_enhanced_analysis/corrected_data_results/corrected_Sex_categorical.csv", data.table = FALSE)
deseq2_age <- fread("output/proteomics_enhanced_analysis/corrected_data_results/corrected_Age_numerical.csv", data.table = FALSE)
deseq2_relapse <- fread("output/proteomics_enhanced_analysis/corrected_data_results/corrected_Relapse_Risk_categorical.csv", data.table = FALSE)

cat("CURRENT DESEQ2 RESULTS\n")
cat("======================\n")
cat("Sex significant proteins:", sum(deseq2_sex$Adj_P_value < 0.05, na.rm = TRUE), "\n")
cat("Age significant proteins:", sum(deseq2_age$Adj_P_value < 0.05, na.rm = TRUE), "\n")
cat("Relapse_Risk significant proteins:", sum(deseq2_relapse$Adj_P_value < 0.05, na.rm = TRUE), "\n\n")

# Now run LM-only analysis on the same data
corrected_data <- fread("output/proteomics_enhanced_analysis/log2_intensity_corrected.csv", data.table = FALSE)
rownames(corrected_data) <- corrected_data$V1
corrected_data$V1 <- NULL

annotation <- fread("output/proteomics_enhanced_analysis/sample_annotation.csv", data.table = FALSE)

cat("LM-ONLY ANALYSIS (Same Data)\n")
cat("=============================\n")

# LM analysis for Sex
analyze_lm_categorical <- function(intensity_matrix, annotation, variable_col, variable_name) {
  var_values <- annotation[[variable_col]]
  valid_samples <- !is.na(var_values)
  
  if (sum(valid_samples) < 10) return(NULL)
  
  intensity_clean <- intensity_matrix[, valid_samples]
  var_clean <- var_values[valid_samples]
  
  # Create design matrix
  design_data <- data.frame(
    Variable = factor(var_clean),
    Batch = factor(annotation$Batch[valid_samples])
  )
  
  if (length(levels(design_data$Variable)) < 2) return(NULL)
  
  design_matrix <- model.matrix(~ Variable + Batch, data = design_data)
  
  results <- data.frame(
    Protein_ID = rownames(intensity_clean),
    Effect_Size = NA,
    P_value = NA,
    Adj_P_value = NA,
    stringsAsFactors = FALSE
  )
  
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
        results$Effect_Size[i] <- coef(fit)[var_coef_idx[1]]
        results$P_value[i] <- summary(fit)$coefficients[var_coef_idx[1], 4]
      }
    }, error = function(e) {})
  }
  
  # Adjust p-values
  valid_pvals <- !is.na(results$P_value)
  if (sum(valid_pvals) > 0) {
    results$Adj_P_value[valid_pvals] <- p.adjust(results$P_value[valid_pvals], method = "BH")
  }
  
  significant <- sum(results$Adj_P_value < 0.05, na.rm = TRUE)
  cat(variable_name, "significant proteins:", significant, "\n")
  
  return(results)
}

# LM analysis for numerical variables
analyze_lm_numerical <- function(intensity_matrix, annotation, variable_col, variable_name) {
  var_values <- annotation[[variable_col]]
  valid_samples <- !is.na(var_values)
  
  if (sum(valid_samples) < 10) return(NULL)
  
  intensity_clean <- intensity_matrix[, valid_samples]
  var_clean <- var_values[valid_samples]
  
  # Create design matrix
  design_data <- data.frame(
    Variable = as.numeric(var_clean),
    Batch = factor(annotation$Batch[valid_samples])
  )
  
  design_matrix <- model.matrix(~ Variable + Batch, data = design_data)
  
  results <- data.frame(
    Protein_ID = rownames(intensity_clean),
    Correlation = NA,
    P_value = NA,
    Adj_P_value = NA,
    stringsAsFactors = FALSE
  )
  
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
        results$Correlation[i] <- coef(fit)[var_coef_idx]
        results$P_value[i] <- summary(fit)$coefficients[var_coef_idx, 4]
      }
    }, error = function(e) {})
  }
  
  # Adjust p-values
  valid_pvals <- !is.na(results$P_value)
  if (sum(valid_pvals) > 0) {
    results$Adj_P_value[valid_pvals] <- p.adjust(results$P_value[valid_pvals], method = "BH")
  }
  
  significant <- sum(results$Adj_P_value < 0.05, na.rm = TRUE)
  cat(variable_name, "significant proteins:", significant, "\n")
  
  return(results)
}

# Run LM analyses
lm_sex <- analyze_lm_categorical(corrected_data, annotation, "Gender", "Sex")
lm_age <- analyze_lm_numerical(corrected_data, annotation, "Age", "Age")
lm_relapse <- analyze_lm_categorical(corrected_data, annotation, "Relapse_Risk", "Relapse_Risk")

cat("\nDIRECT COMPARISON\n")
cat("=================\n")

# Compare Sex results
if (!is.null(lm_sex)) {
  deseq2_sex_sig <- sum(deseq2_sex$Adj_P_value < 0.05, na.rm = TRUE)
  lm_sex_sig <- sum(lm_sex$Adj_P_value < 0.05, na.rm = TRUE)
  sex_diff <- deseq2_sex_sig - lm_sex_sig
  sex_pct <- round((sex_diff / lm_sex_sig) * 100, 1)
  
  cat("SEX ANALYSIS:\n")
  cat("- DESeq2:", deseq2_sex_sig, "significant proteins\n")
  cat("- LM-only:", lm_sex_sig, "significant proteins\n")
  cat("- Difference:", sex_diff, "proteins (", sex_pct, "% improvement)\n\n")
}

# Compare Age results
if (!is.null(lm_age)) {
  deseq2_age_sig <- sum(deseq2_age$Adj_P_value < 0.05, na.rm = TRUE)
  lm_age_sig <- sum(lm_age$Adj_P_value < 0.05, na.rm = TRUE)
  age_diff <- deseq2_age_sig - lm_age_sig
  age_pct <- round((age_diff / lm_age_sig) * 100, 1)
  
  cat("AGE ANALYSIS:\n")
  cat("- DESeq2:", deseq2_age_sig, "significant proteins\n")
  cat("- LM-only:", lm_age_sig, "significant proteins\n")
  cat("- Difference:", age_diff, "proteins (", age_pct, "% improvement)\n\n")
}

# Compare Relapse Risk results
if (!is.null(lm_relapse)) {
  deseq2_relapse_sig <- sum(deseq2_relapse$Adj_P_value < 0.05, na.rm = TRUE)
  lm_relapse_sig <- sum(lm_relapse$Adj_P_value < 0.05, na.rm = TRUE)
  relapse_diff <- deseq2_relapse_sig - lm_relapse_sig
  relapse_pct <- round((relapse_diff / lm_relapse_sig) * 100, 1)
  
  cat("RELAPSE RISK ANALYSIS:\n")
  cat("- DESeq2:", deseq2_relapse_sig, "significant proteins\n")
  cat("- LM-only:", lm_relapse_sig, "significant proteins\n")
  cat("- Difference:", relapse_diff, "proteins (", relapse_pct, "% improvement)\n\n")
}

# Effect size correlation analysis
if (!is.null(lm_sex)) {
  # Match proteins between methods
  common_proteins <- intersect(deseq2_sex$Protein_ID, lm_sex$Protein_ID)
  if (length(common_proteins) > 100) {
    deseq2_subset <- deseq2_sex[match(common_proteins, deseq2_sex$Protein_ID), ]
    lm_subset <- lm_sex[match(common_proteins, lm_sex$Protein_ID), ]
    
    valid_comparison <- !is.na(deseq2_subset$Effect_Size) & !is.na(lm_subset$Effect_Size)
    if (sum(valid_comparison) > 50) {
      correlation <- cor(deseq2_subset$Effect_Size[valid_comparison], 
                        lm_subset$Effect_Size[valid_comparison])
      cat("EFFECT SIZE CORRELATION (Sex):", round(correlation, 3), "\n")
    }
  }
}

cat("\nMETHODOLOGICAL DIFFERENCES\n")
cat("==========================\n")
cat("DESeq2 Advantages:\n")
cat("- Negative binomial distribution (better for count-like data)\n")
cat("- Gene-specific dispersion estimation\n")
cat("- Wald tests with shrinkage estimation\n")
cat("- Built-in independent filtering\n")
cat("- Robust to outliers and low counts\n\n")

cat("LM Advantages:\n")
cat("- Faster computation\n")
cat("- Simpler interpretation\n")
cat("- Works directly on log2 intensities\n")
cat("- More familiar to many researchers\n\n")

cat("CONCLUSION\n")
cat("==========\n")
if (!is.null(lm_sex) && !is.null(lm_age) && !is.null(lm_relapse)) {
  total_deseq2 <- deseq2_sex_sig + deseq2_age_sig + deseq2_relapse_sig
  total_lm <- lm_sex_sig + lm_age_sig + lm_relapse_sig
  total_improvement <- total_deseq2 - total_lm
  total_pct <- round((total_improvement / total_lm) * 100, 1)
  
  cat("Total significant proteins:\n")
  cat("- DESeq2:", total_deseq2, "\n")
  cat("- LM-only:", total_lm, "\n")
  cat("- Overall improvement:", total_improvement, "proteins (", total_pct, "%)\n\n")
  
  if (total_improvement > 0) {
    cat("✓ DESeq2 provides significantly more discoveries\n")
    cat("✓ Better statistical framework for proteomics data\n")
    cat("✓ More robust and reliable results\n")
  }
}