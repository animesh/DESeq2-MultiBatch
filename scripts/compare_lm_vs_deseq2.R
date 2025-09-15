#!/usr/bin/env Rscript
# Comparison of Linear Model vs DESeq2 approaches

suppressMessages({
  library(data.table)
  library(dplyr)
  library(DESeq2)
})

# Read the corrected data from our pipeline
corrected_data <- fread("output/proteomics_enhanced_analysis/log2_intensity_corrected.csv", data.table = FALSE)
rownames(corrected_data) <- corrected_data$V1
corrected_data$V1 <- NULL

annotation <- fread("output/proteomics_enhanced_analysis/sample_annotation.csv", data.table = FALSE)
rownames(annotation) <- annotation$Sample_ID

# Focus on Sex analysis as an example
cat("=== COMPARISON: Linear Model vs DESeq2 for Sex Analysis ===\n\n")

# Method 1: Current Linear Model Approach
cat("METHOD 1: Linear Model (lm) Approach\n")
cat("=====================================\n")

# Prepare design matrix for lm
design_data <- annotation
design_data$Sex <- factor(annotation$Gender)  # Use Gender column
design_data$Batch <- factor(design_data$Batch)

# Create design matrix
design_matrix <- model.matrix(~ Sex + Batch, data = design_data)

lm_results <- data.frame(
  Protein_ID = rownames(corrected_data),
  LM_Effect_Size = NA,
  LM_P_value = NA,
  stringsAsFactors = FALSE
)

# Analyze first 100 proteins with lm
for (i in 1:min(100, nrow(corrected_data))) {
  protein_data <- as.numeric(corrected_data[i, ])
  valid_idx <- !is.na(protein_data)
  
  if (sum(valid_idx) < 10) next
  
  y <- protein_data[valid_idx]
  X <- design_matrix[valid_idx, ]
  
  tryCatch({
    fit <- lm(y ~ X - 1)
    sex_coef_idx <- grep("SexMale", colnames(X))
    if (length(sex_coef_idx) > 0) {
      lm_results$LM_Effect_Size[i] <- coef(fit)[sex_coef_idx[1]]
      lm_results$LM_P_value[i] <- summary(fit)$coefficients[sex_coef_idx[1], 4]
    }
  }, error = function(e) {})
}

lm_significant <- sum(p.adjust(lm_results$LM_P_value, method = "BH") < 0.05, na.rm = TRUE)
cat("LM significant proteins (first 100):", lm_significant, "\n")
cat("LM mean effect size:", round(mean(abs(lm_results$LM_Effect_Size), na.rm = TRUE), 3), "\n\n")

# Method 2: DESeq2 Approach
cat("METHOD 2: DESeq2 Approach\n")
cat("==========================\n")

# Convert log2 intensities back to approximate counts for DESeq2
pseudo_counts <- round(2^corrected_data)
pseudo_counts[is.na(pseudo_counts)] <- 0
pseudo_counts[pseudo_counts < 1] <- 1

# Take first 100 proteins for comparison
pseudo_counts_subset <- pseudo_counts[1:100, ]

# Prepare colData
coldata <- annotation
coldata$Sex <- factor(annotation$Gender, levels = c("Female", "Male"))  # Use Gender column
coldata$Batch <- factor(coldata$Batch)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = pseudo_counts_subset,
  colData = coldata,
  design = ~ Batch + Sex
)

# Run DESeq2
dds <- DESeq(dds, quiet = TRUE)

# Get results for Sex
res_sex <- results(dds, name = "Sex_Male_vs_Female", independentFiltering = FALSE)

deseq2_results <- data.frame(
  Protein_ID = rownames(res_sex),
  DESeq2_Effect_Size = res_sex$log2FoldChange,
  DESeq2_P_value = res_sex$pvalue,
  DESeq2_Adj_P_value = res_sex$padj,
  stringsAsFactors = FALSE
)

deseq2_significant <- sum(deseq2_results$DESeq2_Adj_P_value < 0.05, na.rm = TRUE)
cat("DESeq2 significant proteins (first 100):", deseq2_significant, "\n")
cat("DESeq2 mean effect size:", round(mean(abs(deseq2_results$DESeq2_Effect_Size), na.rm = TRUE), 3), "\n\n")

# Compare the two approaches
cat("COMPARISON SUMMARY\n")
cat("==================\n")
cat("Linear Model approach:\n")
cat("- Uses normal distribution assumption\n")
cat("- Simple least squares fitting\n")
cat("- Effect size: difference in means\n")
cat("- P-values from t-distribution\n\n")

cat("DESeq2 approach:\n")
cat("- Uses negative binomial distribution\n")
cat("- Accounts for count data overdispersion\n")
cat("- Effect size: log2 fold change\n")
cat("- P-values from Wald test or likelihood ratio test\n")
cat("- Built-in normalization and dispersion estimation\n\n")

# Merge results for direct comparison
comparison <- merge(lm_results[1:100, ], deseq2_results, by = "Protein_ID", all = TRUE)
comparison <- comparison[!is.na(comparison$LM_Effect_Size) & !is.na(comparison$DESeq2_Effect_Size), ]

if (nrow(comparison) > 0) {
  correlation <- cor(comparison$LM_Effect_Size, comparison$DESeq2_Effect_Size, use = "complete.obs")
  cat("Correlation between LM and DESeq2 effect sizes:", round(correlation, 3), "\n")
  
  # Count agreements in significance
  lm_sig <- p.adjust(comparison$LM_P_value, method = "BH") < 0.05
  deseq2_sig <- comparison$DESeq2_Adj_P_value < 0.05
  agreement <- sum(lm_sig == deseq2_sig, na.rm = TRUE) / sum(!is.na(lm_sig) & !is.na(deseq2_sig))
  cat("Agreement in significance calls:", round(agreement * 100, 1), "%\n")
}

cat("\nKEY DIFFERENCES:\n")
cat("1. Statistical Model: LM assumes normal distribution, DESeq2 uses negative binomial\n")
cat("2. Data Type: LM works on log2 intensities, DESeq2 designed for count data\n")
cat("3. Normalization: LM uses pre-normalized data, DESeq2 has built-in normalization\n")
cat("4. Dispersion: LM assumes constant variance, DESeq2 estimates gene-specific dispersion\n")
cat("5. Multiple Testing: Both use Benjamini-Hochberg, but DESeq2 has additional filtering\n")