#!/usr/bin/env Rscript
# Detailed comparison of current LM approach vs full DESeq2 approach

suppressMessages({
  library(data.table)
  library(dplyr)
  library(DESeq2)
  library(ggplot2)
})

cat("=== DETAILED COMPARISON: Current Pipeline vs Full DESeq2 ===\n\n")

# Read data
corrected_data <- fread("output/proteomics_enhanced_analysis/log2_intensity_corrected.csv", data.table = FALSE)
rownames(corrected_data) <- corrected_data$V1
corrected_data$V1 <- NULL

annotation <- fread("output/proteomics_enhanced_analysis/sample_annotation.csv", data.table = FALSE)

# Current pipeline results for Sex
current_results <- fread("output/proteomics_enhanced_analysis/corrected_data_results/corrected_Sex_categorical.csv", data.table = FALSE)

cat("CURRENT PIPELINE RESULTS (LM-based)\n")
cat("====================================\n")
current_sig <- sum(current_results$Adj_P_value < 0.05, na.rm = TRUE)
cat("Significant proteins:", current_sig, "out of", nrow(current_results), "\n")
cat("Mean effect size:", round(mean(abs(current_results$Effect_Size), na.rm = TRUE), 3), "\n")
cat("Method: Linear models on log2-transformed intensities\n\n")

# Full DESeq2 approach on all proteins
cat("FULL DESeq2 APPROACH\n")
cat("====================\n")

# Convert to pseudo-counts for DESeq2
pseudo_counts <- round(2^corrected_data)
pseudo_counts[is.na(pseudo_counts)] <- 0
pseudo_counts[pseudo_counts < 1] <- 1

# Prepare metadata
coldata <- data.frame(
  Sample_ID = colnames(pseudo_counts),
  Sex = factor(annotation$Gender, levels = c("Female", "Male")),
  Batch = factor(annotation$Batch),
  Age = annotation$Age,
  Location = factor(annotation$Location),
  stringsAsFactors = FALSE
)
rownames(coldata) <- coldata$Sample_ID

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = pseudo_counts,
  colData = coldata,
  design = ~ Batch + Sex
)

# Filter low count genes (DESeq2 best practice)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

cat("Proteins retained after DESeq2 filtering:", nrow(dds), "out of", nrow(pseudo_counts), "\n")

# Run DESeq2
dds <- DESeq(dds, quiet = TRUE)

# Get Sex results
res_sex <- results(dds, name = "Sex_Male_vs_Female", independentFiltering = FALSE)
deseq2_results <- as.data.frame(res_sex)

deseq2_sig <- sum(deseq2_results$padj < 0.05, na.rm = TRUE)
cat("Significant proteins:", deseq2_sig, "out of", nrow(deseq2_results), "\n")
cat("Mean effect size:", round(mean(abs(deseq2_results$log2FoldChange), na.rm = TRUE), 3), "\n")
cat("Method: Negative binomial GLM with dispersion estimation\n\n")

# Compare overlapping proteins
common_proteins <- intersect(current_results$Protein_ID, rownames(deseq2_results))
cat("DIRECT COMPARISON ON", length(common_proteins), "COMMON PROTEINS\n")
cat("=======================================================\n")

current_subset <- current_results[current_results$Protein_ID %in% common_proteins, ]
deseq2_subset <- deseq2_results[common_proteins, ]

# Match by protein ID
current_subset <- current_subset[match(common_proteins, current_subset$Protein_ID), ]
deseq2_subset <- deseq2_subset[common_proteins, ]

# Compare significance calls
current_sig_subset <- current_subset$Adj_P_value < 0.05
deseq2_sig_subset <- deseq2_subset$padj < 0.05

# Remove NAs for comparison
valid_comparison <- !is.na(current_sig_subset) & !is.na(deseq2_sig_subset)
current_sig_clean <- current_sig_subset[valid_comparison]
deseq2_sig_clean <- deseq2_sig_subset[valid_comparison]

agreement <- sum(current_sig_clean == deseq2_sig_clean) / length(current_sig_clean)
cat("Agreement in significance calls:", round(agreement * 100, 1), "%\n")

# Count different categories
both_sig <- sum(current_sig_clean & deseq2_sig_clean)
lm_only <- sum(current_sig_clean & !deseq2_sig_clean)
deseq2_only <- sum(!current_sig_clean & deseq2_sig_clean)
neither_sig <- sum(!current_sig_clean & !deseq2_sig_clean)

cat("Both significant:", both_sig, "\n")
cat("LM only significant:", lm_only, "\n")
cat("DESeq2 only significant:", deseq2_only, "\n")
cat("Neither significant:", neither_sig, "\n\n")

# Effect size correlation
valid_effects <- !is.na(current_subset$Effect_Size) & !is.na(deseq2_subset$log2FoldChange)
if (sum(valid_effects) > 10) {
  correlation <- cor(current_subset$Effect_Size[valid_effects], 
                    deseq2_subset$log2FoldChange[valid_effects])
  cat("Effect size correlation:", round(correlation, 3), "\n\n")
}

cat("PRACTICAL IMPLICATIONS\n")
cat("======================\n")
cat("1. SENSITIVITY:\n")
cat("   - DESeq2 found", deseq2_only, "additional significant proteins\n")
cat("   - LM found", lm_only, "proteins that DESeq2 missed\n\n")

cat("2. STATISTICAL POWER:\n")
cat("   - DESeq2: Uses gene-specific dispersion → better power for low-count genes\n")
cat("   - LM: Assumes constant variance → may miss true positives\n\n")

cat("3. FALSE DISCOVERY CONTROL:\n")
cat("   - DESeq2: Built-in independent filtering reduces false positives\n")
cat("   - LM: Standard BH correction only\n\n")

cat("4. BIOLOGICAL INTERPRETATION:\n")
cat("   - DESeq2: Log2 fold changes directly interpretable\n")
cat("   - LM: Effect sizes are differences in log2 space\n\n")

cat("RECOMMENDATIONS\n")
cat("===============\n")
if (deseq2_only > lm_only) {
  cat("✓ DESeq2 approach finds more significant proteins\n")
  cat("✓ Better suited for count-like proteomics data\n")
  cat("✓ More robust statistical framework\n")
} else {
  cat("• Both methods perform similarly on this dataset\n")
  cat("• Choice depends on downstream analysis needs\n")
}

cat("\nCURRENT PIPELINE STATUS:\n")
cat("- Uses DESeq2 library for batch correction framework\n")
cat("- Uses LM for annotation analysis (fast, simple)\n")
cat("- Could be enhanced to use full DESeq2 statistical pipeline\n")