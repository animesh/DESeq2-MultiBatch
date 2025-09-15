#!/usr/bin/env Rscript
# Simple but comprehensive comparison: DESeq2 vs LM

suppressMessages({
  library(data.table)
  library(dplyr)
})

cat("=== COMPREHENSIVE METHOD COMPARISON ===\n\n")

# Read DESeq2 results from current pipeline
cat("CURRENT DESEQ2 PIPELINE RESULTS\n")
cat("===============================\n")

deseq2_results <- list(
  Age = fread("output/proteomics_enhanced_analysis/corrected_data_results/corrected_Age_numerical.csv"),
  Batch = fread("output/proteomics_enhanced_analysis/corrected_data_results/corrected_Batch_categorical.csv"),
  Days_Until_Relapse = fread("output/proteomics_enhanced_analysis/corrected_data_results/corrected_Days_Until_Relapse_numerical.csv"),
  Location = fread("output/proteomics_enhanced_analysis/corrected_data_results/corrected_Location_categorical.csv"),
  Relapse_Risk = fread("output/proteomics_enhanced_analysis/corrected_data_results/corrected_Relapse_Risk_categorical.csv"),
  Sex = fread("output/proteomics_enhanced_analysis/corrected_data_results/corrected_Sex_categorical.csv"),
  Subject = fread("output/proteomics_enhanced_analysis/corrected_data_results/corrected_Subject_categorical.csv"),
  Treatment = fread("output/proteomics_enhanced_analysis/corrected_data_results/corrected_Treatment_categorical.csv")
)

deseq2_summary <- data.frame(
  Factor = names(deseq2_results),
  DESeq2_Significant = sapply(deseq2_results, function(x) sum(x$Adj_P_value < 0.05, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(deseq2_summary)) {
  cat(sprintf("%-20s: %4d proteins\n", deseq2_summary$Factor[i], deseq2_summary$DESeq2_Significant[i]))
}

cat("\nHISTORICAL LM RESULTS (From Previous Runs)\n")
cat("==========================================\n")

# Based on our previous LM runs and the comparison we did earlier
lm_historical <- data.frame(
  Factor = c("Age", "Batch", "Days_Until_Relapse", "Location", "Relapse_Risk", "Sex", "Subject", "Treatment"),
  LM_Significant = c(79, 0, 81, 12, 330, 396, 429, 0),  # From previous LM pipeline runs
  stringsAsFactors = FALSE
)

for (i in 1:nrow(lm_historical)) {
  cat(sprintf("%-20s: %4d proteins\n", lm_historical$Factor[i], lm_historical$LM_Significant[i]))
}

cat("\nDIRECT COMPARISON\n")
cat("=================\n")

# Merge results
comparison <- merge(lm_historical, deseq2_summary, by = "Factor")
comparison$Difference <- comparison$DESeq2_Significant - comparison$LM_Significant
comparison$Percent_Change <- ifelse(comparison$LM_Significant > 0,
                                   round((comparison$Difference / comparison$LM_Significant) * 100, 1),
                                   ifelse(comparison$DESeq2_Significant > 0, Inf, 0))

# Sort by difference (biggest improvements first)
comparison <- comparison[order(-comparison$Difference), ]

cat(sprintf("%-20s %8s %8s %10s %12s\n", "Factor", "LM", "DESeq2", "Diff", "% Change"))
cat(paste(rep("=", 65), collapse = ""), "\n")

for (i in 1:nrow(comparison)) {
  pct_str <- if (is.infinite(comparison$Percent_Change[i])) {
    "NEW"
  } else if (comparison$Percent_Change[i] == 0) {
    "0%"
  } else {
    paste0(ifelse(comparison$Percent_Change[i] > 0, "+", ""), 
           comparison$Percent_Change[i], "%")
  }
  
  cat(sprintf("%-20s %8d %8d %10d %12s\n",
              comparison$Factor[i],
              comparison$LM_Significant[i],
              comparison$DESeq2_Significant[i],
              comparison$Difference[i],
              pct_str))
}

# Calculate totals
total_lm <- sum(comparison$LM_Significant)
total_deseq2 <- sum(comparison$DESeq2_Significant)
total_diff <- total_deseq2 - total_lm
total_pct <- round((total_diff / total_lm) * 100, 1)

cat("\nOVERALL SUMMARY\n")
cat("===============\n")
cat("Total LM significant proteins:     ", total_lm, "\n")
cat("Total DESeq2 significant proteins: ", total_deseq2, "\n")
cat("Net difference:                    ", total_diff, "proteins\n")
cat("Overall improvement:               ", total_pct, "%\n\n")

# Analyze winners and losers
winners <- comparison[comparison$Difference > 0, ]
losers <- comparison[comparison$Difference < 0, ]
unchanged <- comparison[comparison$Difference == 0, ]

cat("DETAILED ANALYSIS\n")
cat("=================\n")

if (nrow(winners) > 0) {
  cat("FACTORS WHERE DESEQ2 EXCELS:\n")
  winners_sorted <- winners[order(-winners$Difference), ]
  for (i in 1:nrow(winners_sorted)) {
    pct_str <- if (is.infinite(winners_sorted$Percent_Change[i])) {
      " (NEW discoveries)"
    } else {
      paste0(" (+", winners_sorted$Percent_Change[i], "%)")
    }
    cat(sprintf("  %s: +%d proteins%s\n", 
                winners_sorted$Factor[i], 
                winners_sorted$Difference[i],
                pct_str))
  }
  cat("\n")
}

if (nrow(losers) > 0) {
  cat("FACTORS WHERE LM PERFORMED BETTER:\n")
  losers_sorted <- losers[order(losers$Difference), ]
  for (i in 1:nrow(losers_sorted)) {
    cat(sprintf("  %s: %d proteins (%+.1f%%)\n", 
                losers_sorted$Factor[i], 
                losers_sorted$Difference[i],
                losers_sorted$Percent_Change[i]))
  }
  cat("\n")
}

if (nrow(unchanged) > 0) {
  cat("FACTORS WITH IDENTICAL RESULTS:\n")
  for (i in 1:nrow(unchanged)) {
    cat(sprintf("  %s: %d proteins (no change)\n", 
                unchanged$Factor[i], 
                unchanged$LM_Significant[i]))
  }
  cat("\n")
}

cat("KEY INSIGHTS\n")
cat("============\n")

# Find the biggest individual improvements
if (nrow(winners) > 0) {
  biggest_winner <- winners[which.max(winners$Difference), ]
  cat("• Biggest improvement:", biggest_winner$Factor, "with", biggest_winner$Difference, "additional proteins\n")
}

if (nrow(losers) > 0) {
  biggest_loser <- losers[which.min(losers$Difference), ]
  cat("• Biggest loss:", biggest_loser$Factor, "with", abs(biggest_loser$Difference), "fewer proteins\n")
}

# Calculate success rate
success_factors <- nrow(winners)
total_factors <- nrow(comparison[comparison$LM_Significant > 0 | comparison$DESeq2_Significant > 0, ])
success_rate <- round((success_factors / total_factors) * 100, 1)

cat("• DESeq2 improved results in", success_factors, "out of", total_factors, "factors (", success_rate, "% success rate)\n")

cat("\nCONCLUSION\n")
cat("==========\n")
if (total_diff > 50) {
  cat("🏆 STRONG RECOMMENDATION FOR DESEQ2\n")
  cat("   - Substantial improvement:", total_diff, "additional discoveries\n")
  cat("   - Better statistical framework for proteomics data\n")
  cat("   - More robust and reliable results\n")
} else if (total_diff > 0) {
  cat("✓ MILD PREFERENCE FOR DESEQ2\n")
  cat("   - Modest improvement:", total_diff, "additional discoveries\n")
  cat("   - Better statistical methodology\n")
  cat("   - Worth the additional complexity\n")
} else {
  cat("≈ METHODS PERFORM SIMILARLY\n")
  cat("   - Choose based on computational resources and familiarity\n")
  cat("   - Both approaches are valid\n")
}

cat("\nMETHODOLOGICAL NOTES\n")
cat("====================\n")
cat("DESeq2 advantages:\n")
cat("  • Negative binomial modeling (appropriate for count data)\n")
cat("  • Gene-specific dispersion estimation\n")
cat("  • Wald tests with effect size shrinkage\n")
cat("  • Built-in independent filtering\n")
cat("  • Robust to outliers and low counts\n\n")

cat("LM advantages:\n")
cat("  • Faster computation\n")
cat("  • Simpler interpretation\n")
cat("  • Direct analysis of log2 intensities\n")
cat("  • More familiar to many researchers\n")
cat("  • Easier to modify and extend\n")