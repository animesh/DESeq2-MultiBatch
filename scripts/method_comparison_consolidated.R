#!/usr/bin/env Rscript
suppressMessages({
  library(data.table)
  library(dplyr)
})

cat("=== COMPREHENSIVE METHOD COMPARISON ===\n\n")

read_deseq2_results <- function() {
  results_dir <- "output/proteomics_enhanced_analysis/corrected_data_results"
  if (!dir.exists(results_dir)) {
    stop("DESeq2 results directory not found. Run the main pipeline first.")
  }
  
  result_files <- list.files(results_dir, pattern = "corrected_.*\\.csv$", full.names = TRUE)
  deseq2_results <- list()
  
  for (file in result_files) {
    factor_name <- gsub(".*corrected_(.*)_(categorical|numerical)\\.csv", "\\1", basename(file))
    deseq2_results[[factor_name]] <- fread(file)
  }
  
  return(deseq2_results)
}

get_historical_lm_results <- function() {
  data.frame(
    Factor = c("Age", "Batch", "Days_Until_Relapse", "Location", "Relapse_Risk", "Sex", "Subject", "Treatment"),
    LM_Significant = c(79, 0, 81, 12, 330, 396, 429, 0),
    stringsAsFactors = FALSE
  )
}

summarize_results <- function(results, method_name) {
  cat(paste(method_name, "RESULTS\n"))
  cat(paste(rep("=", nchar(method_name) + 8), collapse = ""), "\n")
  
  if (method_name == "DESEQ2") {
    summary_data <- data.frame(
      Factor = names(results),
      Significant = sapply(results, function(x) sum(x$Adj_P_value < 0.05, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  } else {
    summary_data <- results
    names(summary_data)[2] <- "Significant"
  }
  
  for (i in 1:nrow(summary_data)) {
    cat(sprintf("%-20s: %4d proteins\n", summary_data$Factor[i], summary_data$Significant[i]))
  }
  cat("\n")
  
  return(summary_data)
}

compare_methods <- function(lm_results, deseq2_results) {
  comparison <- merge(lm_results, deseq2_results, by = "Factor")
  names(comparison) <- c("Factor", "LM_Significant", "DESeq2_Significant")
  
  comparison$Difference <- comparison$DESeq2_Significant - comparison$LM_Significant
  comparison$Percent_Change <- ifelse(comparison$LM_Significant > 0,
                                     round((comparison$Difference / comparison$LM_Significant) * 100, 1),
                                     ifelse(comparison$DESeq2_Significant > 0, Inf, 0))
  
  comparison[order(-comparison$Difference), ]
}

print_comparison_table <- function(comparison) {
  cat("DIRECT COMPARISON\n")
  cat("=================\n")
  cat(sprintf("%-20s %8s %8s %10s %12s\n", "Factor", "LM", "DESeq2", "Diff", "% Change"))
  cat(paste(rep("=", 65), collapse = ""), "\n")
  
  for (i in 1:nrow(comparison)) {
    pct_str <- if (is.infinite(comparison$Percent_Change[i])) "NEW"
               else if (comparison$Percent_Change[i] == 0) "0%"
               else paste0(ifelse(comparison$Percent_Change[i] > 0, "+", ""), comparison$Percent_Change[i], "%")
    
    cat(sprintf("%-20s %8d %8d %10d %12s\n",
                comparison$Factor[i], comparison$LM_Significant[i], comparison$DESeq2_Significant[i],
                comparison$Difference[i], pct_str))
  }
  cat("\n")
}

print_summary_stats <- function(comparison) {
  total_lm <- sum(comparison$LM_Significant)
  total_deseq2 <- sum(comparison$DESeq2_Significant)
  total_diff <- total_deseq2 - total_lm
  total_pct <- round((total_diff / total_lm) * 100, 1)
  
  cat("OVERALL SUMMARY\n")
  cat("===============\n")
  cat("Total LM significant proteins:     ", total_lm, "\n")
  cat("Total DESeq2 significant proteins: ", total_deseq2, "\n")
  cat("Net difference:                    ", total_diff, "proteins\n")
  cat("Overall improvement:               ", total_pct, "%\n\n")
  
  return(list(total_lm = total_lm, total_deseq2 = total_deseq2, total_diff = total_diff, total_pct = total_pct))
}

analyze_winners_losers <- function(comparison) {
  winners <- comparison[comparison$Difference > 0, ]
  losers <- comparison[comparison$Difference < 0, ]
  unchanged <- comparison[comparison$Difference == 0, ]
  
  cat("DETAILED ANALYSIS\n")
  cat("=================\n")
  
  if (nrow(winners) > 0) {
    cat("FACTORS WHERE DESEQ2 EXCELS:\n")
    winners_sorted <- winners[order(-winners$Difference), ]
    for (i in 1:nrow(winners_sorted)) {
      pct_str <- if (is.infinite(winners_sorted$Percent_Change[i])) " (NEW discoveries)"
                 else paste0(" (+", winners_sorted$Percent_Change[i], "%)")
      cat(sprintf("  %s: +%d proteins%s\n", winners_sorted$Factor[i], winners_sorted$Difference[i], pct_str))
    }
    cat("\n")
  }
  
  if (nrow(losers) > 0) {
    cat("FACTORS WHERE LM PERFORMED BETTER:\n")
    losers_sorted <- losers[order(losers$Difference), ]
    for (i in 1:nrow(losers_sorted)) {
      cat(sprintf("  %s: %d proteins (%+.1f%%)\n", losers_sorted$Factor[i], losers_sorted$Difference[i], losers_sorted$Percent_Change[i]))
    }
    cat("\n")
  }
  
  if (nrow(unchanged) > 0) {
    cat("FACTORS WITH IDENTICAL RESULTS:\n")
    for (i in 1:nrow(unchanged)) {
      cat(sprintf("  %s: %d proteins (no change)\n", unchanged$Factor[i], unchanged$LM_Significant[i]))
    }
    cat("\n")
  }
  
  return(list(winners = winners, losers = losers, unchanged = unchanged))
}

print_insights <- function(analysis_results, summary_stats) {
  winners <- analysis_results$winners
  losers <- analysis_results$losers
  
  cat("KEY INSIGHTS\n")
  cat("============\n")
  
  if (nrow(winners) > 0) {
    biggest_winner <- winners[which.max(winners$Difference), ]
    cat("• Biggest improvement:", biggest_winner$Factor, "with", biggest_winner$Difference, "additional proteins\n")
  }
  
  if (nrow(losers) > 0) {
    biggest_loser <- losers[which.min(losers$Difference), ]
    cat("• Biggest loss:", biggest_loser$Factor, "with", abs(biggest_loser$Difference), "fewer proteins\n")
  }
  
  success_factors <- nrow(winners)
  total_factors <- nrow(winners) + nrow(losers) + nrow(analysis_results$unchanged)
  success_rate <- round((success_factors / total_factors) * 100, 1)
  
  cat("• DESeq2 improved results in", success_factors, "out of", total_factors, "factors (", success_rate, "% success rate)\n")
}

print_conclusion <- function(summary_stats) {
  total_diff <- summary_stats$total_diff
  
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
}

print_methodological_notes <- function() {
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
}

main <- function() {
  deseq2_results <- read_deseq2_results()
  lm_historical <- get_historical_lm_results()
  
  deseq2_summary <- summarize_results(deseq2_results, "DESEQ2")
  lm_summary <- summarize_results(lm_historical, "LM (HISTORICAL)")
  
  comparison <- compare_methods(lm_historical, deseq2_summary)
  print_comparison_table(comparison)
  
  summary_stats <- print_summary_stats(comparison)
  analysis_results <- analyze_winners_losers(comparison)
  print_insights(analysis_results, summary_stats)
  print_conclusion(summary_stats)
  print_methodological_notes()
}

if (!interactive()) {
  main()
}