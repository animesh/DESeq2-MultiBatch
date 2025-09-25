#!/usr/bin/env Rscript
suppressMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
})

RESULTS_DIR <- "output/proteomics_enhanced_analysis"
OUTPUT_DIR <- "output/overlap_analysis_results"

create_output_dir <- function() {
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  cat("Overlap analysis output directory:", OUTPUT_DIR, "\n")
}

print_section <- function(title) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat(paste0("  ", title), "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
}

dark_theme <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "#1e1e1e", color = NA),
      panel.background = element_rect(fill = "#1e1e1e", color = NA),
      panel.grid.major = element_line(color = "#404040", linewidth = 0.3),
      panel.grid.minor = element_line(color = "#303030", linewidth = 0.2),
      text = element_text(color = "#ffffff"),
      axis.text = element_text(color = "#cccccc"),
      axis.title = element_text(color = "#ffffff"),
      plot.title = element_text(color = "#ffffff", hjust = 0.5),
      legend.background = element_rect(fill = "#1e1e1e", color = NA),
      legend.text = element_text(color = "#cccccc"),
      legend.title = element_text(color = "#ffffff"),
      strip.background = element_rect(fill = "#2d2d2d", color = NA),
      strip.text = element_text(color = "#ffffff")
    )
}

load_differential_proteins <- function(results_dir, data_type = "corrected") {
  data_dir <- file.path(results_dir, paste0(data_type, "_data_results"))
  if (!dir.exists(data_dir)) stop("Results directory not found: ", data_dir)
  
  result_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
  differential_lists <- list()
  
  for (file in result_files) {
    basename_file <- basename(file)
    factor_name <- gsub(paste0(data_type, "_"), "", basename_file)
    factor_name <- gsub("_categorical\\.csv|_numerical\\.csv", "", factor_name)
    
    cat("Processing", factor_name, "...\n")
    results <- fread(file)
    
    significant_proteins <- if ("Adj_P_value" %in% colnames(results)) {
      results[Adj_P_value < 0.05 & !is.na(Adj_P_value), Protein_ID]
    } else character(0)
    
    differential_lists[[factor_name]] <- significant_proteins
    cat("  -", factor_name, ":", length(significant_proteins), "significant proteins\n")
  }
  
  return(differential_lists)
}

calculate_pairwise_overlaps <- function(diff_lists) {
  factor_names <- names(diff_lists)
  n_factors <- length(factor_names)
  
  overlap_matrix <- matrix(0, nrow = n_factors, ncol = n_factors)
  rownames(overlap_matrix) <- colnames(overlap_matrix) <- factor_names
  
  overlap_details <- data.frame(
    Factor1 = character(), Factor2 = character(), Factor1_Count = integer(),
    Factor2_Count = integer(), Overlap_Count = integer(), Jaccard_Index = numeric(),
    Overlap_Percentage_F1 = numeric(), Overlap_Percentage_F2 = numeric(),
    stringsAsFactors = FALSE)
  
  for (i in 1:n_factors) {
    for (j in 1:n_factors) {
      if (i <= j) {
        list1 <- diff_lists[[factor_names[i]]]
        list2 <- diff_lists[[factor_names[j]]]
        
        if (i == j) {
          overlap_count <- length(list1)
          jaccard <- perc1 <- perc2 <- if (length(list1) > 0) 100.0 else 0.0
        } else {
          overlap_count <- length(intersect(list1, list2))
          union_count <- length(union(list1, list2))
          jaccard <- if (union_count > 0) overlap_count / union_count else 0
          perc1 <- if (length(list1) > 0) 100 * overlap_count / length(list1) else 0
          perc2 <- if (length(list2) > 0) 100 * overlap_count / length(list2) else 0
        }
        
        overlap_matrix[i, j] <- overlap_matrix[j, i] <- overlap_count
        
        if (i != j) {
          overlap_details <- rbind(overlap_details, data.frame(
            Factor1 = factor_names[i], Factor2 = factor_names[j],
            Factor1_Count = length(list1), Factor2_Count = length(list2),
            Overlap_Count = overlap_count, Jaccard_Index = jaccard,
            Overlap_Percentage_F1 = perc1, Overlap_Percentage_F2 = perc2,
            stringsAsFactors = FALSE))
        }
      }
    }
  }
  
  list(matrix = overlap_matrix, details = overlap_details)
}

create_overlap_heatmap <- function(overlap_matrix, title) {
  heatmap_data <- expand.grid(Factor1 = rownames(overlap_matrix), Factor2 = colnames(overlap_matrix))
  heatmap_data$Overlap <- as.vector(overlap_matrix)
  
  ggplot(heatmap_data, aes(x = Factor1, y = Factor2, fill = Overlap)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = Overlap), color = "white", fontface = "bold", size = 3) +
    scale_fill_gradient(low = "navy", high = "red", name = "Overlap\nCount") +
    labs(title = title, x = "Factor", y = "Factor") +
    dark_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(angle = 0), panel.grid = element_blank())
}

create_jaccard_heatmap <- function(overlap_details, title) {
  jaccard_matrix <- matrix(1, nrow = length(unique(c(overlap_details$Factor1, overlap_details$Factor2))),
                          ncol = length(unique(c(overlap_details$Factor1, overlap_details$Factor2))))
  factor_names <- sort(unique(c(overlap_details$Factor1, overlap_details$Factor2)))
  rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- factor_names
  
  for (i in 1:nrow(overlap_details)) {
    f1 <- overlap_details$Factor1[i]
    f2 <- overlap_details$Factor2[i]
    jaccard <- overlap_details$Jaccard_Index[i]
    jaccard_matrix[f1, f2] <- jaccard_matrix[f2, f1] <- jaccard
  }
  
  heatmap_data <- expand.grid(Factor1 = factor_names, Factor2 = factor_names)
  heatmap_data$Jaccard <- as.vector(jaccard_matrix)
  
  ggplot(heatmap_data, aes(x = Factor1, y = Factor2, fill = Jaccard)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", Jaccard)), color = "white", fontface = "bold", size = 3) +
    scale_fill_gradient(low = "navy", high = "orange", name = "Jaccard\nIndex") +
    labs(title = title, x = "Factor", y = "Factor") +
    dark_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(angle = 0), panel.grid = element_blank())
}

create_factor_sizes_plot <- function(diff_lists) {
  factor_data <- data.frame(
    Factor = names(diff_lists),
    Count = sapply(diff_lists, length),
    stringsAsFactors = FALSE)
  
  factor_data <- factor_data[order(-factor_data$Count), ]
  factor_data$Factor <- factor(factor_data$Factor, levels = factor_data$Factor)
  
  ggplot(factor_data, aes(x = Factor, y = Count)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = Count), vjust = -0.5, color = "white", fontface = "bold") +
    labs(title = "Significant Proteins per Factor", x = "Factor", y = "Number of Proteins") +
    dark_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

create_overlap_summary_plot <- function(overlap_details) {
  ggplot(overlap_details, aes(x = Overlap_Count, y = Jaccard_Index)) +
    geom_point(aes(size = Overlap_Count), alpha = 0.7, color = "lightblue") +
    geom_text(aes(label = paste(Factor1, "∩", Factor2)), vjust = -0.8, size = 2.5, color = "white") +
    scale_size_continuous(range = c(2, 8), name = "Overlap\nCount") +
    labs(title = "Overlap Count vs Jaccard Similarity", x = "Overlap Count", y = "Jaccard Index") +
    dark_theme()
}

create_comprehensive_protein_table <- function(diff_lists, protein_info = NULL) {
  all_proteins <- unique(unlist(diff_lists))
  cat("Total unique proteins across all factors:", length(all_proteins), "\n")
  
  protein_table <- data.frame(Protein_ID = all_proteins, stringsAsFactors = FALSE)
  
  for (factor_name in names(diff_lists)) {
    protein_table[[factor_name]] <- protein_table$Protein_ID %in% diff_lists[[factor_name]]
  }
  
  protein_table$Total_Factors <- rowSums(protein_table[, names(diff_lists), drop = FALSE])
  
  if (!is.null(protein_info)) {
    protein_table <- merge(protein_table, protein_info, by.x = "Protein_ID", by.y = "Protein IDs", all.x = TRUE)
  }
  
  protein_table[order(-protein_table$Total_Factors), ]
}

analyze_specific_overlaps <- function(diff_lists) {
  specific_overlaps <- list()
  factor_names <- names(diff_lists)
  
  specific_pairs <- list(
    c("Sex", "Age"), c("Sex", "Days_Until_Relapse"), c("Sex", "Relapse_Risk"),
    c("Sex", "Subject"), c("Sex", "Treatment"), c("Age", "Relapse_Risk"),
    c("Age", "Days_Until_Relapse"), c("Relapse_Risk", "Days_Until_Relapse"))
  
  for (pair in specific_pairs) {
    if (all(pair %in% factor_names)) {
      overlap_proteins <- intersect(diff_lists[[pair[1]]], diff_lists[[pair[2]]])
      specific_overlaps[[paste(pair, collapse = " and ")]] <- overlap_proteins
    }
  }
  
  if (all(c("Sex", "Age", "Days_Until_Relapse") %in% factor_names)) {
    all_clinical <- Reduce(intersect, diff_lists[c("Sex", "Age", "Days_Until_Relapse")])
    specific_overlaps[["All Clinical"]] <- all_clinical
  }
  
  comprehensive_table <- create_comprehensive_protein_table(diff_lists)
  specific_overlaps[["Core 3plus factors"]] <- comprehensive_table$Protein_ID[comprehensive_table$Total_Factors >= 3]
  specific_overlaps[["Core 4plus factors"]] <- comprehensive_table$Protein_ID[comprehensive_table$Total_Factors >= 4]
  
  specific_overlaps
}

save_overlap_results <- function(diff_lists, overlap_results, comprehensive_table, specific_overlaps) {
  fwrite(comprehensive_table, file.path(OUTPUT_DIR, "comprehensive_protein_table.csv"))
  fwrite(overlap_results$details, file.path(OUTPUT_DIR, "pairwise_overlaps.csv"))
  fwrite(as.data.frame(overlap_results$matrix), file.path(OUTPUT_DIR, "overlap_matrix.csv"), row.names = TRUE)
  
  factor_summary <- data.frame(
    Factor = names(diff_lists),
    Protein_Count = sapply(diff_lists, length),
    Percentage_of_Total = round(100 * sapply(diff_lists, length) / length(unique(unlist(diff_lists))), 2),
    stringsAsFactors = FALSE)
  fwrite(factor_summary, file.path(OUTPUT_DIR, "factor_summary_stats.csv"))
  
  specific_summary <- data.frame(
    Overlap_Type = names(specific_overlaps),
    Protein_Count = sapply(specific_overlaps, length),
    stringsAsFactors = FALSE)
  fwrite(specific_summary, file.path(OUTPUT_DIR, "specific_overlaps_summary.csv"))
  
  for (overlap_name in names(specific_overlaps)) {
    if (length(specific_overlaps[[overlap_name]]) > 0) {
      filename <- paste0("proteins_", gsub("[^A-Za-z0-9]", "_", overlap_name), ".csv")
      protein_df <- data.frame(Protein_ID = specific_overlaps[[overlap_name]], stringsAsFactors = FALSE)
      fwrite(protein_df, file.path(OUTPUT_DIR, filename))
    }
  }
}

create_overlap_report <- function(diff_lists, overlap_results, comprehensive_table, specific_overlaps) {
  total_unique <- length(unique(unlist(diff_lists)))
  factor_summary <- data.frame(Factor = names(diff_lists), Count = sapply(diff_lists, length))
  factor_summary <- factor_summary[order(-factor_summary$Count), ]
  
  top_overlaps <- head(overlap_results$details[order(-overlap_results$details$Overlap_Count), ], 10)
  top_jaccard <- head(overlap_results$details[order(-overlap_results$details$Jaccard_Index), ], 5)
  
  multi_factor_counts <- table(comprehensive_table$Total_Factors)
  
  report_lines <- c(
    "# Protein Overlap Analysis Report", "", "## Executive Summary", "",
    paste("- **Total unique proteins**:", total_unique),
    paste("- **Factors analyzed**:", length(diff_lists)),
    paste("- **Pairwise comparisons**:", nrow(overlap_results$details)), "", "## Factor Summary", ""
  )
  
  for (i in 1:nrow(factor_summary)) {
    percentage <- round(100 * factor_summary$Count[i] / total_unique, 2)
    report_lines <- c(report_lines, paste("- **", factor_summary$Factor[i], "**:", factor_summary$Count[i],
                                         "proteins (", percentage, "% of total unique proteins)"))
  }
  
  report_lines <- c(report_lines, "", "## Top Overlaps (by count)", "")
  for (i in 1:min(10, nrow(top_overlaps))) {
    row <- top_overlaps[i, ]
    perc1 <- round(row$Overlap_Percentage_F1, 1)
    perc2 <- round(row$Overlap_Percentage_F2, 1)
    report_lines <- c(report_lines, paste("- **", row$Factor1, " ∩ ", row$Factor2, "**:", row$Overlap_Count,
                                         "proteins (", perc1, "% of ", row$Factor1, ", ", perc2, "% of ", row$Factor2, ")"))
  }
  
  report_lines <- c(report_lines, "", "## Highest Similarity (by Jaccard Index)", "")
  for (i in 1:min(5, nrow(top_jaccard))) {
    row <- top_jaccard[i, ]
    report_lines <- c(report_lines, paste("- **", row$Factor1, " ∩ ", row$Factor2, "**: Jaccard =",
                                         round(row$Jaccard_Index, 3), " (", row$Overlap_Count, "proteins)"))
  }
  
  report_lines <- c(report_lines, "", "## Specific Overlaps of Interest", "")
  for (overlap_name in names(specific_overlaps)) {
    count <- length(specific_overlaps[[overlap_name]])
    report_lines <- c(report_lines, paste("- **", overlap_name, "**:", count, "proteins"))
  }
  
  report_lines <- c(report_lines, "", "## Core Proteins (Multi-Factor)", "")
  for (factors in names(multi_factor_counts)[names(multi_factor_counts) >= 3]) {
    if (as.numeric(factors) >= 3) {
      count <- multi_factor_counts[factors]
      report_lines <- c(report_lines, paste("- **", factors, "+ factors**:", count, "proteins"))
    }
  }
  
  writeLines(report_lines, file.path(OUTPUT_DIR, "OVERLAP_ANALYSIS_REPORT.md"))
  cat("Overlap analysis report saved to:", file.path(OUTPUT_DIR, "OVERLAP_ANALYSIS_REPORT.md"), "\n")
}

main <- function() {
  create_output_dir()
  
  print_section("PROTEIN OVERLAP ANALYSIS")
  print_section("LOADING DIFFERENTIAL PROTEINS")
  
  differential_lists <- load_differential_proteins(RESULTS_DIR, "corrected")
  significant_factors <- differential_lists[sapply(differential_lists, length) > 0]
  
  if (length(significant_factors) == 0) {
    stop("No factors with significant proteins found!")
  }
  
  cat("Factors with significant proteins:\n")
  for (factor in names(significant_factors)) {
    cat("  -", factor, ":", length(significant_factors[[factor]]), "proteins\n")
  }
  
  print_section("CALCULATING OVERLAPS")
  overlap_results <- calculate_pairwise_overlaps(significant_factors)
  
  print_section("CREATING VISUALIZATIONS")
  
  plot1 <- create_factor_sizes_plot(significant_factors)
  ggsave(file.path(OUTPUT_DIR, "factor_sizes.png"), plot1, width = 12, height = 8, dpi = 300)
  
  plot2 <- create_overlap_heatmap(overlap_results$matrix, "Protein Overlap Counts")
  ggsave(file.path(OUTPUT_DIR, "overlap_heatmap.png"), plot2, width = 10, height = 8, dpi = 300)
  
  plot3 <- create_jaccard_heatmap(overlap_results$details, "Jaccard Similarity Index")
  ggsave(file.path(OUTPUT_DIR, "jaccard_heatmap.png"), plot3, width = 10, height = 8, dpi = 300)
  
  plot4 <- create_overlap_summary_plot(overlap_results$details)
  ggsave(file.path(OUTPUT_DIR, "overlap_summary.png"), plot4, width = 12, height = 8, dpi = 300)
  
  print_section("CREATING COMPREHENSIVE PROTEIN TABLE")
  comprehensive_table <- create_comprehensive_protein_table(significant_factors)
  
  print_section("ANALYZING SPECIFIC OVERLAPS")
  specific_overlaps <- analyze_specific_overlaps(significant_factors)
  
  print_section("SAVING RESULTS")
  save_overlap_results(significant_factors, overlap_results, comprehensive_table, specific_overlaps)
  create_overlap_report(significant_factors, overlap_results, comprehensive_table, specific_overlaps)
  
  print_section("OVERLAP ANALYSIS COMPLETE")
  cat("Results saved to:", OUTPUT_DIR, "\n")
  cat("Key files generated:\n")
  cat("  - Visualizations: factor_sizes.png, overlap_heatmap.png, jaccard_heatmap.png\n")
  cat("  - Comprehensive table: comprehensive_protein_table.csv\n")
  cat("  - Overlap details: pairwise_overlaps.csv\n")
  cat("  - Specific overlaps: specific_overlaps_summary.csv + individual files\n")
  cat("  - Summary report: OVERLAP_ANALYSIS_REPORT.md\n")
  
  return(TRUE)
}

if (!interactive()) {
  main()
}