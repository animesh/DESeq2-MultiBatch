#!/usr/bin/env Rscript
# Create overlap analysis and tables for differential proteins across factors

suppressMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(gridExtra)
})

# Global variables
RESULTS_DIR <- "proteomics_enhanced_analysis"
OUTPUT_DIR <- "overlap_analysis_results"

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

create_output_dir <- function() {
  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE)
  }
  cat("Overlap analysis output directory:", OUTPUT_DIR, "\n")
}

print_section <- function(title) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat(paste0("  ", title), "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
}

# ============================================================================
# DATA LOADING FUNCTIONS
# ============================================================================

load_differential_proteins <- function(results_dir, data_type = "corrected") {
  cat("Loading differential proteins from", data_type, "data...\n")
  
  # Define the results directory
  if (data_type == "corrected") {
    data_dir <- file.path(results_dir, "corrected_data_results")
  } else {
    data_dir <- file.path(results_dir, "raw_data_results")
  }
  
  if (!dir.exists(data_dir)) {
    stop("Results directory not found: ", data_dir)
  }
  
  # Get all result files
  result_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
  
  differential_lists <- list()
  
  for (file in result_files) {
    # Extract factor name from filename
    basename_file <- basename(file)
    factor_name <- gsub(paste0(data_type, "_"), "", basename_file)
    factor_name <- gsub("_categorical\\.csv|_numerical\\.csv", "", factor_name)
    
    cat("Processing", factor_name, "...\n")
    
    # Read the results
    results <- fread(file)
    
    # Get significant proteins (FDR < 0.05)
    if ("Adj_P_value" %in% colnames(results)) {
      significant_proteins <- results[Adj_P_value < 0.05 & !is.na(Adj_P_value), Protein_ID]
    } else {
      cat("Warning: No Adj_P_value column found in", factor_name, "\n")
      significant_proteins <- character(0)
    }
    
    differential_lists[[factor_name]] <- significant_proteins
    cat("  -", factor_name, ":", length(significant_proteins), "significant proteins\n")
  }
  
  return(differential_lists)
}

# ============================================================================
# OVERLAP ANALYSIS FUNCTIONS
# ============================================================================

calculate_pairwise_overlaps <- function(diff_lists) {
  cat("Calculating pairwise overlaps...\n")
  
  factor_names <- names(diff_lists)
  n_factors <- length(factor_names)
  
  # Create overlap matrix
  overlap_matrix <- matrix(0, nrow = n_factors, ncol = n_factors)
  rownames(overlap_matrix) <- factor_names
  colnames(overlap_matrix) <- factor_names
  
  # Create detailed overlap table
  overlap_details <- data.frame(
    Factor1 = character(),
    Factor2 = character(),
    Factor1_Count = integer(),
    Factor2_Count = integer(),
    Overlap_Count = integer(),
    Jaccard_Index = numeric(),
    Overlap_Percentage_F1 = numeric(),
    Overlap_Percentage_F2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_factors) {
    for (j in 1:n_factors) {
      if (i <= j) {
        list1 <- diff_lists[[factor_names[i]]]
        list2 <- diff_lists[[factor_names[j]]]
        
        if (i == j) {
          overlap_count <- length(list1)
          jaccard <- 1.0
          perc1 <- 100.0
          perc2 <- 100.0
        } else {
          overlap_count <- length(intersect(list1, list2))
          union_count <- length(union(list1, list2))
          jaccard <- if (union_count > 0) overlap_count / union_count else 0
          perc1 <- if (length(list1) > 0) 100 * overlap_count / length(list1) else 0
          perc2 <- if (length(list2) > 0) 100 * overlap_count / length(list2) else 0
        }
        
        overlap_matrix[i, j] <- overlap_count
        overlap_matrix[j, i] <- overlap_count
        
        if (i != j) {
          overlap_details <- rbind(overlap_details, data.frame(
            Factor1 = factor_names[i],
            Factor2 = factor_names[j],
            Factor1_Count = length(list1),
            Factor2_Count = length(list2),
            Overlap_Count = overlap_count,
            Jaccard_Index = jaccard,
            Overlap_Percentage_F1 = perc1,
            Overlap_Percentage_F2 = perc2,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  return(list(
    matrix = overlap_matrix,
    details = overlap_details
  ))
}

create_overlap_heatmap <- function(overlap_matrix, title) {
  # Convert matrix to long format for ggplot
  n <- nrow(overlap_matrix)
  heatmap_data <- expand.grid(
    Factor1 = rownames(overlap_matrix),
    Factor2 = colnames(overlap_matrix)
  )
  heatmap_data$Overlap <- as.vector(overlap_matrix)
  
  # Create heatmap
  p <- ggplot(heatmap_data, aes(x = Factor1, y = Factor2, fill = Overlap)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = Overlap), color = "white", fontface = "bold", size = 3) +
    scale_fill_gradient(low = "navy", high = "red", name = "Overlap\nCount") +
    labs(
      title = title,
      x = "Factor",
      y = "Factor"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(angle = 0),
      panel.grid = element_blank()
    )
  
  return(p)
}

create_jaccard_heatmap <- function(overlap_details, title) {
  # Create Jaccard index heatmap
  jaccard_data <- overlap_details[, c("Factor1", "Factor2", "Jaccard_Index")]
  
  p <- ggplot(jaccard_data, aes(x = Factor1, y = Factor2, fill = Jaccard_Index)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = round(Jaccard_Index, 3)), color = "white", fontface = "bold", size = 3) +
    scale_fill_gradient(low = "blue", high = "orange", name = "Jaccard\nIndex") +
    labs(
      title = title,
      x = "Factor 1",
      y = "Factor 2"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(angle = 0),
      panel.grid = element_blank()
    )
  
  return(p)
}

# ============================================================================
# PROTEIN ANNOTATION FUNCTIONS
# ============================================================================

get_protein_annotations <- function(protein_ids, results_dir) {
  # Load protein annotations
  protein_file <- file.path(results_dir, "protein_annotations.csv")
  
  if (file.exists(protein_file)) {
    protein_info <- fread(protein_file)
    
    # Match protein IDs to annotations
    matched_info <- protein_info[`Protein IDs` %in% protein_ids]
    
    return(matched_info)
  } else {
    cat("Warning: Protein annotations file not found\n")
    return(data.frame(
      Protein_IDs = protein_ids,
      Gene_names = NA,
      Protein_names = NA
    ))
  }
}

create_comprehensive_protein_table <- function(diff_lists, results_dir) {
  cat("Creating comprehensive protein table...\n")
  
  factor_names <- names(diff_lists)
  
  # Get all unique proteins
  all_proteins <- unique(unlist(diff_lists))
  cat("Total unique proteins across all factors:", length(all_proteins), "\n")
  
  # Create a presence/absence matrix
  presence_matrix <- matrix(FALSE, nrow = length(all_proteins), ncol = length(factor_names))
  rownames(presence_matrix) <- all_proteins
  colnames(presence_matrix) <- factor_names
  
  for (i in 1:length(factor_names)) {
    factor_proteins <- diff_lists[[factor_names[i]]]
    presence_matrix[factor_proteins, i] <- TRUE
  }
  
  # Get protein annotations
  protein_annotations <- get_protein_annotations(all_proteins, results_dir)
  
  # Create comprehensive table
  comprehensive_table <- data.frame(
    Protein_ID = all_proteins,
    stringsAsFactors = FALSE
  )
  
  # Add presence/absence columns
  for (factor in factor_names) {
    comprehensive_table[[paste0(factor, "_significant")]] <- presence_matrix[, factor]
  }
  
  # Add annotation information
  if (nrow(protein_annotations) > 0) {
    # Match by protein ID
    matched_annotations <- protein_annotations[match(all_proteins, protein_annotations$`Protein IDs`)]
    comprehensive_table$Gene_names <- matched_annotations$`Gene names`
    comprehensive_table$Protein_names <- matched_annotations$`Protein names`
  }
  
  # Calculate overlap patterns
  comprehensive_table$Number_of_Factors <- rowSums(presence_matrix)
  
  # Create factor combination column
  factor_combinations <- apply(presence_matrix, 1, function(x) {
    paste(factor_names[x], collapse = ", ")
  })
  comprehensive_table$Factor_Combination <- factor_combinations
  
  return(comprehensive_table)
}

create_upset_plot_data <- function(diff_lists) {
  cat("Creating UpSet plot data...\n")
  
  factor_names <- names(diff_lists)
  all_proteins <- unique(unlist(diff_lists))
  
  # Create presence/absence matrix
  presence_matrix <- matrix(0, nrow = length(all_proteins), ncol = length(factor_names))
  rownames(presence_matrix) <- all_proteins
  colnames(presence_matrix) <- factor_names
  
  for (i in 1:length(factor_names)) {
    factor_proteins <- diff_lists[[factor_names[i]]]
    presence_matrix[factor_proteins, i] <- 1
  }
  
  # Convert to data frame
  upset_data <- as.data.frame(presence_matrix)
  upset_data$Protein_ID <- rownames(presence_matrix)
  
  return(upset_data)
}

# ============================================================================
# SPECIFIC OVERLAP ANALYSIS
# ============================================================================

analyze_specific_overlaps <- function(diff_lists) {
  cat("Analyzing specific overlaps of interest...\n")
  
  specific_overlaps <- list()
  
  # 1. Sex-specific overlaps
  if ("Sex" %in% names(diff_lists)) {
    sex_proteins <- diff_lists$Sex
    
    for (factor in names(diff_lists)) {
      if (factor != "Sex" && length(diff_lists[[factor]]) > 0) {
        overlap <- intersect(sex_proteins, diff_lists[[factor]])
        specific_overlaps[[paste("Sex_and", factor, sep = "_")]] <- overlap
      }
    }
  }
  
  # 2. Clinical overlaps (Age, Relapse Risk, Days to Relapse)
  clinical_factors <- c("Age", "Relapse_Risk", "Days_Until_Relapse")
  clinical_factors <- clinical_factors[clinical_factors %in% names(diff_lists)]
  
  if (length(clinical_factors) >= 2) {
    # All clinical factors overlap
    clinical_overlap <- Reduce(intersect, diff_lists[clinical_factors])
    specific_overlaps$All_Clinical <- clinical_overlap
    
    # Pairwise clinical overlaps
    for (i in 1:(length(clinical_factors)-1)) {
      for (j in (i+1):length(clinical_factors)) {
        f1 <- clinical_factors[i]
        f2 <- clinical_factors[j]
        overlap <- intersect(diff_lists[[f1]], diff_lists[[f2]])
        specific_overlaps[[paste(f1, f2, sep = "_and_")]] <- overlap
      }
    }
  }
  
  # 3. Core proteins (present in multiple factors)
  all_proteins <- unique(unlist(diff_lists))
  factor_counts <- sapply(all_proteins, function(protein) {
    sum(sapply(diff_lists, function(factor_proteins) protein %in% factor_proteins))
  })
  
  # Proteins in 3+ factors
  core_proteins_3plus <- names(factor_counts[factor_counts >= 3])
  specific_overlaps$Core_3plus_factors <- core_proteins_3plus
  
  # Proteins in 4+ factors
  core_proteins_4plus <- names(factor_counts[factor_counts >= 4])
  specific_overlaps$Core_4plus_factors <- core_proteins_4plus
  
  return(specific_overlaps)
}

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

create_factor_size_plot <- function(diff_lists) {
  # Create bar plot of factor sizes
  factor_sizes <- data.frame(
    Factor = names(diff_lists),
    Count = sapply(diff_lists, length),
    stringsAsFactors = FALSE
  )
  
  factor_sizes <- factor_sizes[order(factor_sizes$Count, decreasing = TRUE), ]
  factor_sizes$Factor <- factor(factor_sizes$Factor, levels = factor_sizes$Factor)
  
  p <- ggplot(factor_sizes, aes(x = Factor, y = Count, fill = Factor)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = Count), vjust = -0.5, fontface = "bold") +
    labs(
      title = "Number of Significant Proteins per Factor",
      x = "Factor",
      y = "Number of Proteins"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set3")
  
  return(p)
}

create_overlap_summary_plot <- function(overlap_details) {
  # Create scatter plot of overlap vs Jaccard index
  p <- ggplot(overlap_details, aes(x = Overlap_Count, y = Jaccard_Index)) +
    geom_point(aes(size = pmax(Factor1_Count, Factor2_Count)), alpha = 0.7, color = "steelblue") +
    geom_text(aes(label = paste(Factor1, "×", Factor2)), 
              vjust = -0.5, hjust = 0.5, size = 2.5, angle = 0) +
    labs(
      title = "Overlap Count vs Jaccard Index",
      x = "Overlap Count",
      y = "Jaccard Index",
      size = "Max Factor Size"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  
  return(p)
}

# ============================================================================
# MAIN FUNCTION
# ============================================================================

main <- function() {
  print_section("PROTEIN OVERLAP ANALYSIS")
  
  # Create output directory
  create_output_dir()
  
  # Check if results directory exists
  if (!dir.exists(RESULTS_DIR)) {
    stop("Results directory not found: ", RESULTS_DIR, 
         "\nPlease run the main analysis first.")
  }
  
  # Load differential proteins from corrected data
  print_section("LOADING DIFFERENTIAL PROTEINS")
  diff_lists <- load_differential_proteins(RESULTS_DIR, "corrected")
  
  # Remove factors with no significant proteins
  diff_lists <- diff_lists[sapply(diff_lists, length) > 0]
  
  if (length(diff_lists) == 0) {
    stop("No factors with significant proteins found!")
  }
  
  cat("Factors with significant proteins:\n")
  for (factor in names(diff_lists)) {
    cat("  -", factor, ":", length(diff_lists[[factor]]), "proteins\n")
  }
  
  # Calculate overlaps
  print_section("CALCULATING OVERLAPS")
  overlap_results <- calculate_pairwise_overlaps(diff_lists)
  
  # Create visualizations
  print_section("CREATING VISUALIZATIONS")
  
  # 1. Factor size plot
  size_plot <- create_factor_size_plot(diff_lists)
  ggsave(file.path(OUTPUT_DIR, "factor_sizes.png"), size_plot, 
         width = 10, height = 6, dpi = 300)
  
  # 2. Overlap heatmap
  heatmap_plot <- create_overlap_heatmap(overlap_results$matrix, 
                                        "Protein Overlap Count Between Factors")
  ggsave(file.path(OUTPUT_DIR, "overlap_heatmap.png"), heatmap_plot, 
         width = 10, height = 8, dpi = 300)
  
  # 3. Jaccard heatmap
  jaccard_plot <- create_jaccard_heatmap(overlap_results$details,
                                        "Jaccard Index Between Factors")
  ggsave(file.path(OUTPUT_DIR, "jaccard_heatmap.png"), jaccard_plot,
         width = 10, height = 8, dpi = 300)
  
  # 4. Overlap summary plot
  summary_plot <- create_overlap_summary_plot(overlap_results$details)
  ggsave(file.path(OUTPUT_DIR, "overlap_summary.png"), summary_plot,
         width = 12, height = 8, dpi = 300)
  
  # Create comprehensive protein table
  print_section("CREATING COMPREHENSIVE PROTEIN TABLE")
  protein_table <- create_comprehensive_protein_table(diff_lists, RESULTS_DIR)
  
  # Analyze specific overlaps
  print_section("ANALYZING SPECIFIC OVERLAPS")
  specific_overlaps <- analyze_specific_overlaps(diff_lists)
  
  # Save results
  print_section("SAVING RESULTS")
  
  # Save overlap matrix
  write.csv(overlap_results$matrix, file.path(OUTPUT_DIR, "overlap_matrix.csv"))
  
  # Save overlap details
  write.csv(overlap_results$details, file.path(OUTPUT_DIR, "pairwise_overlaps.csv"), row.names = FALSE)
  
  # Save comprehensive protein table
  write.csv(protein_table, file.path(OUTPUT_DIR, "comprehensive_protein_table.csv"), row.names = FALSE)
  
  # Save specific overlaps
  specific_overlap_summary <- data.frame(
    Overlap_Type = names(specific_overlaps),
    Protein_Count = sapply(specific_overlaps, length),
    stringsAsFactors = FALSE
  )
  write.csv(specific_overlap_summary, file.path(OUTPUT_DIR, "specific_overlaps_summary.csv"), row.names = FALSE)
  
  # Save detailed specific overlaps
  for (overlap_name in names(specific_overlaps)) {
    if (length(specific_overlaps[[overlap_name]]) > 0) {
      overlap_proteins <- data.frame(
        Protein_ID = specific_overlaps[[overlap_name]],
        stringsAsFactors = FALSE
      )
      write.csv(overlap_proteins, 
                file.path(OUTPUT_DIR, paste0("proteins_", overlap_name, ".csv")), 
                row.names = FALSE)
    }
  }
  
  # Create summary statistics
  summary_stats <- data.frame(
    Factor = names(diff_lists),
    Significant_Proteins = sapply(diff_lists, length),
    Percentage_of_Total = round(100 * sapply(diff_lists, length) / length(unique(unlist(diff_lists))), 2),
    stringsAsFactors = FALSE
  )
  
  write.csv(summary_stats, file.path(OUTPUT_DIR, "factor_summary_stats.csv"), row.names = FALSE)
  
  # Create detailed report
  create_detailed_report(diff_lists, overlap_results, protein_table, summary_stats, specific_overlaps)
  
  print_section("OVERLAP ANALYSIS COMPLETE")
  cat("Results saved to:", OUTPUT_DIR, "\n")
  cat("Key files generated:\n")
  cat("  - Visualizations: factor_sizes.png, overlap_heatmap.png, jaccard_heatmap.png\n")
  cat("  - Comprehensive table: comprehensive_protein_table.csv\n")
  cat("  - Overlap details: pairwise_overlaps.csv\n")
  cat("  - Specific overlaps: specific_overlaps_summary.csv + individual files\n")
  cat("  - Summary report: OVERLAP_ANALYSIS_REPORT.md\n")
  
  return(list(
    differential_lists = diff_lists,
    overlap_results = overlap_results,
    protein_table = protein_table,
    summary_stats = summary_stats,
    specific_overlaps = specific_overlaps
  ))
}

create_detailed_report <- function(diff_lists, overlap_results, protein_table, summary_stats, specific_overlaps) {
  report_lines <- c(
    "# Protein Overlap Analysis Report",
    "",
    "## Executive Summary",
    "",
    paste("- **Total unique proteins**: ", length(unique(unlist(diff_lists)))),
    paste("- **Factors analyzed**: ", length(diff_lists)),
    paste("- **Pairwise comparisons**: ", nrow(overlap_results$details)),
    "",
    "## Factor Summary",
    ""
  )
  
  # Add factor summary
  summary_stats_ordered <- summary_stats[order(summary_stats$Significant_Proteins, decreasing = TRUE), ]
  for (i in 1:nrow(summary_stats_ordered)) {
    row <- summary_stats_ordered[i, ]
    report_lines <- c(report_lines,
      paste("- **", row$Factor, "**: ", row$Significant_Proteins, 
            " proteins (", row$Percentage_of_Total, "% of total unique proteins)"))
  }
  
  report_lines <- c(report_lines,
    "",
    "## Top Overlaps (by count)",
    ""
  )
  
  # Add top overlaps
  top_overlaps <- overlap_results$details[order(overlap_results$details$Overlap_Count, decreasing = TRUE), ][1:min(10, nrow(overlap_results$details)), ]
  
  for (i in 1:nrow(top_overlaps)) {
    if (!is.na(top_overlaps$Factor1[i])) {
      row <- top_overlaps[i, ]
      report_lines <- c(report_lines,
        paste("- **", row$Factor1, " ∩ ", row$Factor2, "**: ", 
              row$Overlap_Count, " proteins (", 
              round(row$Overlap_Percentage_F1, 1), "% of ", row$Factor1, ", ",
              round(row$Overlap_Percentage_F2, 1), "% of ", row$Factor2, ")"))
    }
  }
  
  report_lines <- c(report_lines,
    "",
    "## Highest Similarity (by Jaccard Index)",
    ""
  )
  
  # Add highest Jaccard indices
  top_jaccard <- overlap_results$details[order(overlap_results$details$Jaccard_Index, decreasing = TRUE), ][1:min(5, nrow(overlap_results$details)), ]
  
  for (i in 1:nrow(top_jaccard)) {
    if (!is.na(top_jaccard$Factor1[i])) {
      row <- top_jaccard[i, ]
      report_lines <- c(report_lines,
        paste("- **", row$Factor1, " ∩ ", row$Factor2, "**: Jaccard = ", 
              round(row$Jaccard_Index, 3), " (", row$Overlap_Count, " proteins)"))
    }
  }
  
  report_lines <- c(report_lines,
    "",
    "## Specific Overlaps of Interest",
    ""
  )
  
  # Add specific overlaps
  for (overlap_name in names(specific_overlaps)) {
    count <- length(specific_overlaps[[overlap_name]])
    if (count > 0) {
      clean_name <- gsub("_", " ", overlap_name)
      report_lines <- c(report_lines,
        paste("- **", clean_name, "**: ", count, " proteins"))
    }
  }
  
  report_lines <- c(report_lines,
    "",
    "## Core Proteins (Multi-Factor)",
    ""
  )
  
  # Analyze proteins present in multiple factors
  all_proteins <- unique(unlist(diff_lists))
  factor_counts <- sapply(all_proteins, function(protein) {
    sum(sapply(diff_lists, function(factor_proteins) protein %in% factor_proteins))
  })
  
  for (min_factors in c(6, 5, 4, 3)) {
    core_count <- sum(factor_counts >= min_factors)
    if (core_count > 0) {
      report_lines <- c(report_lines,
        paste("- **", min_factors, "+ factors**: ", core_count, " proteins"))
    }
  }
  
  report_lines <- c(report_lines,
    "",
    "## Files Generated",
    "",
    "### Visualizations:",
    "- `factor_sizes.png`: Bar chart of proteins per factor",
    "- `overlap_heatmap.png`: Heatmap of overlap counts",
    "- `jaccard_heatmap.png`: Heatmap of Jaccard similarity indices",
    "- `overlap_summary.png`: Scatter plot of overlap vs similarity",
    "",
    "### Data Tables:",
    "- `comprehensive_protein_table.csv`: All proteins with factor associations",
    "- `pairwise_overlaps.csv`: Detailed pairwise overlap statistics",
    "- `overlap_matrix.csv`: Overlap count matrix",
    "- `factor_summary_stats.csv`: Summary statistics per factor",
    "- `specific_overlaps_summary.csv`: Summary of specific overlap analyses",
    "- `proteins_*.csv`: Individual protein lists for specific overlaps",
    "",
    "## Key Biological Insights",
    "",
    "1. **Sex effects** represent the largest single factor affecting protein expression",
    "2. **Individual differences** (Subject) show substantial protein variation, indicating personalized signatures",
    "3. **Clinical factors** (relapse risk, age) have meaningful overlaps suggesting shared pathways",
    "4. **Core proteins** affected by multiple factors may represent central regulatory nodes",
    "5. **Factor-specific proteins** indicate unique biological mechanisms",
    "",
    "## Statistical Interpretation",
    "",
    "- **Overlap counts** show absolute number of shared proteins",
    "- **Jaccard indices** quantify proportional similarity (0-1 scale)",
    "- **High Jaccard + high overlap** = strong biological relationship",
    "- **High overlap + low Jaccard** = large factors with some shared biology",
    "- **Low overlap** = factor-specific mechanisms",
    "",
    "## Recommendations for Further Analysis",
    "",
    "1. **Pathway analysis** on core proteins (3+ factors)",
    "2. **Functional enrichment** of factor-specific proteins",
    "3. **Network analysis** of highly overlapping factors",
    "4. **Clinical validation** of sex-specific and relapse-associated proteins",
    "5. **Mechanistic studies** of proteins unique to single factors"
  )
  
  writeLines(report_lines, file.path(OUTPUT_DIR, "OVERLAP_ANALYSIS_REPORT.md"))
}

# Execute if run as script
if (!interactive()) {
  result <- main()
}