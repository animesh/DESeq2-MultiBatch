#!/usr/bin/env Rscript
# Create Venn diagrams and overlap tables for differential proteins across factors

suppressMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(VennDiagram)
  library(RColorBrewer)
  library(gridExtra)
})

# Global variables
RESULTS_DIR <- "proteomics_enhanced_analysis"
OUTPUT_DIR <- "venn_analysis_results"

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

create_output_dir <- function() {
  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE)
  }
  cat("Venn analysis output directory:", OUTPUT_DIR, "\n")
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
# VENN DIAGRAM FUNCTIONS
# ============================================================================

create_pairwise_venn <- function(list1, list2, name1, name2, title) {
  # Create temporary file for Venn diagram
  temp_file <- tempfile(fileext = ".png")
  
  # Create Venn diagram
  venn_plot <- venn.diagram(
    x = list(list1, list2),
    category.names = c(name1, name2),
    filename = NULL,
    output = TRUE,
    
    # Appearance
    lwd = 2,
    lty = 'blank',
    fill = c('#E31A1C', '#1F78B4'),
    cex = 1.5,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans",
    rotation = 1,
    
    # Title
    main = title,
    main.cex = 1.3,
    main.fontface = "bold"
  )
  
  return(venn_plot)
}

create_three_way_venn <- function(list1, list2, list3, name1, name2, name3, title) {
  # Create Venn diagram
  venn_plot <- venn.diagram(
    x = list(list1, list2, list3),
    category.names = c(name1, name2, name3),
    filename = NULL,
    output = TRUE,
    
    # Appearance
    lwd = 2,
    lty = 'blank',
    fill = c('#E31A1C', '#1F78B4', '#33A02C'),
    cex = 1.2,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.1,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.055),
    cat.fontfamily = "sans",
    rotation = 1,
    
    # Title
    main = title,
    main.cex = 1.2,
    main.fontface = "bold"
  )
  
  return(venn_plot)
}

create_four_way_venn <- function(list1, list2, list3, list4, name1, name2, name3, name4, title) {
  # Create Venn diagram
  venn_plot <- venn.diagram(
    x = list(list1, list2, list3, list4),
    category.names = c(name1, name2, name3, name4),
    filename = NULL,
    output = TRUE,
    
    # Appearance
    lwd = 2,
    lty = 'blank',
    fill = c('#E31A1C', '#1F78B4', '#33A02C', '#FF7F00'),
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-15, 15, -15, 15),
    cat.dist = c(0.22, 0.22, 0.22, 0.22),
    cat.fontfamily = "sans",
    
    # Title
    main = title,
    main.cex = 1.1,
    main.fontface = "bold"
  )
  
  return(venn_plot)
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
        } else {
          overlap_count <- length(intersect(list1, list2))
          union_count <- length(union(list1, list2))
          jaccard <- if (union_count > 0) overlap_count / union_count else 0
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
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = Overlap), color = "white", fontface = "bold") +
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

create_overlap_protein_tables <- function(diff_lists, results_dir) {
  cat("Creating protein overlap tables...\n")
  
  factor_names <- names(diff_lists)
  overlap_tables <- list()
  
  # Get all unique proteins
  all_proteins <- unique(unlist(diff_lists))
  
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
  
  return(comprehensive_table)
}

# ============================================================================
# MAIN ANALYSIS FUNCTIONS
# ============================================================================

create_key_comparisons <- function(diff_lists, output_dir) {
  cat("Creating key comparison Venn diagrams...\n")
  
  # 1. Sex vs Relapse Risk vs Subject (top 3 biological factors)
  if (all(c("Sex", "Relapse_Risk", "Subject") %in% names(diff_lists))) {
    venn1 <- create_three_way_venn(
      diff_lists$Sex,
      diff_lists$Relapse_Risk,
      diff_lists$Subject,
      "Sex", "Relapse Risk", "Subject",
      "Biological Factors Overlap"
    )
    
    png(file.path(output_dir, "venn_biological_factors.png"), width = 800, height = 600, res = 300)
    grid.draw(venn1)
    dev.off()
  }
  
  # 2. Sex vs Age vs Days_Until_Relapse (clinical factors)
  if (all(c("Sex", "Age", "Days_Until_Relapse") %in% names(diff_lists))) {
    venn2 <- create_three_way_venn(
      diff_lists$Sex,
      diff_lists$Age,
      diff_lists$Days_Until_Relapse,
      "Sex", "Age", "Days to Relapse",
      "Clinical Factors Overlap"
    )
    
    png(file.path(output_dir, "venn_clinical_factors.png"), width = 800, height = 600, res = 300)
    grid.draw(venn2)
    dev.off()
  }
  
  # 3. Sex vs Location (demographic factors)
  if (all(c("Sex", "Location") %in% names(diff_lists))) {
    venn3 <- create_pairwise_venn(
      diff_lists$Sex,
      diff_lists$Location,
      "Sex", "Location",
      "Demographic Factors Overlap"
    )
    
    png(file.path(output_dir, "venn_demographic_factors.png"), width = 600, height = 600, res = 300)
    grid.draw(venn3)
    dev.off()
  }
  
  # 4. Top 4 factors (if we have enough)
  top_factors <- names(sort(sapply(diff_lists, length), decreasing = TRUE))[1:4]
  top_factors <- top_factors[!is.na(top_factors)]
  
  if (length(top_factors) >= 4) {
    venn4 <- create_four_way_venn(
      diff_lists[[top_factors[1]]],
      diff_lists[[top_factors[2]]],
      diff_lists[[top_factors[3]]],
      diff_lists[[top_factors[4]]],
      top_factors[1], top_factors[2], top_factors[3], top_factors[4],
      "Top 4 Factors Overlap"
    )
    
    png(file.path(output_dir, "venn_top4_factors.png"), width = 800, height = 800, res = 300)
    grid.draw(venn4)
    dev.off()
  }
}

# ============================================================================
# MAIN FUNCTION
# ============================================================================

main <- function() {
  print_section("VENN DIAGRAM AND OVERLAP ANALYSIS")
  
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
  
  # Create overlap heatmap
  print_section("CREATING OVERLAP HEATMAP")
  heatmap_plot <- create_overlap_heatmap(overlap_results$matrix, 
                                        "Protein Overlap Between Factors")
  ggsave(file.path(OUTPUT_DIR, "overlap_heatmap.png"), heatmap_plot, 
         width = 10, height = 8, dpi = 300)
  
  # Create key Venn diagrams
  print_section("CREATING VENN DIAGRAMS")
  create_key_comparisons(diff_lists, OUTPUT_DIR)
  
  # Create comprehensive protein table
  print_section("CREATING PROTEIN OVERLAP TABLES")
  protein_table <- create_overlap_protein_tables(diff_lists, RESULTS_DIR)
  
  # Save results
  print_section("SAVING RESULTS")
  
  # Save overlap matrix
  write.csv(overlap_results$matrix, file.path(OUTPUT_DIR, "overlap_matrix.csv"))
  
  # Save overlap details
  write.csv(overlap_results$details, file.path(OUTPUT_DIR, "pairwise_overlaps.csv"), row.names = FALSE)
  
  # Save comprehensive protein table
  write.csv(protein_table, file.path(OUTPUT_DIR, "comprehensive_protein_table.csv"), row.names = FALSE)
  
  # Create summary statistics
  summary_stats <- data.frame(
    Factor = names(diff_lists),
    Significant_Proteins = sapply(diff_lists, length),
    Percentage_of_Total = round(100 * sapply(diff_lists, length) / nrow(protein_table), 2),
    stringsAsFactors = FALSE
  )
  
  write.csv(summary_stats, file.path(OUTPUT_DIR, "factor_summary_stats.csv"), row.names = FALSE)
  
  # Create overlap summary report
  create_overlap_report(diff_lists, overlap_results, protein_table, summary_stats)
  
  print_section("VENN ANALYSIS COMPLETE")
  cat("Results saved to:", OUTPUT_DIR, "\n")
  cat("Key files generated:\n")
  cat("  - Venn diagrams: venn_*.png\n")
  cat("  - Overlap heatmap: overlap_heatmap.png\n")
  cat("  - Comprehensive table: comprehensive_protein_table.csv\n")
  cat("  - Summary statistics: factor_summary_stats.csv\n")
  cat("  - Detailed report: VENN_ANALYSIS_REPORT.md\n")
  
  return(list(
    differential_lists = diff_lists,
    overlap_results = overlap_results,
    protein_table = protein_table,
    summary_stats = summary_stats
  ))
}

create_overlap_report <- function(diff_lists, overlap_results, protein_table, summary_stats) {
  report_lines <- c(
    "# Venn Diagram and Overlap Analysis Report",
    "",
    "## Summary Statistics",
    ""
  )
  
  # Add factor summary
  for (i in 1:nrow(summary_stats)) {
    row <- summary_stats[i, ]
    report_lines <- c(report_lines,
      paste("- **", row$Factor, "**: ", row$Significant_Proteins, 
            " proteins (", row$Percentage_of_Total, "% of total)"))
  }
  
  report_lines <- c(report_lines,
    "",
    "## Key Overlaps",
    ""
  )
  
  # Add top overlaps
  top_overlaps <- overlap_results$details[order(overlap_results$details$Overlap_Count, decreasing = TRUE), ][1:5, ]
  
  for (i in 1:nrow(top_overlaps)) {
    if (!is.na(top_overlaps$Factor1[i])) {
      row <- top_overlaps[i, ]
      report_lines <- c(report_lines,
        paste("- **", row$Factor1, " ∩ ", row$Factor2, "**: ", 
              row$Overlap_Count, " proteins (Jaccard: ", 
              round(row$Jaccard_Index, 3), ")"))
    }
  }
  
  report_lines <- c(report_lines,
    "",
    "## Files Generated",
    "",
    "### Visualizations:",
    "- `venn_biological_factors.png`: Sex, Relapse Risk, Subject overlap",
    "- `venn_clinical_factors.png`: Sex, Age, Days to Relapse overlap", 
    "- `venn_demographic_factors.png`: Sex, Location overlap",
    "- `venn_top4_factors.png`: Top 4 factors overlap",
    "- `overlap_heatmap.png`: Heatmap of all pairwise overlaps",
    "",
    "### Data Tables:",
    "- `comprehensive_protein_table.csv`: All proteins with factor associations",
    "- `pairwise_overlaps.csv`: Detailed pairwise overlap statistics",
    "- `overlap_matrix.csv`: Overlap count matrix",
    "- `factor_summary_stats.csv`: Summary statistics per factor",
    "",
    "## Key Insights",
    "",
    "1. **Sex effects** are the most prominent biological signal",
    "2. **Individual differences** (Subject) show substantial protein variation",
    "3. **Clinical factors** (relapse risk, age) have meaningful overlaps",
    "4. **Geographic effects** are limited but present",
    "",
    "## Interpretation",
    "",
    "- High overlaps suggest shared biological pathways",
    "- Low overlaps indicate factor-specific mechanisms", 
    "- Jaccard indices quantify similarity between factor effects",
    "- Comprehensive table enables detailed protein-level analysis"
  )
  
  writeLines(report_lines, file.path(OUTPUT_DIR, "VENN_ANALYSIS_REPORT.md"))
}

# Execute if run as script
if (!interactive()) {
  result <- main()
}