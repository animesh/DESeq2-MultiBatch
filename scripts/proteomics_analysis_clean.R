#!/usr/bin/env Rscript
suppressMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(DESeq2)
})

OUTPUT_DIR <- "output/proteomics_enhanced_analysis"

create_output_dir <- function() {
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  cat("Output directory:", OUTPUT_DIR, "\n")
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

print_section <- function(title) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat(paste0("  ", title), "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
}

read_proteingroups_maxquant <- function(file_path) {
  cat("Reading proteinGroups from:", file_path, "\n")
  header <- fread(file_path, nrows = 0, header = TRUE)
  lfq_cols <- grep("^LFQ intensity ", colnames(header), value = TRUE)
  if (length(lfq_cols) == 0) stop("No LFQ intensity columns found")
  
  cols_to_read <- c("Protein IDs", "Gene names", "Protein names", lfq_cols)
  cols_to_read <- cols_to_read[cols_to_read %in% colnames(header)]
  
  data <- fread(file_path, select = cols_to_read, showProgress = FALSE)
  data <- data[!is.na(`Protein IDs`) & `Protein IDs` != "" & !grepl("^REV_|^CON_", `Protein IDs`)]
  
  sample_names <- gsub("^LFQ intensity ", "", lfq_cols)
  setnames(data, old = lfq_cols, new = sample_names)
  
  cat("Found", nrow(data), "proteins with", length(sample_names), "samples\n")
  return(data)
}

read_groups_maxquant <- function(file_path) {
  groups <- fread(file_path)
  if ("Condition" %in% colnames(groups)) groups[, Treatment := Condition]
  else if ("Bio" %in% colnames(groups)) groups[, Treatment := Bio]
  groups[, Name := gsub("\\.raw$|\\.mzML$", "", Name)]
  return(groups)
}

combine_batches_maxquant <- function(protein_b1, protein_b2, groups_b1, groups_b2) {
  groups_b1[, Batch := "B1"]
  groups_b2[, Batch := "B2"]
  combined_groups <- rbind(groups_b1, groups_b2, fill = TRUE)
  
  common_proteins <- intersect(protein_b1$`Protein IDs`, protein_b2$`Protein IDs`)
  cat("Found", length(common_proteins), "common proteins\n")
  if (length(common_proteins) == 0) stop("No common proteins found")
  
  protein_b1 <- protein_b1[`Protein IDs` %in% common_proteins]
  protein_b2 <- protein_b2[`Protein IDs` %in% common_proteins]
  
  sample_cols_b1 <- setdiff(colnames(protein_b1), c("Protein IDs", "Gene names", "Protein names"))
  sample_cols_b2 <- setdiff(colnames(protein_b2), c("Protein IDs", "Gene names", "Protein names"))
  
  intensity_b1 <- as.matrix(protein_b1[, ..sample_cols_b1])
  intensity_b2 <- as.matrix(protein_b2[, ..sample_cols_b2])
  rownames(intensity_b1) <- rownames(intensity_b2) <- sort(common_proteins)
  
  combined_intensity <- cbind(intensity_b1[sort(common_proteins), ], 
                             intensity_b2[sort(common_proteins), ])
  
  all_samples <- c(sample_cols_b1, sample_cols_b2)
  combined_annotation <- combined_groups[Name %in% all_samples]
  sample_order <- match(colnames(combined_intensity), combined_annotation$Name)
  combined_annotation <- combined_annotation[sample_order]
  
  list(intensity = combined_intensity, annotation = combined_annotation,
       protein_info = protein_b1[`Protein IDs` %in% rownames(combined_intensity), 
                                .(`Protein IDs`, `Gene names`, `Protein names`)])
}

preprocess_intensity_data <- function(intensity_matrix, min_valid_ratio = 0.5) {
  intensity_matrix[intensity_matrix == 0] <- NA
  log_intensity <- log2(intensity_matrix)
  min_valid <- ceiling(ncol(log_intensity) * min_valid_ratio)
  valid_proteins <- rowSums(!is.na(log_intensity)) >= min_valid
  log_intensity_filtered <- log_intensity[valid_proteins, ]
  cat("Retained", nrow(log_intensity_filtered), "proteins\n")
  return(log_intensity_filtered)
}

apply_batch_correction_deseq2 <- function(log_intensity, annotation) {
  corrected_intensity <- log_intensity
  batch_effects <- numeric(nrow(log_intensity))
  names(batch_effects) <- rownames(log_intensity)
  
  b1_samples <- which(annotation$Batch == "B1")
  b2_samples <- which(annotation$Batch == "B2")
  
  for (i in 1:nrow(log_intensity)) {
    b1_median <- median(log_intensity[i, b1_samples], na.rm = TRUE)
    b2_median <- median(log_intensity[i, b2_samples], na.rm = TRUE)
    
    if (!is.na(b1_median) && !is.na(b2_median)) {
      batch_effect <- b2_median - b1_median
      batch_effects[i] <- batch_effect
      corrected_intensity[i, b1_samples] <- corrected_intensity[i, b1_samples] + batch_effect/2
      corrected_intensity[i, b2_samples] <- corrected_intensity[i, b2_samples] - batch_effect/2
    }
  }
  
  cat("Mean absolute batch effect:", round(mean(abs(batch_effects), na.rm = TRUE), 3), "log2 units\n")
  list(corrected_intensity = corrected_intensity, batch_effects = batch_effects)
}

prepare_annotation_data <- function(annotation) {
  clean_annotation <- data.frame(Sample_ID = annotation$Name, Batch = annotation$Batch, stringsAsFactors = FALSE)
  
  if ("Gender" %in% colnames(annotation)) {
    clean_annotation$Sex <- ifelse(tolower(annotation$Gender) %in% c("male", "m"), "Male", "Female")
  }
  
  treatment_cols <- c("Treatment", "Bio", "Condition")
  for (col in treatment_cols) {
    if (col %in% colnames(annotation)) {
      clean_annotation$Treatment <- as.character(annotation[[col]])
      break
    }
  }
  
  if ("Age" %in% colnames(annotation)) {
    clean_annotation$Age <- as.numeric(annotation$Age)
    clean_annotation$Age_Group <- cut(clean_annotation$Age, breaks = c(0, 60, 70, 100), 
                                     labels = c("Young", "Middle", "Old"), include.lowest = TRUE)
  }
  
  if ("Location" %in% colnames(annotation)) {
    clean_annotation$Location <- as.character(annotation$Location)
  }
  
  if ("Days until relaps" %in% colnames(annotation)) {
    days_raw <- annotation$`Days until relaps`
    days_numeric <- as.numeric(gsub(" days", "", days_raw))
    clean_annotation$Days_Until_Relapse <- days_numeric
    if (any(!is.na(days_numeric))) {
      clean_annotation$Relapse_Risk <- cut(days_numeric, breaks = c(0, 500, 1000, Inf),
                                          labels = c("High_Risk", "Medium_Risk", "Low_Risk"), include.lowest = TRUE)
    }
  }
  
  if ("Subject" %in% colnames(annotation)) {
    clean_annotation$Subject <- as.character(annotation$Subject)
  }
  
  if ("Cell type" %in% colnames(annotation)) {
    clean_annotation$Cell_Type <- as.character(annotation$`Cell type`)
  }
  
  return(clean_annotation)
}

analyze_all_annotations <- function(intensity_matrix, annotation, data_type = "raw") {
  clean_annotation <- prepare_annotation_data(annotation)
  
  categorical_vars <- c()
  numerical_vars <- c()
  
  for (col in colnames(clean_annotation)) {
    if (col %in% c("Sample_ID")) next
    if (is.numeric(clean_annotation[[col]])) {
      numerical_vars <- c(numerical_vars, col)
    } else {
      var_table <- table(clean_annotation[[col]], useNA = "no")
      if (length(var_table) >= 2 && min(var_table) >= 3) {
        categorical_vars <- c(categorical_vars, col)
      }
    }
  }
  
  categorical_results <- list()
  for (var in categorical_vars) {
    cat("Analyzing", var, "...\n")
    result <- analyze_categorical_variable(intensity_matrix, clean_annotation, var)
    if (!is.null(result)) categorical_results[[var]] <- result
  }
  
  numerical_results <- list()
  for (var in numerical_vars) {
    cat("Analyzing", var, "...\n")
    result <- analyze_numerical_variable(intensity_matrix, clean_annotation, var)
    if (!is.null(result)) numerical_results[[var]] <- result
  }
  
  list(categorical = categorical_results, numerical = numerical_results, annotation = clean_annotation)
}

run_deseq2_analysis <- function(intensity_clean, coldata, design_formula, variable) {
  pseudo_counts <- round(2^intensity_clean)
  pseudo_counts[is.na(pseudo_counts)] <- 0
  pseudo_counts[pseudo_counts < 1] <- 1
  
  protein_results <- data.frame(
    Protein_ID = rownames(intensity_clean), Effect_Size = 0, P_value = 1, Adj_P_value = 1,
    stringsAsFactors = FALSE)
  
  tryCatch({
    dds <- DESeqDataSetFromMatrix(countData = pseudo_counts, colData = coldata, design = design_formula)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    dds <- DESeq(dds, quiet = TRUE)
    
    var_levels <- levels(coldata$Variable)
    res <- tryCatch({
      contrast_name <- paste0("Variable_", var_levels[2], "_vs_", var_levels[1])
      results(dds, name = contrast_name, independentFiltering = FALSE)
    }, error = function(e) {
      results(dds, contrast = c("Variable", var_levels[2], var_levels[1]), independentFiltering = FALSE)
    })
    
    protein_results <- data.frame(
      Protein_ID = rownames(res), Effect_Size = res$log2FoldChange,
      P_value = res$pvalue, Adj_P_value = res$padj, stringsAsFactors = FALSE)
    
    all_proteins <- rownames(intensity_clean)
    missing_proteins <- setdiff(all_proteins, protein_results$Protein_ID)
    if (length(missing_proteins) > 0) {
      missing_results <- data.frame(
        Protein_ID = missing_proteins, Effect_Size = NA, P_value = NA, Adj_P_value = NA,
        stringsAsFactors = FALSE)
      protein_results <- rbind(protein_results, missing_results)
    }
    protein_results <- protein_results[match(all_proteins, protein_results$Protein_ID), ]
    
  }, error = function(e) {
    cat("DESeq2 analysis failed for", variable, "\n")
    return(protein_results)
  })
  
  valid_pvals <- !is.na(protein_results$P_value)
  if (sum(valid_pvals) > 0 && any(is.na(protein_results$Adj_P_value[valid_pvals]))) {
    protein_results$Adj_P_value[valid_pvals] <- p.adjust(protein_results$P_value[valid_pvals], method = "BH")
  }
  
  significant_proteins <- sum(protein_results$Adj_P_value < 0.05, na.rm = TRUE)
  cat("  -", variable, ":", significant_proteins, "significant proteins\n")
  return(protein_results)
}

analyze_categorical_variable <- function(intensity_matrix, annotation, variable) {
  var_values <- annotation[[variable]]
  valid_samples <- !is.na(var_values)
  if (sum(valid_samples) < 10) return(NULL)
  
  intensity_clean <- intensity_matrix[, valid_samples]
  var_clean <- var_values[valid_samples]
  
  coldata <- data.frame(Variable = factor(var_clean), Batch = factor(annotation$Batch[valid_samples]))
  
  if ("Gender" %in% colnames(annotation)) {
    sex_values <- annotation$Gender[valid_samples]
    if (length(unique(sex_values[!is.na(sex_values)])) > 1) {
      coldata$Sex <- factor(sex_values)
    }
  }
  
  if (length(levels(coldata$Variable)) < 2) {
    return(data.frame(Protein_ID = rownames(intensity_clean), Effect_Size = 0, 
                     P_value = 1, Adj_P_value = 1, stringsAsFactors = FALSE))
  }
  
  design_formula <- if ("Sex" %in% colnames(coldata)) ~ Batch + Sex + Variable else ~ Batch + Variable
  return(run_deseq2_analysis(intensity_clean, coldata, design_formula, variable))
}

analyze_numerical_variable <- function(intensity_matrix, annotation, variable) {
  var_values <- annotation[[variable]]
  valid_samples <- !is.na(var_values)
  if (sum(valid_samples) < 10) return(NULL)
  
  intensity_clean <- intensity_matrix[, valid_samples]
  var_clean <- var_values[valid_samples]
  
  coldata <- data.frame(Variable = as.numeric(var_clean), Batch = factor(annotation$Batch[valid_samples]))
  
  if ("Gender" %in% colnames(annotation)) {
    sex_values <- annotation$Gender[valid_samples]
    if (length(unique(sex_values[!is.na(sex_values)])) > 1) {
      coldata$Sex <- factor(sex_values)
    }
  }
  
  design_formula <- if ("Sex" %in% colnames(coldata)) ~ Batch + Sex + Variable else ~ Batch + Variable
  return(run_deseq2_analysis(intensity_clean, coldata, design_formula, variable))
}

compare_annotation_effects <- function(raw_results, corrected_results) {
  all_vars <- unique(c(names(raw_results$categorical), names(raw_results$numerical),
                      names(corrected_results$categorical), names(corrected_results$numerical)))
  
  comparison_data <- data.frame(Variable = character(), Type = character(), 
                               Significant_Raw = integer(), Significant_Corrected = integer(),
                               Change = integer(), Percent_Change = numeric(), stringsAsFactors = FALSE)
  
  for (var in all_vars) {
    var_type <- if (var %in% names(raw_results$categorical)) "Categorical" else "Numerical"
    
    raw_sig <- if (var_type == "Categorical" && var %in% names(raw_results$categorical)) {
      sum(raw_results$categorical[[var]]$Adj_P_value < 0.05, na.rm = TRUE)
    } else if (var_type == "Numerical" && var %in% names(raw_results$numerical)) {
      sum(raw_results$numerical[[var]]$Adj_P_value < 0.05, na.rm = TRUE)
    } else 0
    
    corr_sig <- if (var_type == "Categorical" && var %in% names(corrected_results$categorical)) {
      sum(corrected_results$categorical[[var]]$Adj_P_value < 0.05, na.rm = TRUE)
    } else if (var_type == "Numerical" && var %in% names(corrected_results$numerical)) {
      sum(corrected_results$numerical[[var]]$Adj_P_value < 0.05, na.rm = TRUE)
    } else 0
    
    change <- corr_sig - raw_sig
    percent_change <- if (raw_sig > 0) (change / raw_sig) * 100 else NA
    
    comparison_data <- rbind(comparison_data, data.frame(
      Variable = var, Type = var_type, Significant_Raw = raw_sig, Significant_Corrected = corr_sig,
      Change = change, Percent_Change = percent_change, stringsAsFactors = FALSE))
  }
  
  return(comparison_data)
}

save_enhanced_results <- function(combined_data, raw_intensity, corrected_intensity, 
                                 batch_correction_result, raw_results, corrected_results, 
                                 comparison_summary) {
  fwrite(combined_data$annotation, file.path(OUTPUT_DIR, "sample_annotation.csv"))
  fwrite(as.data.frame(raw_intensity), file.path(OUTPUT_DIR, "log2_intensity_raw.csv"), row.names = TRUE)
  fwrite(as.data.frame(corrected_intensity), file.path(OUTPUT_DIR, "log2_intensity_corrected.csv"), row.names = TRUE)
  fwrite(combined_data$protein_info, file.path(OUTPUT_DIR, "protein_annotations.csv"))
  
  batch_effects_df <- data.frame(Protein_IDs = names(batch_correction_result$batch_effects),
                                Batch_Effect = batch_correction_result$batch_effects, stringsAsFactors = FALSE)
  fwrite(batch_effects_df, file.path(OUTPUT_DIR, "batch_effects.csv"))
  
  raw_dir <- file.path(OUTPUT_DIR, "raw_data_results")
  corr_dir <- file.path(OUTPUT_DIR, "corrected_data_results")
  if (!dir.exists(raw_dir)) dir.create(raw_dir)
  if (!dir.exists(corr_dir)) dir.create(corr_dir)
  
  for (var in names(raw_results$categorical)) {
    fwrite(raw_results$categorical[[var]], file.path(raw_dir, paste0("raw_", var, "_categorical.csv")))
  }
  for (var in names(raw_results$numerical)) {
    fwrite(raw_results$numerical[[var]], file.path(raw_dir, paste0("raw_", var, "_numerical.csv")))
  }
  for (var in names(corrected_results$categorical)) {
    fwrite(corrected_results$categorical[[var]], file.path(corr_dir, paste0("corrected_", var, "_categorical.csv")))
  }
  for (var in names(corrected_results$numerical)) {
    fwrite(corrected_results$numerical[[var]], file.path(corr_dir, paste0("corrected_", var, "_numerical.csv")))
  }
  
  fwrite(comparison_summary, file.path(OUTPUT_DIR, "annotation_effects_comparison.csv"))
  cat("Enhanced results saved to:", OUTPUT_DIR, "\n")
}

create_enhanced_summary_report <- function(combined_data, raw_intensity, corrected_intensity, 
                                          batch_correction_result, raw_results, corrected_results, 
                                          comparison_summary) {
  report_lines <- c(
    "# Enhanced Proteomics Analysis Report", "",
    "## Data Summary",
    paste("- Total proteins analyzed:", nrow(raw_intensity)),
    paste("- Total samples:", ncol(raw_intensity)),
    paste("- Batch 1 samples:", sum(combined_data$annotation$Batch == "B1")),
    paste("- Batch 2 samples:", sum(combined_data$annotation$Batch == "B2")), "",
    "## Batch Correction Results",
    paste("- Mean absolute batch effect:", round(mean(abs(batch_correction_result$batch_effects), na.rm = TRUE), 3), "log2 units"),
    "- Correction method: DESeq2-inspired symmetric median centering",
    "- Statistical analysis: Full DESeq2 pipeline with negative binomial modeling", "",
    "## Comprehensive Annotation Analysis", ""
  )
  
  if (length(raw_results$categorical) > 0) {
    report_lines <- c(report_lines, "### Categorical Variables:")
    for (var in names(raw_results$categorical)) {
      raw_sig <- sum(raw_results$categorical[[var]]$Adj_P_value < 0.05, na.rm = TRUE)
      if (var %in% names(corrected_results$categorical)) {
        corr_sig <- sum(corrected_results$categorical[[var]]$Adj_P_value < 0.05, na.rm = TRUE)
        report_lines <- c(report_lines, paste("- **", var, "**:", raw_sig, "→", corr_sig, "significant proteins"))
      }
    }
    report_lines <- c(report_lines, "")
  }
  
  if (length(raw_results$numerical) > 0) {
    report_lines <- c(report_lines, "### Numerical Variables:")
    for (var in names(raw_results$numerical)) {
      raw_sig <- sum(raw_results$numerical[[var]]$Adj_P_value < 0.05, na.rm = TRUE)
      if (var %in% names(corrected_results$numerical)) {
        corr_sig <- sum(corrected_results$numerical[[var]]$Adj_P_value < 0.05, na.rm = TRUE)
        report_lines <- c(report_lines, paste("- **", var, "**:", raw_sig, "→", corr_sig, "significant proteins"))
      }
    }
  }
  
  writeLines(report_lines, file.path(OUTPUT_DIR, "ENHANCED_ANALYSIS_REPORT.md"))
  cat("Enhanced summary report saved to:", file.path(OUTPUT_DIR, "ENHANCED_ANALYSIS_REPORT.md"), "\n")
}

main <- function() {
  create_output_dir()
  
  print_section("ENHANCED PROTEOMICS ANALYSIS PIPELINE")
  print_section("STEP 1: READING DATA")
  
  protein_b1 <- read_proteingroups_maxquant("txtB1/proteinGroups.txt")
  protein_b2 <- read_proteingroups_maxquant("txtB2/proteinGroups.txt")
  groups_b1 <- read_groups_maxquant("txtB1/Groups.txt")
  groups_b2 <- read_groups_maxquant("txtB2/Groups.txt")
  
  print_section("STEP 2: COMBINING BATCHES")
  combined_data <- combine_batches_maxquant(protein_b1, protein_b2, groups_b1, groups_b2)
  
  print_section("STEP 3: PREPROCESSING DATA")
  raw_intensity <- preprocess_intensity_data(combined_data$intensity)
  
  print_section("STEP 4: BATCH CORRECTION")
  batch_correction_result <- apply_batch_correction_deseq2(raw_intensity, combined_data$annotation)
  corrected_intensity <- batch_correction_result$corrected_intensity
  
  print_section("STEP 5: COMPREHENSIVE ANNOTATION ANALYSIS - RAW DATA")
  raw_results <- analyze_all_annotations(raw_intensity, combined_data$annotation, "raw")
  
  print_section("STEP 6: COMPREHENSIVE ANNOTATION ANALYSIS - CORRECTED DATA")
  corrected_results <- analyze_all_annotations(corrected_intensity, combined_data$annotation, "corrected")
  
  print_section("STEP 7: COMPARING EFFECTS BEFORE AND AFTER CORRECTION")
  comparison_summary <- compare_annotation_effects(raw_results, corrected_results)
  cat("Comparison summary:\n")
  print(comparison_summary)
  
  print_section("STEP 8: SAVING ENHANCED RESULTS")
  save_enhanced_results(combined_data, raw_intensity, corrected_intensity, 
                       batch_correction_result, raw_results, corrected_results, comparison_summary)
  
  print_section("STEP 9: GENERATING ENHANCED SUMMARY REPORT")
  create_enhanced_summary_report(combined_data, raw_intensity, corrected_intensity, 
                                batch_correction_result, raw_results, corrected_results, comparison_summary)
  
  print_section("ENHANCED ANALYSIS COMPLETE")
  cat("All results saved to:", OUTPUT_DIR, "\n")
  cat("Check ENHANCED_ANALYSIS_REPORT.md for a complete summary\n")
  
  return(TRUE)
}

if (!interactive()) {
  main()
}