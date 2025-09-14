#!/usr/bin/env Rscript

main_output_dir <- "output"
if (!dir.exists(main_output_dir)) dir.create(main_output_dir, recursive = TRUE)

subdirs <- c("enhanced_analysis", "overlap_analysis", "reports")
for (subdir in subdirs) {
  subdir_path <- file.path(main_output_dir, subdir)
  if (!dir.exists(subdir_path)) dir.create(subdir_path, recursive = TRUE)
}

source("scripts/complete_proteomics_analysis_enhanced.R")
if (dir.exists("output/proteomics_enhanced_analysis")) {
  file.copy("output/proteomics_enhanced_analysis", file.path(main_output_dir, "enhanced_analysis"), recursive = TRUE, overwrite = TRUE)
}

source("scripts/create_overlap_analysis.R")
if (dir.exists("output/overlap_analysis_results")) {
  file.copy("output/overlap_analysis_results", file.path(main_output_dir, "overlap_analysis"), recursive = TRUE, overwrite = TRUE)
}

report_lines <- c(
  "# Complete Proteomics Analysis Report",
  paste("**Date:** ", Sys.Date()),
  "## Results",
  "- Enhanced analysis: `enhanced_analysis/proteomics_enhanced_analysis/`",
  "- Overlap analysis: `overlap_analysis/overlap_analysis_results/`",
  "## Key Files",
  "- Corrected data: `enhanced_analysis/proteomics_enhanced_analysis/log2_intensity_corrected.csv`",
  "- Protein table: `overlap_analysis/overlap_analysis_results/comprehensive_protein_table.csv`",
  "- Overlaps: `overlap_analysis/overlap_analysis_results/pairwise_overlaps.csv`"
)

writeLines(report_lines, file.path(main_output_dir, "reports", "COMPLETE_PIPELINE_REPORT.md"))
cat("Pipeline completed:", main_output_dir, "\n")