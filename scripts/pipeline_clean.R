#!/usr/bin/env Rscript
main_output_dir <- "output"
subdirs <- c("enhanced_analysis", "overlap_analysis", "reports")

for (subdir in subdirs) {
  subdir_path <- file.path(main_output_dir, subdir)
  if (!dir.exists(subdir_path)) dir.create(subdir_path, recursive = TRUE)
}

cat("Running enhanced proteomics analysis...\n")
source("scripts/proteomics_analysis_clean.R")
if (dir.exists("output/proteomics_enhanced_analysis")) {
  file.copy("output/proteomics_enhanced_analysis", file.path(main_output_dir, "enhanced_analysis"), recursive = TRUE, overwrite = TRUE)
}

cat("Running overlap analysis...\n")
source("scripts/overlap_analysis_clean.R")
if (dir.exists("output/overlap_analysis_results")) {
  file.copy("output/overlap_analysis_results", file.path(main_output_dir, "overlap_analysis"), recursive = TRUE, overwrite = TRUE)
}

report_lines <- c(
  "# Complete Proteomics Analysis Report",
  paste("**Date:**", Sys.Date()),
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