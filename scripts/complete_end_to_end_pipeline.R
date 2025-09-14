#!/usr/bin/env Rscript
# Complete End-to-End Proteomics Analysis Pipeline
# Orchestrates all analysis steps in sequence

cat("=======================================================================\n")
cat("  COMPLETE END-TO-END PROTEOMICS ANALYSIS PIPELINE\n")
cat("=======================================================================\n\n")

# Create main output directory
main_output_dir <- "complete_pipeline_results"
if (!dir.exists(main_output_dir)) {
  dir.create(main_output_dir, recursive = TRUE)
}

# Create subdirectories
subdirs <- c("enhanced_analysis", "overlap_analysis", "reports")
for (subdir in subdirs) {
  subdir_path <- file.path(main_output_dir, subdir)
  if (!dir.exists(subdir_path)) {
    dir.create(subdir_path, recursive = TRUE)
  }
}

cat("Pipeline output directory:", main_output_dir, "\n\n")

# PHASE 1: Enhanced Proteomics Analysis
cat("PHASE 1: ENHANCED PROTEOMICS ANALYSIS\n")
cat("=====================================\n")

# Source and run the enhanced analysis
source("scripts/complete_proteomics_analysis_enhanced.R")

# Copy results to pipeline directory
if (dir.exists("proteomics_enhanced_analysis")) {
  file.copy("proteomics_enhanced_analysis", 
            file.path(main_output_dir, "enhanced_analysis"), 
            recursive = TRUE, overwrite = TRUE)
  cat("Enhanced analysis results copied to pipeline directory\n")
}

cat("\nPHASE 1 COMPLETED ✓\n\n")

# PHASE 2: Overlap Analysis
cat("PHASE 2: OVERLAP ANALYSIS\n")
cat("=========================\n")

# Source and run the overlap analysis
source("scripts/create_overlap_analysis.R")

# Copy results to pipeline directory
if (dir.exists("overlap_analysis_results")) {
  file.copy("overlap_analysis_results", 
            file.path(main_output_dir, "overlap_analysis"), 
            recursive = TRUE, overwrite = TRUE)
  cat("Overlap analysis results copied to pipeline directory\n")
}

cat("\nPHASE 2 COMPLETED ✓\n\n")

# PHASE 3: Generate Comprehensive Report
cat("PHASE 3: COMPREHENSIVE REPORTING\n")
cat("=================================\n")

# Create comprehensive summary report
report_lines <- c(
  "# Complete End-to-End Proteomics Analysis Report",
  "",
  paste("**Analysis Date:** ", Sys.Date()),
  "",
  "## Pipeline Overview",
  "",
  "This report summarizes the complete end-to-end proteomics analysis pipeline that includes:",
  "",
  "1. **Enhanced Proteomics Analysis** - Batch correction and comprehensive annotation analysis",
  "2. **Overlap Analysis** - Protein overlap analysis across different factors",
  "3. **Integrated Results** - Combined insights and recommendations",
  "",
  "## Key Results Summary",
  "",
  "### Enhanced Analysis Results",
  "- Complete batch correction with comprehensive annotation analysis",
  "- Results available in: `enhanced_analysis/proteomics_enhanced_analysis/`",
  "- Main report: `enhanced_analysis/proteomics_enhanced_analysis/ENHANCED_ANALYSIS_REPORT.md`",
  "",
  "### Overlap Analysis Results", 
  "- Comprehensive protein overlap analysis across factors",
  "- Results available in: `overlap_analysis/overlap_analysis_results/`",
  "- Main report: `overlap_analysis/overlap_analysis_results/OVERLAP_ANALYSIS_REPORT.md`",
  "",
  "## Directory Structure",
  "",
  "```",
  "complete_pipeline_results/",
  "├── enhanced_analysis/",
  "│   └── proteomics_enhanced_analysis/",
  "│       ├── ENHANCED_ANALYSIS_REPORT.md",
  "│       ├── log2_intensity_corrected.csv",
  "│       ├── annotation_effects_comparison.csv",
  "│       └── ...",
  "├── overlap_analysis/",
  "│   └── overlap_analysis_results/",
  "│       ├── OVERLAP_ANALYSIS_REPORT.md",
  "│       ├── comprehensive_protein_table.csv",
  "│       ├── pairwise_overlaps.csv",
  "│       └── ...",
  "└── reports/",
  "    └── COMPLETE_PIPELINE_REPORT.md (this file)",
  "```",
  "",
  "## Quick Access to Key Files",
  "",
  "### Main Data Files",
  "- **Batch-corrected data**: `enhanced_analysis/proteomics_enhanced_analysis/log2_intensity_corrected.csv`",
  "- **Sample annotations**: `enhanced_analysis/proteomics_enhanced_analysis/sample_annotation.csv`",
  "- **Protein annotations**: `enhanced_analysis/proteomics_enhanced_analysis/protein_annotations.csv`",
  "",
  "### Analysis Results",
  "- **Before/after comparison**: `enhanced_analysis/proteomics_enhanced_analysis/annotation_effects_comparison.csv`",
  "- **Comprehensive protein table**: `overlap_analysis/overlap_analysis_results/comprehensive_protein_table.csv`",
  "- **Pairwise overlaps**: `overlap_analysis/overlap_analysis_results/pairwise_overlaps.csv`",
  "",
  "### Visualizations",
  "- **Factor sizes**: `overlap_analysis/overlap_analysis_results/factor_sizes.png`",
  "- **Overlap heatmap**: `overlap_analysis/overlap_analysis_results/overlap_heatmap.png`",
  "- **Jaccard heatmap**: `overlap_analysis/overlap_analysis_results/jaccard_heatmap.png`",
  "",
  "## Next Steps",
  "",
  "1. **Review main reports** in each analysis directory",
  "2. **Examine key data files** for downstream analysis",
  "3. **Use batch-corrected data** for further statistical analysis",
  "4. **Investigate protein overlaps** for biological insights",
  "5. **Validate findings** in independent datasets",
  "",
  "## Technical Notes",
  "",
  "- **Batch correction method**: Median centering (DESeq2-MultiBatch adapted)",
  "- **Statistical threshold**: FDR < 0.05 (Benjamini-Hochberg correction)",
  "- **Overlap analysis**: Jaccard indices and absolute counts",
  "- **Missing data handling**: Robust methods with minimum validity requirements",
  "",
  "---",
  "",
  "*This report was generated by the Complete End-to-End Proteomics Analysis Pipeline*",
  paste("*Generated on:", Sys.time(), "*")
)

# Write comprehensive report
writeLines(report_lines, file.path(main_output_dir, "reports", "COMPLETE_PIPELINE_REPORT.md"))

cat("Comprehensive report generated\n")
cat("\nPHASE 3 COMPLETED ✓\n\n")

# FINAL SUMMARY
cat("=======================================================================\n")
cat("  COMPLETE PIPELINE FINISHED SUCCESSFULLY! 🎉\n")
cat("=======================================================================\n\n")

cat("All results saved to:", main_output_dir, "\n")
cat("Main report:", file.path(main_output_dir, "reports", "COMPLETE_PIPELINE_REPORT.md"), "\n\n")

cat("Key outputs:\n")
cat("- Enhanced analysis: enhanced_analysis/proteomics_enhanced_analysis/\n")
cat("- Overlap analysis: overlap_analysis/overlap_analysis_results/\n")
cat("- Comprehensive report: reports/COMPLETE_PIPELINE_REPORT.md\n\n")

cat("Pipeline completed successfully! Ready for downstream analysis.\n")

# Return success status
invisible(0)