#!/usr/bin/env Rscript
# DESeq2-MultiBatch implementation script
# Implements the steps described in the repository README.md

suppressMessages({
  # helper to ensure package is available
  ensure_pkg <- function(pkg) {
    if (!suppressWarnings(requireNamespace(pkg, quietly = TRUE))) {
      if (!suppressWarnings(requireNamespace("BiocManager", quietly = TRUE))) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE)
    }
    invisible(TRUE)
  }
})

# Required packages
ensure_pkg("DESeq2")
library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)

# support optional dry-run flag: call as
# Rscript scripts/run_deseq2_multibatch.R <data_dir> <out_dir> [--dry-run]
dry_run <- FALSE
if (any(args == "--dry-run" | args == "-n")) {
  dry_run <- TRUE
  args <- args[args != "--dry-run" & args != "-n"]
}

data_dir <- if (length(args) >= 1) args[[1]] else "data"
out_dir <- if (length(args) >= 2) args[[2]] else "outputs"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

coldata_path <- file.path(data_dir, "coldata.csv")
counts_path <- file.path(data_dir, "counts.csv")

if (!file.exists(coldata_path) || !file.exists(counts_path)) {
  stop("Required data files not found in '", data_dir, "'. Expected 'coldata.csv' and 'counts.csv'.")
}

message("Reading data from: ", data_dir)
# read coldata: expect first column to be sample names -> use as rownames
coldata <- read.csv(coldata_path, row.names = 1, check.names = FALSE)

# convert character columns to factors (as in README)
coldata[] <- lapply(coldata, function(x) if (is.character(x) || is.factor(x)) as.factor(x) else x)

# read counts: rows = genes, cols = samples
counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)

# ensure samples order matches (robust matching)
cn <- colnames(counts)
rn <- rownames(coldata)

if (all(cn %in% rn)) {
  # exact match
  coldata <- coldata[cn, , drop = FALSE]
} else {
  message("Counts/sample names do not match exactly. Attempting normalization and smart mapping...")

  # 1) try stripping leading 'X' that sometimes appears when R made syntactic names
  cn_stripX <- sub("^X", "", cn)

  # 2) if column names are numeric (e.g. '1','2',...) try prefixing with 'Sample'
  cn_pref_sample <- ifelse(grepl('^[0-9]+$', cn_stripX), paste0('Sample', cn_stripX), cn_stripX)

  use_names <- NULL
  if (all(cn_stripX %in% rn)) {
    use_names <- cn_stripX
  } else if (all(cn_pref_sample %in% rn)) {
    use_names <- cn_pref_sample
  } else {
    # 3) try a relaxed matching by removing non-alphanumeric and lowercasing
    clean <- function(x) tolower(gsub('[^0-9A-Za-z]', '', x))
    rn_clean_map <- setNames(rn, clean(rn))
    cn_clean <- clean(cn)
    if (all(cn_clean %in% names(rn_clean_map))) {
      # map cleaned column-names to original rownames
      use_names <- rn_clean_map[cn_clean]
    }
  }

  if (is.null(use_names)) {
    stop("Sample names in counts and coldata do not match after attempted normalization.\n",
         "counts columns (first 10): ", paste(head(cn, 10), collapse = ", "), "\n",
         "coldata rownames (first 10): ", paste(head(rn, 10), collapse = ", "), "\n",
         "Please ensure sample IDs match (rownames in coldata and column names in counts).")
  }

  # reorder coldata to match counts; if we generated mapped names that differ from the
  # original column names, also rename the counts columns so they match rownames(coldata)
  coldata <- coldata[use_names, , drop = FALSE]
  if (!all(use_names == cn)) {
    colnames(counts) <- use_names
  }
}

# if dry-run requested, write an alignment marker and exit (useful for tests)
if (exists('dry_run') && isTRUE(dry_run)) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  writeLines(c('alignment: OK', paste('mapped_columns:', paste(colnames(counts), collapse = ','))), con = file.path(out_dir, 'alignment_ok.txt'))
  message('Dry-run: alignment successful. Wrote alignment marker to: ', file.path(out_dir, 'alignment_ok.txt'))
  quit(status = 0)
}

### SINGLE BATCH DESIGN
message("Running DESeq2 - single batch design")
design_single <- as.formula("~ Day + Batch + Genotype + Sex + Sex:Day + Treatment")
dds_single <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = coldata, design = design_single)
dds_single <- DESeq(dds_single)
saveRDS(dds_single, file = file.path(out_dir, "dds_single.rds"))

### DOUBLE BATCH DESIGN
message("Running DESeq2 - double batch design (interaction with Sex)")
design_double <- as.formula("~ Day + Batch + Genotype + Sex + Sex:Batch + Sex:Day + Treatment")
dds_double <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = coldata, design = design_double)
dds_double <- DESeq(dds_double)
saveRDS(dds_double, file = file.path(out_dir, "dds_double.rds"))

### Extract contrasts and compute scaling factors
message("Extracting contrasts and computing scaling factors")

# helper to compute scaling factor dataframe from a results object
compute_scaling_df <- function(res_df) {
  res_df$log2FoldChange[is.na(res_df$log2FoldChange)] <- 0
  scaling_A <- sqrt(2^(res_df$log2FoldChange))
  scaling_B <- 1 / scaling_A
  out <- res_df
  out$Scaling_Factor_Batch_A <- scaling_A
  out$Scaling_Factor_Batch_B <- scaling_B
  return(out)
}

# single: try to use coefficient name if present else use contrast
rn_single <- resultsNames(dds_single)
if ("Batch_b_vs_a" %in% rn_single) {
  batch_res <- results(dds_single, name = "Batch_b_vs_a", independentFiltering = FALSE)
} else {
  batch_res <- results(dds_single, contrast = c("Batch", "b", "a"), independentFiltering = FALSE)
}
batch_df <- as.data.frame(batch_res)
batch_df <- compute_scaling_df(batch_df)
write.csv(batch_df, file = file.path(out_dir, "batch_single_results.csv"), row.names = TRUE)

# double: get base batch and the sex-specific (interaction) batch
rn_double <- resultsNames(dds_double)
if ("Batch_b_vs_a" %in% rn_double) {
  batch_xx_res <- results(dds_double, name = "Batch_b_vs_a", independentFiltering = FALSE)
} else {
  batch_xx_res <- results(dds_double, contrast = c("Batch", "b", "a"), independentFiltering = FALSE)
}
batch_xx_df <- as.data.frame(batch_xx_res)

# attempt to compute the sex-specific batch contrast following README
batch_xy_df <- NULL
if ("Batchb.Sexxy" %in% rn_double) {
  batch_xy_res <- results(dds_double, contrast = list(c("Batch_b_vs_a", "Batchb.Sexxy")), independentFiltering = FALSE)
  batch_xy_df <- as.data.frame(batch_xy_res)
} else {
  # warn and fallback: use the base batch for both sexes (user should inspect resultsNames)
  warning("Could not find Sex-specific batch interaction coefficient name in resultsNames(dds_double).\n",
          "Falling back to the base batch effect for both sexes. Check resultsNames(dds_double) for available coefficients.")
  batch_xy_df <- batch_xx_df
}

batch_xx_df <- compute_scaling_df(batch_xx_df)
batch_xy_df <- compute_scaling_df(batch_xy_df)
write.csv(batch_xx_df, file = file.path(out_dir, "batch_double_xx_results.csv"), row.names = TRUE)
write.csv(batch_xy_df, file = file.path(out_dir, "batch_double_xy_results.csv"), row.names = TRUE)

### Apply corrections to normalized counts
message("Computing normalized counts and applying scaling (single and double)")
norm_single <- as.data.frame(counts(dds_single, normalized = TRUE))
scaled_single <- norm_single

for (sample_id in colnames(norm_single)) {
  Batch <- as.character(coldata[sample_id, "Batch"])
  if (is.na(Batch)) stop("Batch is NA for sample: ", sample_id)
  if (Batch == "a") {
    scaling_factors <- batch_df$Scaling_Factor_Batch_A
  } else if (Batch == "b") {
    scaling_factors <- batch_df$Scaling_Factor_Batch_B
  } else {
    stop("Unexpected Batch value for sample: ", sample_id)
  }
  # multiply per-gene scaling factors
  scaled_single[, sample_id] <- norm_single[, sample_id] * scaling_factors
}

norm_double <- as.data.frame(counts(dds_double, normalized = TRUE))
scaled_double <- norm_double
for (sample_id in colnames(norm_double)) {
  Sex <- as.character(coldata[sample_id, "Sex"])
  Batch <- as.character(coldata[sample_id, "Batch"])
  if (is.na(Sex) || is.na(Batch)) stop("Sex or Batch is NA for sample: ", sample_id)

  if (Sex == "xx" && Batch == "a") {
    scaling_factors <- batch_xx_df$Scaling_Factor_Batch_A
  } else if (Sex == "xx" && Batch == "b") {
    scaling_factors <- batch_xx_df$Scaling_Factor_Batch_B
  } else if (Sex == "xy" && Batch == "a") {
    scaling_factors <- batch_xy_df$Scaling_Factor_Batch_A
  } else if (Sex == "xy" && Batch == "b") {
    scaling_factors <- batch_xy_df$Scaling_Factor_Batch_B
  } else {
    stop("Unexpected Sex/Batch combination for sample: ", sample_id)
  }

  scaled_double[, sample_id] <- norm_double[, sample_id] * scaling_factors
}

## write outputs
write.csv(norm_single, file = file.path(out_dir, "normalized_counts_single.csv"), row.names = TRUE)
write.csv(scaled_single, file = file.path(out_dir, "scaled_counts_single.csv"), row.names = TRUE)

write.csv(norm_double, file = file.path(out_dir, "normalized_counts_double.csv"), row.names = TRUE)
write.csv(scaled_double, file = file.path(out_dir, "scaled_counts_double.csv"), row.names = TRUE)

message("Finished. Outputs written to: ", normalizePath(out_dir))