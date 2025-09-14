# The steps here are quite simple : quantify the effect of batch on each gene
# and the effect of both the batch and the interaction term between batch and sex on each gene.

# The first will be compensated for in XX samples
# while the second will be compensated for on XY samples
# (since reference level is XX, interaction term is XY-related)
#install.packages("BiocManager")
#BiocManager::install(c("DESeq2", "HTSFilter", "dplyr", "ggplot2"))
library(DESeq2)
library(HTSFilter)
library(dplyr)
library(ggplot2)

# Define each part of the path as a separate component
drive <- "C:"
users <- "Users"
user_name <- "julie"
onedrive <- "OneDrive - Université Laval"
folder <- "Stage recherche/WGCNA_directory/DESeq2_batch_correction"

# Use file.path() to construct the full path
orig_full_path <- "/home/ash022/DESeq2-MultiBatch/data/"

# Candidate data directories to try (project-local and known absolute paths)
candidate_paths <- c(
  orig_full_path,
  file.path(getwd(), "data"),
  "/mnt/z/Download/DESeq2-MultiBatch/data",
  file.path(dirname(getwd()), "data")
)

# Pick the first existing candidate
existing <- candidate_paths[dir.exists(candidate_paths) | file.exists(candidate_paths)]
if (length(existing) == 0) {
  stop(sprintf("data directory not found. Tried: %s\nPlease create a 'data' directory in the repo root or adjust the path in scripts/orig.r",
               paste(candidate_paths, collapse = ", ")))
}

full_path <- existing[1]

# Set the working directory for knitting
knitr::opts_knit$set(root.dir = full_path)

# Set the working directory
setwd(full_path)
getwd()

# Load count data
counts <- read.csv("Z:/Download/DESeq2-MultiBatch/data/counts.csv", sep = ",", row.names = 1, check.names = FALSE)
plot(counts[,1], counts[,2]) # Quick check to see if data looks ok
# Load metadata
colData <- read.csv("Z:/Download/DESeq2-MultiBatch/data/metadata.csv")
summary(colData)
# Convert the condition columns to factors
# Ensure sample IDs are rownames so later indexing by sample name works
if (!"Sample_ID" %in% colnames(colData)) {
  stop('metadata.csv must contain a Sample_ID column')
}
rownames(colData) <- colData$Sample_ID

# Normalize simple casing issues so comparisons later match expected levels
colData$Sex <- toupper(as.character(colData$Sex))
colData$Batch <- tolower(as.character(colData$Batch))

colData$Treatment <- as.factor(colData$Treatment)
colData$Cultivar <- as.factor(colData$Cultivar)
colData$Day <- as.factor(colData$Day)
colData$Sex <- as.factor(colData$Sex)
colData$Batch <- as.factor(colData$Batch)
colData$Sample_ID <- as.factor(colData$Sample_ID)

# Check the result
head(colData)

# DESIGN
# Construct the design
# For batch correction purposes, the design was simplified while remaining accurate (only the interaction between treatment and cultivar was removed)
design <- ~ Treatment + Cultivar + Sex + Day + Sex:Day + Sex:Batch + Batch

# To construct the object, use counts as matrix, convert first
count_matrix <- as.matrix(counts)

# If count columns are numeric indices (e.g. "1","2",...), map by position to metadata Sample_IDs
if (all(grepl('^\\d+$', colnames(count_matrix)))) {
  if (ncol(count_matrix) == nrow(colData)) {
    cat('Detected numeric column names in counts; mapping by position to metadata Sample_IDs for initial DESeqDataSet.\n')
    colnames(count_matrix) <- rownames(colData)
  } else {
    stop('Counts have numeric column names and number of columns does not match metadata rows; cannot map by position for initial dataset.')
  }
}

# construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = design)

# Define reference value for each factor (needed for simplicity, helps calling for contrasts)
dds$Treatment <- relevel(dds$Treatment, ref = "ctrl")
if ("DK" %in% levels(dds$Cultivar)) {
  dds$Cultivar <- relevel(dds$Cultivar, ref = "DK")
} else {
  # Fallback to the first available level if 'DK' is not present in the metadata
  dds$Cultivar <- relevel(dds$Cultivar, ref = levels(dds$Cultivar)[1])
}
if ("XX" %in% levels(dds$Sex)) {
  dds$Sex <- relevel(dds$Sex, ref = "XX")
} else {
  dds$Sex <- relevel(dds$Sex, ref = levels(dds$Sex)[1])
}
if ("a" %in% levels(dds$Batch)) {
  dds$Batch <- relevel(dds$Batch, ref = "a")
} else {
  dds$Batch <- relevel(dds$Batch, ref = levels(dds$Batch)[1])
}
if ("0" %in% levels(dds$Day)) {
  dds$Day <- relevel(dds$Day, ref = "0")
} else {
  dds$Day <- relevel(dds$Day, ref = levels(dds$Day)[1])
}

# Run DESeq
dds <-DESeq(dds)

# Remove low count genes with HTSFilter
dds_filtered <- HTSFilter(dds, s.min = 15, s.max = 50, s.len=25, plot = TRUE)$filteredData

# Since there are no reason to believe any batch (winter v. spring) is more representative of the biological nature of gene expression, the batches will be "averaged"
# For example, if a gene count is around 1000 for XX plants in batch A, and around 2000 for XX plants in batch B, the LFC of the "batch" factor will be 1
# Instead of aiming for the average expression count, which would be 1500 (batch A * 1.5, batch B / 1.333), requiring a different factor for each (1.5 v. 1.333),
# it seems more intuitive and justified to scale both batches by a equivalent factor
# This adjusts counts symmetrically in log space, maintaining consistency in how adjustments are applied across all genes
# Both batches are brought closer to each other in a proportional manner on the logarithmic scale
# All normalized counts for XX samples (all days) from batch A will be timed by √(batch B/batch A) = √(2) = 1.4142,
# while all normalized counts for XX samples (all days) from batch B will be divided by 1.4142, or timed by 0.7071
# For any LFC = X, factor for batch A is √(2^x) and factor for batch B is 1/(factor for batch A)

# Get the batch effect and the interaction term batch effect (robust selection)
rn <- resultsNames(dds_filtered)
cat("resultsNames(dds_filtered):", paste(rn, collapse = ", "), "\n")

# Try to find the batch main effect name (prefer exact 'Batch_b_vs_a')
batch_name <- if ("Batch_b_vs_a" %in% rn) {
  "Batch_b_vs_a"
} else {
  bn <- rn[grep("Batch.*b.*vs.*a|Batch_b_vs_a|Batchb.*vs.*a|Batch_b_vs_a", rn, ignore.case = TRUE)]
  if (length(bn) > 0) bn[1] else rn[grep("Batch", rn, ignore.case = TRUE)][1]
}
cat("Selected batch_name:", batch_name, "\n")

# Try to find an interaction name involving Sex and Batch
## Prefer interaction-like names that mention Batch but exclude simple main-effect "_vs_" patterns
interaction_candidates <- rn[grepl("batch", rn, ignore.case = TRUE) & grepl("sex|sexxy|sex_xy|sex.xy|sex:", rn, ignore.case = TRUE) & !grepl("_vs_", rn, ignore.case = TRUE)]
if (length(interaction_candidates) == 0) {
  ## fallback: any sex-related name that is not a main-effect _vs_ entry
  interaction_candidates <- rn[grepl("sex", rn, ignore.case = TRUE) & !grepl("_vs_", rn, ignore.case = TRUE)]
}
interaction_name <- if (length(interaction_candidates) > 0) interaction_candidates[1] else NA
cat("Selected interaction_name:", interaction_name, "\n")

res_batch_XX <- tryCatch({
  results(dds_filtered, name = batch_name, independentFiltering = FALSE)
}, error = function(e) {
  warning('Failed to fetch batch main effect by name: ', conditionMessage(e))
  NULL
})

res_batch_XY <- NULL
if (!is.na(interaction_name)) {
  res_batch_XY <- tryCatch({
  ## Use the exact names discovered in resultsNames() to form the contrast list
  results(dds_filtered, contrast = list(c(batch_name, interaction_name)), independentFiltering = FALSE)
  }, error = function(e) {
    warning('Failed to fetch batch+interaction contrast: ', conditionMessage(e))
    NULL
  })
}

if (is.null(res_batch_XX) && is.null(res_batch_XY)) {
  stop('Could not extract batch effect results; please inspect resultsNames(dds_filtered) and available contrasts')
}

# Change to dataframe
df_batch_XX <- as.data.frame(res_batch_XX)
df_batch_XY <- as.data.frame(res_batch_XY)

# Add columns for scaling factors
df_batch_XX$Scaling_Factor_Batch_A <- sqrt(2^(df_batch_XX$log2FoldChange))
df_batch_XX$Scaling_Factor_Batch_B <- 1 / df_batch_XX$Scaling_Factor_Batch_A

df_batch_XY$Scaling_Factor_Batch_A <- sqrt(2^(df_batch_XY$log2FoldChange))
df_batch_XY$Scaling_Factor_Batch_B <- 1 / df_batch_XY$Scaling_Factor_Batch_A

# Extract normalized counts of genes
normalized_counts <- counts(dds_filtered, normalized = TRUE)
normalized_counts_df <- as.data.frame(normalized_counts)

# Create a copy of normalized counts to apply scaling
scaled_counts <- normalized_counts_df

# Iterate through each sample in normalized_counts_df
for (Sample_ID in colnames(normalized_counts_df)) {
  # Reset per-sample variable
  scaling_factors <- NULL
  
  # Get the 'sex' and 'batch' information for the current sample
  Sex <- as.character(colData[Sample_ID, "Sex"])
  Batch <- as.character(colData[Sample_ID, "Batch"])
  
  # Normalize casing to match the factor conventions used earlier
  Sex <- toupper(Sex)
  Batch <- tolower(Batch)
  
  # Determine the scaling factor column based on 'sex' and 'batch'
  if (Sex == "XX" && Batch == "a") {
    scaling_factors <- df_batch_XX$Scaling_Factor_Batch_A
  } else if (Sex == "XX" && Batch == "b") {
    scaling_factors <- df_batch_XX$Scaling_Factor_Batch_B
  } else if (Sex == "XY" && Batch == "a") {
    scaling_factors <- df_batch_XY$Scaling_Factor_Batch_A
  } else if (Sex == "XY" && Batch == "b") {
    scaling_factors <- df_batch_XY$Scaling_Factor_Batch_B
  }
  
  # If scaling_factors wasn't set due to unexpected sex/batch values, stop with a helpful message
  if (!exists("scaling_factors") || is.null(scaling_factors)) {
    stop(sprintf("Could not determine scaling factors for sample '%s' (Sex='%s', Batch='%s'). Check metadata and factor levels.", Sample_ID, Sex, Batch))
  }
  
  # Apply the scaling factor to each gene for the current sample
  # Match genes by their rownames in normalized_counts_df and the scaling factors
  scaled_counts[, Sample_ID] <- normalized_counts_df[, Sample_ID] * scaling_factors
}

# QA SECTION
# Two steps :
# 1 : re-run DESeq2 with the new scaled counts, we expect LFC for Batch and Batch interaction to be close to 0 for all genes.
# 2 : Plot PCA of day 0 data.

# Re-run DESeq2
# Since it expects integer, convert scaled counts to integers
scaled_counts_rounded <- round(scaled_counts)

# Ensure scaled counts are a numeric matrix and columns match colData rownames
scaled_counts_rounded <- as.matrix(scaled_counts_rounded)
mode(scaled_counts_rounded) <- "numeric"

# Check column/sample names
sample_cols <- colnames(scaled_counts_rounded)
meta_samples <- rownames(colData)
if (!all(sample_cols %in% meta_samples)) {
  missing_in_meta <- setdiff(sample_cols, meta_samples)
  missing_in_counts <- setdiff(meta_samples, sample_cols)
  if (length(missing_in_meta) > 0) {
    stop(sprintf("The following samples are in the counts but missing from metadata: %s", paste(missing_in_meta, collapse=", ")))
  }
  # If metadata has extra samples, reorder and subset counts to match metadata order
  if (length(missing_in_counts) > 0) {
    warning(sprintf("Metadata contains samples not present in counts; subsetting metadata to counts columns. Missing in counts: %s", paste(missing_in_counts, collapse=", ")))
    colData <- colData[sample_cols, , drop=FALSE]
  }
  # Recompute meta_samples after possible subsetting
  meta_samples <- rownames(colData)
}

# Reorder count columns to metadata order
scaled_counts_rounded <- scaled_counts_rounded[, meta_samples, drop=FALSE]

## Diagnostics: print dimensions and name checks to help debug SummarizedExperiment construction
cat('scaled_counts_rounded dims:', paste(dim(scaled_counts_rounded), collapse=' x '), '\n')
cat('first columns of counts:', paste(head(colnames(scaled_counts_rounded)), collapse=', '), '\n')
cat('metadata samples (rownames):', paste(head(rownames(colData)), collapse=', '), '\n')
cat('All count columns in metadata?', all(colnames(scaled_counts_rounded) %in% rownames(colData)), '\n')
cat('All metadata samples in count columns?', all(rownames(colData) %in% colnames(scaled_counts_rounded)), '\n')
if (!identical(colnames(scaled_counts_rounded), rownames(colData))) {
  mism1 <- setdiff(colnames(scaled_counts_rounded), rownames(colData))
  mism2 <- setdiff(rownames(colData), colnames(scaled_counts_rounded))
  cat('mismatch counts->meta:', paste(mism1, collapse=','), '\n')
  cat('mismatch meta->counts:', paste(mism2, collapse=','), '\n')
}

## If count columns look like simple numeric indices (e.g. "1","2",...), map them by position to metadata Sample_IDs
if (all(grepl('^\\d+$', colnames(scaled_counts_rounded)))) {
  if (ncol(scaled_counts_rounded) == nrow(colData)) {
    cat('Detected numeric column names in counts; mapping by position to metadata Sample_IDs.\n')
    colnames(scaled_counts_rounded) <- rownames(colData)
  } else {
    stop('Counts have numeric column names and number of columns does not match metadata rows; cannot map by position.')
  }
}

dds_scaled <- DESeqDataSetFromMatrix(countData = scaled_counts_rounded,
                              colData = colData,
                              design = design)

# Run DESeq2
dds_scaled <-DESeq(dds_scaled) # Remember : the result does not need to be filtered with HTSFilter, since the dataframe comes from previously filtered DESeq object

# Check LFC with corrected data
res_batch_XX_corrected <- results(dds_scaled, name = "Batch_b_vs_a", independentFiltering = FALSE) # batch effect on XX
res_batch_XY_corrected <- results(dds_scaled, contrast=list(c("Batch_b_vs_a","SexXY.Batchb")), independentFiltering = FALSE) # batch effect on XY (batch + interaction term)

# Change to dataframe
df_batch_XX_corrected <- as.data.frame(res_batch_XX_corrected)
df_batch_XY_corrected <- as.data.frame(res_batch_XY_corrected)

# All genes are 0.999 adjusted p value. Some are p-value significant but with very low counts

# PCA plot
# Identify samples to keep (PCA done on day 0 only)
samples_to_keep_day0 <- colData(dds_scaled)$Day == "0"

# Subset the dds_filtered object to include only the samples needed
dds_scaled_0 <- dds_scaled[, samples_to_keep_day0]
dds_filtered_0 <- dds_filtered[, samples_to_keep_day0]

# Create function for PCA
create_pca <- function(dds, plot_title) {
  
  # Perform VST transformation
  vsd <- vst(dds, blind = TRUE)
  # Blind has to be true since vst() requires a model matrix to be fitted, which can't be done with missing day 14 data. PCA do not change much either way
  
  # Create a grouping column in the colData for sample visualization
  colData_vsd <- as.data.frame(colData(vsd))
  colData_vsd$Group <- paste(colData_vsd$Treatment_Name, colData_vsd$Sex, sep = "_")
  colData(vsd) <- DataFrame(colData_vsd)
  
  # Generate PCA data
  pcaData <- plotPCA(vsd, intgroup = c("Group", "Cultivar"), returnData = TRUE)
  
  # Extract PCA coordinates
  pca_coords <- as.data.frame(pcaData[, 1:2])
  
  # Add Cultivar as a grouping variable
  pca_coords$Cultivar <-colData(vsd)$Cultivar
  pca_coords$Group <-colData(vsd)$Group
  
  # Calculate the percentage of variance explained for PC1 and PC2
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  # Define custom colors for groups
  custom_colors <- c("FF_XX" = "#33cc33",
                     "IFF_XY" = "#ff6600", 
                     "MF_XY" = "#993399", 
                     "IMF_XX" = "#0066ff")
  
  # Create PCA plot using ggplot2
  pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Group, group = Group)) +
    geom_point(aes(shape = Cultivar), size = 2) +
    stat_ellipse(type = "t", level = 0.95, linetype = 2) +
    scale_color_manual(values = custom_colors) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.3))
  
  return(pca_plot)
}

# Use the function for dds_filtered and dds_filtered_veg
pca_plot_scaled <- create_pca(dds_scaled_0, "PCA of VSD for all genes, day 0, batch corrected")
pca_plot_filtered <- create_pca(dds_filtered_0, "PCA of VSD for all genes, day 0, before correction")

# Print PCAs
print(pca_plot_scaled)
print(pca_plot_filtered)

# Final part to get VST Data for WGCNA
# Create VST data for all samples
vsd <- vst(dds_scaled, blind = FALSE) # By setting blind to FALSE, the dispersions already estimated will be used to perform transformations
# This seems to be the better method
# Setting blind = TRUE yields slightly different values
# Please refer to this section of the bioconductor package :
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#blind-dispersion-estimation

# The assay function is used to extract the matrix of normalized values
vsd_assay <- assay(vsd)

# DISCLAIMER BEFORE SAVING DATA
# The batch correction comes from effect at day 0 (only time point with controls for both batches for XX and XY)
# Thus, it assumes the batch effect is the same for day 1 and day 14
# There are no other ways to extract batch from these days, as any difference between batches can be attributed to treatments or photoperiod/maturity effect
# While this is reasonable for data at day 1, which is the same tissue type and was done only 24h later, it might be a bold assumption for day 14 data
# There are ways to check for this assumption, which I have looked into a bit, with unclear results
# There could be an argument for keeping counts uncorrected at day 14, and correcting at day 1
# This is simply a reminder about the structure of the experimental design
# DISCLAIMER OVER thanks for using this script

# Saving data
save(vsd_assay, file = "VSD_data_all_days_batch_corrected_26_11_V1.RData")
