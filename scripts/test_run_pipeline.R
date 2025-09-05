#!/usr/bin/env Rscript
# Integration test: create temp CSVs, run run_deseq2_multibatch.R, validate outputs

message('Starting integration test for run_deseq2_multibatch.R')

tmp_root <- tempfile('deseq2_integration_')
dir.create(tmp_root)
data_dir <- file.path(tmp_root, 'data')
out_dir <- file.path(tmp_root, 'outputs')
dir.create(data_dir)
dir.create(out_dir)

# create tiny counts matrix (genes x samples)
genes <- paste0('Gtest', 1:6)
# create 4 samples, columns named 1..4 to test numeric->Sample mapping
counts_mat <- matrix(c(
  10,20,30,40,
  5,10,5,10,
  50,60,55,65,
  0,1,0,2,
  100,110,95,120,
  3,3,2,1
), nrow = length(genes), byrow = TRUE)
rownames(counts_mat) <- genes
colnames(counts_mat) <- as.character(1:4)
counts_df <- as.data.frame(counts_mat)
write.csv(counts_df, file = file.path(data_dir, 'counts.csv'), row.names = TRUE)

# create coldata with Sample1..Sample4 as rownames and design columns
coldata <- data.frame(
  Genotype = factor(c('g1','g1','g2','g2')),
  Sex = factor(c('xx','xy','xx','xy')),
  Day = factor(c('0','1','0','1')),
  Batch = factor(c('a','b','a','b')),
  Treatment = factor(c('ctrl','t1','t1','ctrl'))
)
rownames(coldata) <- paste0('Sample', 1:4)
write.csv(coldata, file = file.path(data_dir, 'coldata.csv'), row.names = TRUE)

# run the pipeline script in dry-run mode to validate alignment only
cmd <- c('scripts/run_deseq2_multibatch.R', data_dir, out_dir, '--dry-run')
res <- system2('Rscript', args = cmd, stdout = TRUE, stderr = TRUE)
cat(paste(res, collapse = '\n'), '\n')

# validate that alignment marker exists
marker <- file.path(out_dir, 'alignment_ok.txt')
if (!file.exists(marker)) stop('Expected alignment marker missing: ', marker)
if (file.info(marker)$size <= 0) stop('Alignment marker is empty: ', marker)

message('Integration dry-run test passed: alignment marker present at ', marker)

invisible(NULL)
