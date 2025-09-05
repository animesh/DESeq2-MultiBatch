#!/usr/bin/env Rscript
# Full end-to-end integration test: generate balanced dataset, run pipeline, validate outputs

message('Starting full integration test (this runs DESeq2 and may take some time)')

# generate data
source('scripts/generate_full_integration_data.R')
data_dir <- file.path('tests', 'fixtures', 'full_integration', 'data')
out_dir <- file.path(tempdir(), 'deseq2_full_integration_outputs')
if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)
dir.create(out_dir, recursive = TRUE)

# run pipeline (full run)
res <- system2('Rscript', args = c('scripts/run_deseq2_multibatch.R', data_dir, out_dir), stdout = TRUE, stderr = TRUE)
cat(paste(res, collapse = '\n'), '\n')

# check expected outputs
expected <- c('dds_single.rds', 'dds_double.rds', 'batch_single_results.csv', 'normalized_counts_single.csv')
for (f in expected) {
  fp <- file.path(out_dir, f)
  if (!file.exists(fp)) stop('Missing expected output: ', fp)
  if (file.info(fp)$size == 0) stop('Empty output file: ', fp)
}

message('Full integration test passed. Outputs at: ', normalizePath(out_dir))
