#!/usr/bin/env Rscript
# Generate a small balanced, non-collinear synthetic dataset for full integration tests

set.seed(42)

out_root <- file.path('tests', 'fixtures', 'full_integration')
data_dir <- file.path(out_root, 'data')
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# design: full factorial of Genotype x Sex x Day x Batch with 2 replicates (2x2x2x2x2 = 32 samples)
Genotype <- c('g1', 'g2')
Sex <- c('xx', 'xy')
Day <- c('0', '1')
Batch <- c('a', 'b')
Replicate <- 1:2

design <- expand.grid(Genotype = Genotype, Sex = Sex, Day = Day, Batch = Batch, Replicate = Replicate, stringsAsFactors = FALSE)
# assign Treatment randomly but balanced: half ctrl, half treat
n <- nrow(design)
tr <- rep(c('ctrl','treat'), length.out = n)
tr <- sample(tr, n) # shuffle to avoid deterministic alignment
design$Treatment <- tr
rownames(design) <- paste0('Sample', seq_len(nrow(design)))

coldata <- design[, c('Genotype','Sex','Day','Batch','Treatment')]
write.csv(cbind(SampleID = rownames(coldata), coldata), file = file.path(data_dir, 'coldata.csv'), row.names = FALSE)

# create counts: 50 genes, deterministic but includes batch/sex effects to be realistic
n_genes <- 50
genes <- paste0('G', seq_len(n_genes))
nsamples <- nrow(coldata)
counts <- matrix(0, nrow = n_genes, ncol = nsamples)
colnames(counts) <- rownames(coldata)
rownames(counts) <- genes

for (g in seq_len(n_genes)) {
  base <- 100 + g * 3
  for (s in seq_len(nsamples)) {
    samp <- rownames(coldata)[s]
    sex <- coldata[s, 'Sex']
    batch <- coldata[s, 'Batch']
    # add simple effects: odd genes affected by batch, genes divisible by 3 affected by sex
    val <- base
    if ((g %% 2) == 1 && batch == 'b') val <- val + 20
    if ((g %% 3) == 0 && sex == 'xy') val <- val + 15
    # small gene-specific random noise (Poisson around val)
    counts[g, s] <- rpois(1, lambda = max(1, val))
  }
}

write.csv(as.data.frame(counts), file = file.path(data_dir, 'counts.csv'), row.names = TRUE)

cat('Wrote synthetic full-integration data to:', normalizePath(data_dir), '\n')