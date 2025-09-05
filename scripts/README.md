# Scripts

This folder contains a runnable R script implementing the workflow described in the repository `README.md`.

Script: `run_deseq2_multibatch.R`

Usage:

1. Place `coldata.csv` and `counts.csv` files in the `data/` directory (this repo already includes those files).
2. From the repository root run the script with Rscript:

```bash
Rscript scripts/run_deseq2_multibatch.R data outputs
```

The script will:

- Install (via BiocManager) and load `DESeq2` if it is not available.
- Run two DESeq2 analyses (single and double batch designs).
- Compute scaling factors from batch contrasts and apply per-gene scaling to normalized counts.
- Write outputs into the `outputs/` folder: normalized counts, scaled counts, saved DDS objects and CSVs with batch results.

Notes:

- The scaled counts are for inspection and plotting only; do not reuse them for DE analysis without careful consideration (see README).
- If `resultsNames(dds_double)` uses different coefficient names than assumed, the script will fallback to using base batch effect and will print a warning; inspect `resultsNames(dds_double)` for custom names.
