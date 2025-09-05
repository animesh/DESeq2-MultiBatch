batch_xx$Scaling_Factor_Batch_A <- sqrt(2^(batch_xx$log2FoldChange))
batch_xy$Scaling_Factor_Batch_A <- sqrt(2^(batch_xy$log2FoldChange))
# DESeq2-MultiBatch

Batch Correction for Multi-Factorial RNA-seq Experiments

Overview
--------
This repository implements the DESeq2-MultiBatch pipeline described in the paper and provides a runnable R script plus tests that exercise the important data-alignment and analysis steps.

This README was updated to document additions made in the repository:
- A runnable pipeline: `scripts/run_deseq2_multibatch.R` (reads `data/coldata.csv` and `data/counts.csv`, runs DESeq2, computes batch contrasts and scaling factors, applies per-gene scaling, and writes outputs).
- Robust sample-name alignment logic: handles common mismatches (numeric column names, leading `X`, punctuation/case differences) and reports helpful errors.
- Dry-run support for quickly validating sample alignment without running DESeq2 (`--dry-run` flag).
- Unit and integration tests:
  - `scripts/test_align_names.R` — unit tests for alignment function.
  - `scripts/test_run_pipeline.R` — integration dry-run test that writes small temp CSVs, runs the pipeline in dry-run mode, and verifies an alignment marker is produced.

Quick links
-----------
- Main pipeline: `scripts/run_deseq2_multibatch.R`
- Unit tests: `scripts/test_align_names.R`
- Integration test: `scripts/test_run_pipeline.R`

Requirements
------------
- Linux/macOS/Windows with R (>= 4.0 recommended).
- Internet access to install Bioconductor packages (unless pre-installed).
- R packages: `DESeq2` and its Bioconductor dependencies. The pipeline will attempt to install missing packages via `BiocManager::install()` when run.

Recommended environment setup (fast path)
----------------------------------------
To avoid long package installation times during testing, pre-install the main dependencies interactively once:

1) Start R, then:

```r
install.packages('BiocManager')
BiocManager::install(c('DESeq2'))
```

Or use your institution/container image with DESeq2 already installed.

Reproducing the analysis (step-by-step)
--------------------------------------
Assume you have a working R environment and this repository checked out. The repository contains a sample `data/` with `coldata.csv` and `counts.csv`.

1) Dry-run (alignment check only — fast)

This verifies sample-name alignment and mapping logic without fitting DESeq2 (suitable for CI and quick checks):

```bash
Rscript scripts/run_deseq2_multibatch.R data outputs --dry-run
```

On success a small marker file will be written:

```
outputs/alignment_ok.txt
```

2) Run the full pipeline (may install packages and take longer)

```bash
Rscript scripts/run_deseq2_multibatch.R data outputs
```

This will:
- Read `data/coldata.csv` (expects sample IDs in the first column which become rownames; the script is tolerant to common naming mismatches).
- Read `data/counts.csv` (rows = genes, columns = samples).
- Run DESeq2 with two designs (single and double batch interaction), save `dds_single.rds` and `dds_double.rds` in `outputs/`.
- Extract batch contrasts and compute per-gene scaling factors.
- Produce normalized and scaled counts CSVs and result CSVs in `outputs/`.

Produced files (written to the `outputs/` directory):

- `dds_single.rds` — saved DESeqDataSet after DESeq() for the single design.
- `dds_double.rds` — saved DESeqDataSet after DESeq() for the double design.
- `batch_single_results.csv` — batch contrast results and scaling factors for the single design.
- `batch_double_xx_results.csv`, `batch_double_xy_results.csv` — sex-specific/double-design batch results and scaling factors.
- `normalized_counts_single.csv`, `normalized_counts_double.csv` — normalized counts from DESeq2 for each design.
- `scaled_counts_single.csv`, `scaled_counts_double.csv` — normalized counts multiplied by per-gene scaling factors.
- `alignment_ok.txt` — (only if run with `--dry-run`) contains a short alignment report.

Notes, caveats and implementation specifics
-----------------------------------------
- Sample-name alignment: the pipeline will attempt multiple mapping strategies when column names in `counts.csv` do not exactly match `coldata` rownames:
  1) Strip a leading `X`.
  2) If a column name is purely numeric (e.g., `1`), it will try `Sample1` mapping.
  3) Relaxed matching by lowercasing and stripping non-alphanumeric characters.
  If no mapping succeeds, the script stops with a helpful diagnostic message.

- Design variables: the README's original example used `Day` which may be numeric. DESeq2 will treat numeric covariates as continuous by default — if Day is categorical/timepoint, convert to factor in your `coldata` (the script converts character columns to factors but will not coerce integer columns automatically). If you pass a design that is collinear (e.g., small toy datasets), DESeq2 will error with "model matrix not full rank". Use the `--dry-run` mode for alignment checks and create a larger, non-collinear synthetic dataset for full pipeline tests.

- Scaling approach: per-gene scaling factors are computed from the batch log2 fold change (LFC):

  Scaling_Factor_Batch_A = sqrt(2^(log2FoldChange))
  Scaling_Factor_Batch_B = 1 / Scaling_Factor_Batch_A

  This follows the approach described in the original README/paper. Scaled counts are intended as a corrected table for downstream inspection or visualization, not for re-running DESeq2 without careful consideration — p-values will be invalid if you naively re-run DESeq2 on scaled counts.

Tests included
--------------
- Unit tests: `scripts/test_align_names.R` — run with:

```bash
Rscript scripts/test_align_names.R
```

- Integration dry-run: `scripts/test_run_pipeline.R` — runs a tiny synthetic dataset through the pipeline in `--dry-run` mode and checks for `alignment_ok.txt`:

```bash
Rscript scripts/test_run_pipeline.R
```

CI suggestions
--------------
- Add a small GitHub Actions workflow that runs the unit test and the integration dry-run test on push/PR. Avoid running the full DESeq2 pipeline in shared CI unless the runner has DESeq2 preinstalled or caching is configured.

How to contribute
-----------------
- If you want to improve the alignment heuristics, move the helper into `scripts/utils.R` and source it from both the pipeline and the tests.
- If you want true end-to-end integration tests, provide a small non-collinear synthetic dataset in `tests/fixtures/` and update the integration test to run the full pipeline (without `--dry-run`).

Contact / Citation
------------------
If used in research, please cite the original preprint:

Roy, J., Monthony, A. S., Torkamaneh, D. (2025). DESeq2-MultiBatch: Batch Correction for Multi-Factorial RNA-seq Experiments. https://doi.org/10.1101/2025.04.20.649392

For more detail and figures, see the HTML file in this repository (being updated).
