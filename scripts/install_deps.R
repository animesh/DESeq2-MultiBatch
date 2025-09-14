#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://cloud.r-project.org"))
cat("Starting dependency installation...\n")

tryCatch({
  # CRAN packages used in scripts/orig.r
  cran_pkgs <- c("dplyr", "ggplot2")
  to_install_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
  if (length(to_install_cran) > 0) {
    cat("Installing CRAN packages:", paste(to_install_cran, collapse = ", "), "\n")
    install.packages(to_install_cran)
  } else {
    cat("All CRAN packages already installed.\n")
  }

  # Ensure BiocManager is available
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager")
  }

  # Bioconductor packages needed by the analysis
  bioc_pkgs <- c("DESeq2", "HTSFilter")
  cat("Installing Bioconductor packages:", paste(bioc_pkgs, collapse = ", "), "\n")
  BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

  cat("Dependency installation completed.\n")
}, error = function(e) {
  cat("ERROR during installation:\n")
  cat(conditionMessage(e), "\n")
  quit(status = 1)
})

invisible(NULL)
options(repos = c(CRAN="https://cloud.r-project.org"))
pkgs <- c("dplyr","ggplot2")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# Install Bioconductor packages non-interactively
BiocManager::install(c("DESeq2","HTSFilter"), ask=FALSE, update=FALSE)
# Print session info to help debugging
sessionInfo()
