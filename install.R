# Install CRAN packages
install.packages(c("tidyverse", "ggplot2", "plotly"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("DESeq2", "GenomicRanges", "mixOmics"))
