<!-- badges: start -->
  [![Launch Rstudio Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SAADAT-Abu/RNASeq_Lecture/master?urlpath=rstudio)
  <!-- badges: end -->
# RNA-Seq Analysis with DESeq2

This repository contains the script and data required for RNA-Seq differential expression analysis using the DESeq2 package in R. The script guides you through the steps of data preprocessing, normalization, and differential expression analysis.

---

## Contents

1. `RNA_sequencing.Rmd` - R Markdown script for RNA-Seq analysis.
2. `/data/` - Folder containing input data files:
   - `raw_counts.csv`: Raw count matrix (genes × samples).
   - `metadata.csv`: Metadata for samples, including experimental conditions.

---

## Requirements

### Software
- R (version ≥ 4.0)
- RStudio (optional, for easier execution)

### R Packages
The following R packages are required:

- tidyverse
- DESeq2
- ggplot2
- plotly
- GenomicRanges
- mixOmics

Install these packages using the following commands in R:
```R

install.packages(c("tidyverse", "ggplot2", "plotly"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("DESeq2", "GenomicRanges", "mixOmics"))

```

### Clone this Repository:

```bash

git clone https://github.com/your-username/RNA_Seq_Analysis.git
cd RNA_Seq_Analysis

```
