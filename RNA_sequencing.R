onLoad <- function(libname, pkgname) {
  # List of Bioconductor packages package depends on
  bioc_packages <- c("GenomicRanges","mixOmics","DESeq2") 
  # Check if BiocManager is available; install it if not
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Install missing Bioconductor dependencies
  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE)
    }
  }
  
}

## ----message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------

library(tidyverse)      # for data manipulation
library(DESeq2)         # for Differential gene expression analysis
library(ggplot2)        # for plotting
library(plotly)         # for interactive plots
library(GenomicRanges)  # for creating and manipulating genomic ranges
library(mixOmics)       # for PCA


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

assay<- read.csv("assay.csv", row.names = 1)  # Importing gene expression assay data
colData <- read.csv("colData.csv")            # Importing clinical information
temp<- read.csv("rawRange.csv")               # Importing rawRanges as dataframe

rowRanges <- GRanges(                                             # creating rawRanges from dataframge
             seqnames = Rle(temp$chr),                            # Chromosome names
             ranges = IRanges(start = temp$start, end =temp$end), # Start and End positions
             strand = Rle(temp$strand),                           # Strand information (if available)
             data.frame(Gene = temp$Gene, Length = temp$length))  # Additional columns



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

print("The count data as Assay: ")
head(assay)

print("The clinical data data as colData: ")
head(colData)

print("The gene coordinates as rawRanges ")
head(rowRanges)



## ----message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------

nreads <- data.frame(reads = colSums(assay)) %>%                   # Calculating library size by summing read counts for each sample
  rownames_to_column("sample") %>%
  mutate(Sample_Group = case_when(grepl("HK", sample) ~ "HK",
                         grepl("WT", sample) ~ "WT")) 

lsize <- ggplot(nreads, aes(x=Sample_Group, y=reads, fill=Sample_Group,label=sample)) +
  geom_boxplot() +
  geom_jitter(width = 0.25) +
  scale_fill_manual(values = c("forest green", "steel blue")) +
  theme_minimal() + 
  ggtitle("Library size")


ggplotly(lsize)




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

assay.filtered <- assay[rowSums(assay > 0) >= 2 & rowSums(assay) >= 20,]

paste0("Number of genes before filtering: ", nrow(assay))
paste0("Number of genes after filtering: ", nrow(assay.filtered))



## ----message=FALSE , warning=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------

# Add pseudocount and calculate log2-transformed counts
log_counts <- log2(assay.filtered + 1)

# Calculate relative log expression
medians <- apply(log_counts, 1, median)         # Median for each gene across samples
rle <- sweep(log_counts, 1, medians, FUN = "-") # Subtract medians (log scale)

# Reshape data for ggplot
rle_long <- reshape2::melt(rle)
colnames(rle_long) <- c("Sample", "Deviation")
rle_long$Sample_group <-  rep(c("HK", "WT"), c(1005,1005))

# Plot RLE as a boxplot colored by Sample_Group
ggplot(rle_long, aes(x = Sample, y = Deviation, fill = Sample_group)) +
  geom_boxplot(outlier.size = 0.5, outlier.colour = "red") +
  scale_fill_manual(values = c("forest green", "steel blue")) +
  theme_minimal() +
  labs(
    title = "RLE Plot before Normalization",
    x = "Sample",
    y = "Relative Log Expression (Deviation)",
    fill = "Sample_group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for clarity



## ----message=FALSE , warning=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------

DE.obj <- DESeqDataSetFromMatrix(countData = assay.filtered,
                                 colData   = colData,
                                 design    = ~ Sample_Group)



## ----message=FALSE , warning=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------

DE.obj <- estimateSizeFactors(DE.obj)            # Estimating size factor
 
DE.normalized <- counts(DE.obj, normalized = T)  # Normalizing counts

# Add pseudocount and calculate log2-transformed counts
log_counts <- log2(DE.normalized + 1) %>% as.data.frame()

# Calculate relative log expression
medians <- apply(log_counts, 1, median)         # Median for each gene across samples
rle <- sweep(log_counts, 1, medians, FUN = "-") # Subtract medians (log scale)

# Reshape data for ggplot
rle_long <- reshape2::melt(rle)
colnames(rle_long) <- c("Sample", "Deviation")
rle_long$Sample_group <-  rep(c("HK", "WT"), c(1005,1005))

# Plot RLE as a boxplot colored by Sample_Group
ggplot(rle_long, aes(x = Sample, y = Deviation, fill = Sample_group)) +
  geom_boxplot(outlier.size = 0.5, outlier.colour = "red") +
  scale_fill_manual(values = c("forest green", "steel blue")) +
  theme_minimal() +
  labs(
    title = "RLE Plot after Normalization",
    x = "Sample",
    y = "Relative Log Expression (Deviation)",
    fill = "Sample_group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for clarity



## ----message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------

pca1<- pca(t(DE.normalized), center = T, scale = T)       # Performing PCA

plotIndiv(pca1, group = colData$Sample_Group, ind.names = F, 
          title = "Total expression profile",
          legend = T, legend.position = "bottom", col.per.group = c("forest green", "steel blue"))



## ----message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------

DE.obj$Sample_Group <- relevel(DE.obj$Sample_Group, ref = "HK")  # Setting reference group

DE.obj <- DESeq(DE.obj)                                          #Fitting the model

res.con <- results(DE.obj,                                       # Extracting result table
                   contrast=c("Sample_Group", "WT", "HK"),       # Setting Contrast 
                   pAdjustMethod = "BH",                         # Using 'BH' for p-value adjustment
                   alpha = 0.1)                                  # Adjusted p-value threshold
                   
res.con <- res.con[!is.na(res.con$padj),]                        # Removing genes with NA adj p-values                  




## ----message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------

volcano <- as.data.frame(res.con) %>%
  rownames_to_column("gene") %>%
  mutate(Logpadj=-log10(padj)) %>%
  mutate(sig=case_when(log2FoldChange >= 1 & padj < 0.05 ~ "Up",
                       log2FoldChange <= -1 & padj < 0.05 ~ "Down")) %>%
  ggplot(aes(x=log2FoldChange, y=Logpadj, color=sig, label=gene)) + 
  geom_point(size=1) +
  theme_bw() + 
  ylab("-Log10(FDR)") + 
  theme(legend.position = "none") +
  scale_color_manual(values = c("firebrick", "steelblue")) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') + 
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  ggtitle("Volcano DEGs")
  

ggplotly(volcano)



