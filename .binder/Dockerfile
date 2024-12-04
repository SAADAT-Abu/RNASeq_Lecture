FROM rocker/binder:4.4.1
LABEL maintainer='Abu_Saadat'
COPY --chown=${NB_USER} . ${HOME}
USER ${NB_USER}



RUN R -e "options(repos = list(CRAN = 'https://cloud.r-project.org')); \
          install.packages(c('matrixStats', 'tidyverse', 'ggplot2', 'plotly', 'ellipse', 'corpcor', 'reshape2', 'matrixStats', 'rARPACK', 'gridExtra', 'ggrepel')); \
          if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
          BiocManager::install(c('MatrixGenerics','SparseArray','BiocParallel', 'DESeq2', 'GenomicRanges', 'mixOmics'))"


