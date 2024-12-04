FROM rocker/binder:4.4.1
LABEL maintainer='Abu_Saadat'
COPY --chown=${NB_USER} . ${HOME}
USER ${NB_USER}



RUN R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2024-12-04/')); \
          install.packages(c('tidyverse', 'ggplot2', 'plotly')); \
          if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
          BiocManager::install(c('DESeq2', 'GenomicRanges', 'mixOmics'))"

RUN rm DESCRIPTION.1; exit 0
