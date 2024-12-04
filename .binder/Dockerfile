FROM rocker/binder:4.4.2
LABEL maintainer='Abu_Saadat'
COPY --chown=${NB_USER} . ${HOME}
USER ${NB_USER}



RUN wget https://github.com/SAADAT-Abu/RNASeq_Lecture/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2024-12-04/')); devtools::install_deps()"

RUN rm DESCRIPTION.1; exit 0
