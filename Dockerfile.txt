FROM rocker/binder:4.4.2

USER root

RUN install2.r --error \
    BiocManager \
    remotes \
    tidyverse \
    vegan \
    zCompositions \
    ggforce \
    cowplot \
    RColorBrewer \
    broom \
    pheatmap \
    rmdformats

RUN R -e "BiocManager::install(c('phyloseq','rhdf5'), ask = FALSE, update = FALSE)"

RUN R -e "remotes::install_version('microeco', version = '2.0.0', repos = 'https://cloud.r-project.org')"

RUN R -e "remotes::install_github('ChiLiubio/file2meco')"

RUN R -q -e "stopifnot( \
  requireNamespace('phyloseq', quietly = TRUE), \
  requireNamespace('rhdf5', quietly = TRUE), \
  requireNamespace('microeco', quietly = TRUE), \
  requireNamespace('file2meco', quietly = TRUE), \
  requireNamespace('rmdformats', quietly = TRUE) \
)"

COPY . /home/rstudio
RUN chown -R rstudio:rstudio /home/rstudio

USER rstudio
WORKDIR /home/rstudio
