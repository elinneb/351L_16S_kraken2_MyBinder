FROM rocker/binder:4.4.2

USER root

# system packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    nano \
    tzdata \
    libncurses5-dev \
    libncursesw5-dev \
    parallel \
    fastqc \
    trimmomatic \
    multiqc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# conda packages already available in rocker/binder via mamba
RUN mamba install -y -c conda-forge -c bioconda \
    kraken2 \
    kraken-biom \
    numpy \
    pandas \
    biom-format \
    bioconductor-phyloseq \
    bioconductor-biomformat \
    bioconductor-rhdf5 \
    r-cowplot \
    r-rcolorbrewer \
    r-ggplot2 \
    && mamba clean -afy

# R packages not conveniently available above
RUN install2.r --error \
    remotes

RUN R -e "options(repos='https://cloud.r-project.org'); install.packages(c('microeco'))"
RUN R -e "remotes::install_github('ChiLiubio/file2meco')"

# copy repo contents into RStudio user's home
COPY . /home/rstudio
RUN chown -R rstudio:rstudio /home/rstudio

USER rstudio
WORKDIR /home/rstudio
