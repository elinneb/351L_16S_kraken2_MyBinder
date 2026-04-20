FROM rocker/binder:4.4.2

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    nano \
    tzdata \
    libncurses5-dev \
    libncursesw5-dev \
    parallel \
    fastqc \
    trimmomatic \
    multiqc \
    kraken2 \
    python3 \
    python3-pip \
    curl \
    file \
    wget \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir \
    numpy \
    pandas \
    biom-format \
    kraken-biom \
    notebook \
    jupyterlab

RUN install2.r --error --skipinstalled \
    cowplot \
    RColorBrewer \
    ggplot2 \
    remotes

RUN R -q -e "install.packages('microeco', repos='https://cloud.r-project.org')"
RUN R -q -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"
RUN R -q -e "BiocManager::install(c('phyloseq','biomformat','rhdf5'), ask=FALSE, update=FALSE)"
RUN R -q -e "remotes::install_github('ChiLiubio/file2meco')"

COPY . /home/rstudio
RUN chown -R rstudio:rstudio /home/rstudio

USER rstudio
WORKDIR /home/rstudio

RUN mkdir -p kraken_outputs SILVA_DB

RUN cat > /home/rstudio/fetch_silva_db.sh <<'EOS'
#!/bin/bash
set -euo pipefail

cd /home/rstudio
mkdir -p kraken_outputs SILVA_DB

UA="Mozilla/5.0"
fetch() {
  local out="$1"
  local url="$2"

  rm -f "$out"

  curl -A "$UA" -fL --retry 5 --retry-all-errors --connect-timeout 30 -o "$out" "$url" || \
  wget --user-agent="$UA" --tries=5 --timeout=30 -O "$out" "$url"

  file "$out" | grep -qi 'html' && { echo "Bad download (HTML): $out"; exit 1; }
  [ -s "$out" ] || { echo "Empty file: $out"; exit 1; }
}

fetch SILVA_DB/database50mers.kmer_distrib  "https://figshare.com/ndownloader/files/40054801?download=1"
fetch SILVA_DB/database75mers.kmer_distrib  "https://figshare.com/ndownloader/files/40054804?download=1"
fetch SILVA_DB/database100mers.kmer_distrib "https://figshare.com/ndownloader/files/40054807?download=1"
fetch SILVA_DB/database150mers.kmer_distrib "https://figshare.com/ndownloader/files/40054810?download=1"
fetch SILVA_DB/database200mers.kmer_distrib "https://figshare.com/ndownloader/files/40054813?download=1"
fetch SILVA_DB/database250mers.kmer_distrib "https://figshare.com/ndownloader/files/40054816?download=1"
fetch SILVA_DB/hash.k2d "https://figshare.com/ndownloader/files/40054819?download=1"
fetch SILVA_DB/opts.k2d "https://figshare.com/ndownloader/files/40054822?download=1"
fetch SILVA_DB/README.md "https://figshare.com/ndownloader/files/40054825?download=1"
fetch SILVA_DB/taxo.k2d "https://figshare.com/ndownloader/files/40054828?download=1"

kraken2-inspect --db SILVA_DB | head -n 3 >/dev/null
EOS

RUN chmod +x /home/rstudio/fetch_silva_db.sh

# Try to prefetch during image build, but do not fail the whole image if Figshare is flaky.
RUN /home/rstudio/fetch_silva_db.sh || echo "Figshare download skipped during build; run fetch_silva_db.sh inside the Binder session if needed."
