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

RUN cat > /home/rstudio/fetch_silva_db.sh <<'EOF'
#!/bin/bash
set -u

mkdir -p /home/rstudio/kraken_outputs
mkdir -p /home/rstudio/SILVA_DB
cd /home/rstudio

UA="Mozilla/5.0"
FAILED=0

download() {
  local out="$1"
  local url="$2"

  echo "Downloading $out ..."
  if curl -A "$UA" -fL --retry 5 --retry-delay 2 -o "$out" "$url"; then
    if [ ! -s "$out" ]; then
      echo "Empty file: $out"
      FAILED=1
    elif file "$out" | grep -qi html; then
      echo "Bad download (HTML): $out"
      FAILED=1
    else
      echo "OK: $out"
    fi
  else
    echo "Failed: $out"
    FAILED=1
  fi
}

download SILVA_DB/database50mers.kmer_distrib  "https://figshare.com/ndownloader/files/40054801?download=1"
download SILVA_DB/database75mers.kmer_distrib  "https://figshare.com/ndownloader/files/40054804?download=1"
download SILVA_DB/database100mers.kmer_distrib "https://figshare.com/ndownloader/files/40054807?download=1"
download SILVA_DB/database150mers.kmer_distrib "https://figshare.com/ndownloader/files/40054810?download=1"
download SILVA_DB/database200mers.kmer_distrib "https://figshare.com/ndownloader/files/40054813?download=1"
download SILVA_DB/database250mers.kmer_distrib "https://figshare.com/ndownloader/files/40054816?download=1"

download SILVA_DB/hash.k2d   "https://figshare.com/ndownloader/files/40054819?download=1"
download SILVA_DB/opts.k2d   "https://figshare.com/ndownloader/files/40054822?download=1"
download SILVA_DB/README.md  "https://figshare.com/ndownloader/files/40054825?download=1"
download SILVA_DB/taxo.k2d   "https://figshare.com/ndownloader/files/40054828?download=1"

if [ -s SILVA_DB/hash.k2d ] && [ -s SILVA_DB/opts.k2d ] && [ -s SILVA_DB/taxo.k2d ]; then
  kraken2-inspect --db SILVA_DB | head -n 3 >/dev/null || true
fi

if [ "$FAILED" -ne 0 ]; then
  echo
  echo "Some SILVA_DB files failed to download."
  echo "Re-run this after Binder starts:"
  echo "  bash /home/rstudio/fetch_silva_db.sh"
  exit 1
fi

echo "All SILVA_DB files downloaded successfully."
EOF

RUN chmod +x /home/rstudio/fetch_silva_db.sh

RUN /home/rstudio/fetch_silva_db.sh || echo "Partial/failed DB download during build; run bash /home/rstudio/fetch_silva_db.sh after Binder starts."
