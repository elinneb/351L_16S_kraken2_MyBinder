FROM rocker/binder:4.4.2

USER root

SHELL ["/bin/bash", "-c"]

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
    jq \
    bash \
    bash-completion \
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

RUN cat > /usr/local/bin/fetch_silva_db.sh <<'EOF'
#!/bin/bash
set -u

mkdir -p /home/rstudio/kraken_outputs
mkdir -p /home/rstudio/SILVA_DB
cd /home/rstudio

API_URL="https://api.figshare.com/v2/articles/22572970"
TMP_JSON="/tmp/figshare_22572970.json"

echo "Fetching Figshare article metadata..."
if ! curl -fL --retry 5 --retry-delay 2 -H "Accept: application/json" -o "$TMP_JSON" "$API_URL"; then
  echo "Could not fetch Figshare metadata."
  exit 0
fi

echo "Downloading files listed in article 22572970 ..."
jq -r '.files[] | [.name, .download_url] | @tsv' "$TMP_JSON" | while IFS=$'\t' read -r name url; do
  out="/home/rstudio/SILVA_DB/$name"
  echo "Downloading $name ..."
  curl -fL --retry 5 --retry-delay 2 -o "$out" "$url" || true
done

if [ -s /home/rstudio/SILVA_DB/hash.k2d ] && \
   [ -s /home/rstudio/SILVA_DB/opts.k2d ] && \
   [ -s /home/rstudio/SILVA_DB/taxo.k2d ]; then
  kraken2-inspect --db /home/rstudio/SILVA_DB | head -n 3 >/dev/null || true
fi

echo "SILVA_DB contents:"
ls -lh /home/rstudio/SILVA_DB || true
EOF

RUN chmod +x /usr/local/bin/fetch_silva_db.sh

RUN cat > /usr/local/bin/binder-entrypoint.sh <<'EOF'
#!/bin/bash
set -e

export SHELL=/bin/bash
export TERM=xterm-256color
export HOME=/home/rstudio

if [ ! -f "$HOME/.bashrc" ]; then
  cat > "$HOME/.bashrc" <<'EOBRC'
export SHELL=/bin/bash
export TERM=xterm-256color
export HISTFILE=~/.bash_history
export HISTSIZE=10000
export HISTFILESIZE=20000
export PROMPT_COMMAND='history -a'
[ -f /etc/bash_completion ] && . /etc/bash_completion
EOBRC
fi

touch "$HOME/.bash_history"
chmod 600 "$HOME/.bash_history" || true

/usr/local/bin/fetch_silva_db.sh || true

exec "$@"
EOF

RUN chmod +x /usr/local/bin/binder-entrypoint.sh
RUN chown -R rstudio:rstudio /home/rstudio
RUN chsh -s /bin/bash rstudio

USER rstudio
WORKDIR /home/rstudio

ENTRYPOINT ["/usr/local/bin/binder-entrypoint.sh"]
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser"]
