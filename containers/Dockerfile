FROM rocker/r-ver:4.3

WORKDIR /app

COPY . /app

RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); pkgs <- trimws(readLines('r-packages.txt')); pkgs <- pkgs[pkgs != '' & !startsWith(pkgs, '#')]; BiocManager::install(pkgs)"

LABEL tool="phosr_maxquant_tmt_modular.R"
