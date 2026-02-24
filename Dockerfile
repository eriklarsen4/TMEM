# TMEM Dockerfile
# Base image: Bioconductor + R (stable, multi-arch compatible)
# Provides BiocManager, annotation repos, and system libs needed for TMEM.
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install system dependencies required by TMEM and its Bioconductor packages
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy the TMEM source tree into the container
# (Binary artifacts are excluded via .dockerignore)
COPY TMEM /TMEM
WORKDIR /TMEM

# Install Bioconductor annotation packages required by TMEM
RUN R -e "install.packages('BiocManager'); \
          BiocManager::install(version = '3.18', ask = FALSE); \
          BiocManager::install(c('GO.db', 'AnnotationDbi', \
                                 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Dm.eg.db'));"

# Install TMEM's CRAN dependencies
RUN R -e "install.packages(c('dplyr','tidyr','purrr','stringr','rlang','assertthat'), \
                           repos='https://cloud.r-project.org')"

# Install orthogene from CRAN (grr dependency removed upstream)
RUN R -e "install.packages('orthogene', repos='https://cloud.r-project.org')"

# Install TMEM from the local source directory
# (This replaces the old .tar.gz install step)
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org'); \
          remotes::install_local('.', upgrade = 'never')"

# Default command (optional)
CMD ["bash"]