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

# Copy the TMEM package source archive into the container
# (Assumes TMEM_1.0.0.zip is built before docker build)
COPY TMEM_1.0.0.zip /tmp/TMEM_1.0.0.zip

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

# Install TMEM from the copied source archive
RUN R -e \"install.packages('/tmp/TMEM_1.0.0.zip', repos = NULL, type = 'source')\"

# Remove temporary files
RUN rm -f /tmp/TMEM_1.0.0.zip

# Default command (optional)
CMD ["bash"]