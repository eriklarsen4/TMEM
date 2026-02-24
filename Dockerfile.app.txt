# App Dockerfile
# Base image: TMEM (contains TMEM + Bioconductor deps)
FROM swingdoc/tmem:latest

# Copy Endo and Hippo source archives
COPY Endo_*.zip /tmp/
# COPY Hippo_*.zip /tmp/

# Install Endo and Hippo
# RUN R -e "install.packages('/tmp/Endo_1.0.0.zip', repos = NULL, type = 'source')"
RUN R -e "install.packages('/tmp/Hippo_1.0.0.zip', repos = NULL, type = 'source')"

# Copy the Shiny app into the image
# COPY TMEMapp/ /srv/shiny-server/

# Default command: run the Shiny app
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server')"]