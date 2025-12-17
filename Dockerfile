
# copy the base image
FROM rocker/rstudio:4.3.2

# copy the local tmem directory for the container
COPY TMEM_1.0.0.zip /C:/Users/Erik/Desktop/Programming/R/Bio/TMEM/TMEM_1.0.0.zip

# install the TMEM package
RUN R -e "install.packages('/C:/Users/Erik/Desktop/Programming/R/Bio/TMEM/TMEM_1.0.0.zip', repos = NULL, type = 'source')"
