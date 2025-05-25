FROM rocker/rstudio:4.3.2

ENV USER=rstudio
ENV PASSWORD=tmem

COPY . /C:/Users/Erik/Desktop/Programming/R/Bio/TMEM
WORKDIR /C:/Users/Erik/Desktop/Programming/R/Bio/TMEM

EXPOSE 8787