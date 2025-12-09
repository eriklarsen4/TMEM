  <!-- badges: start -->
  [![R-CMD-check](https://github.com/eriklarsen4/TMEM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/eriklarsen4/TMEM/actions/workflows/R-CMD-check.yaml)
  ![Static Badge](https://img.shields.io/badge/MBNeuroLab-darkblue)
  ![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
  <!-- badges: end -->

**TMEM**

Repository and package for bioinformatics- and molecular biology-related
functions, facilitating biological data exploration and easier plotting of
downstream NGS analyses (heatmaps, bargraphs, volcano plots, MA plots, scatter
plots)

Adapted from source code used to analyze data in **Dr. Martha Bhattacharya's lab**
at the **University of Arizona**

  + The lab studies transmembrane protein 184b (TMEM184B), a protein involved in
  axon degeneration and cancer, across multiple neuroscience model systems

This package contains functions for:

  + gathering `gene ontology` information for streamlined downstream
  bioinformatics plotting
  + gathering orthologous gene and alias information
  + gathering the genes and specific `gene ontology` terms associated with a
  generic `gene ontology`
  + extracting a Z-score matrix from a transcriptional profile dataframe for
  plotting heatmaps

Three publications involved the use of these functions in this package (`TMEM`) 
and each of these publications have (or will have) their own repositories
([eriklarsen4/Itch](<https://github.com/eriklarsen4/Itch>),
[eriklarsen4/Endo](<https://github.com/eriklarsen4/Endo>), and 
`eriklarsen4/Hippo` coming soon). Please see those repositories for each paper's
respective URL.

Included data for the `TMEM` package is derived from the publicly available
"[itch paper](<https://pmc.ncbi.nlm.nih.gov/articles/PMC8854445>)." All of that 
publication's genomics and Calcium imaging data, analyses, and code are 
contained within the [eriklarsen4/Itch](<https://github.com/eriklarsen4/Itch>)
repository. Please visit that repository for more detail.

These functions were developed, in-house, years before any `Bioconductor`
APIs were developed. The goal was to integrate downstream analytical results (pathway
analysis, gene ontology, networks, etc., from e.g. [PANTHER](<https://geneontology.org>)
or KEGG) to DEA and GSEA results files for mapping-- for visualizing results
without downloading more files, manually transcribing results, or creating extra
plots. I have not yet compared the functionality between `Bioconductor`/GO.db with this package. 

Please see the [vignette](<https://github.com/eriklarsen4/TMEM/blob/main/vignettes/TMEM.md>) for more detail.

Coming soon to this repository are whole-transcriptome sequencing 
(WTS) and whole-exome sequencing (WES) tutorials, along with a Markdown with 
examples for building a `Shiny` app, with a link to the example app.
