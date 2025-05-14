 <!-- badges: start -->
  [![R-CMD-check](https://github.com/eriklarsen4/TMEM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/eriklarsen4/TMEM/actions/workflows/R-CMD-check.yaml)
  [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/TMEM)](https://cran.r-project.org/package=TMEM)
  <!-- badges: end -->

**TMEM**

Repository and package for bioinformatics- and molecular biology-related functions, facilitating biological data exploration and easier plotting of downstream NGS analyses (heatmaps, bargraphs, volcano plots, MA plots, scatter plots)

Adapted from source code used to analyze data in **Dr. Martha Bhattacharya's lab** at the **University of Arizona**

  + The lab studies transmembrane protein 184b (TMEM184B), a protein involved in axon degeneration, across multiple neuroscience model systems

This package contains functions for:

  + gathering `gene ontology` information for streamlined downstream bioinformatics plotting
  + gathering orthologous gene and alias information
  + gathering the genes and specific `gene ontology` terms associated with a generic `gene ontology`
  + extracting a Z-score matrix from a transcriptional profile dataframe for plotting heatmaps

Included data is derived from the publicly available dataset from the "[itch paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8854445/)",
also contained within the [Itch](https://github.com/eriklarsen4/Itch) repository

Three separate publications involved the use of these functions (see [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC8854445/), [here](https://pubmed.ncbi.nlm.nih.gov/37730546/), and [here](https://pubmed.ncbi.nlm.nih.gov/39975166/)); these publications
have their own repositories and packages (again, see [here](https://github.com/eriklarsen4/Itch)).
These functions were developed, in-house, years before any Bioconductor updates (APIs) to enable R/Python users of GSEA tools to access the same data from the command line instead of downloading and parsing from websites.

A vignette has been made available [here](https://github.com/eriklarsen4/TMEM/blob/main/vignettes/TMEM.md)
