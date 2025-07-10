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

These functions were developed, in-house, years before any `Bioconductor`
updates (APIs) that now contain similar functionality.
The goal was to enable direct analysis and plotting of data
acquired from websites providing downstream GSEA analysis (visualizing
[usegalaxy.org](<https://usegalaxy.org>) DESeq2 results,
[geneontology.org/PANTHER](<https://geneontology.org>)/KEGG results, etc.).
I have not yet compared the functionality between `Bioconductor`/GO.db with this package.

Three publications involved the use of these functions. Each of these
publications have (or will have) their own repositories
([eriklarsen4/Itch](<https://github.com/eriklarsen4/Itch>),
[eriklarsen4/Endo](<https://github.com/eriklarsen4/Endo>), and 
`eriklarsen4/Hippo` coming soon).
Please see those repositories for vignettes including transcript expression [heatmaps, differential gene set analyses (volcano plots, MA plots)](https://github.com/eriklarsen4/Itch/blob/main/Code/DEA/ItchDEA.md),  [gene set enrichment analysis results](<https://github.com/eriklarsen4/Itch/blob/main/Code/GSEA/ItchGSEA.md>), [general statistical analyses (hypothesis testing and exploratory data analysis)](https://github.com/eriklarsen4/Endo/blob/master/vignettes/FIREpHly.md), and informal time series analyses. 
They also contain each paper's respective URL.

The data included in this package is derived from the publicly available dataset from the
"[itch paper](<https://pmc.ncbi.nlm.nih.gov/articles/PMC8854445>)". This data is also
contained within the [eriklarsen4/Itch](<https://github.com/eriklarsen4/Itch>)
repository.

Please see this repository's [vignette](<https://github.com/eriklarsen4/TMEM/blob/main/vignettes/TMEM.md>) for more detail about the **TMEM** package's functions.
