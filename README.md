 <!-- badges: start -->
  [![R-CMD-check](https://github.com/eriklarsen4/TMEM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/eriklarsen4/TMEM/actions/workflows/R-CMD-check.yaml)
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

Included data is derived from the publicly available dataset from the "[itch paper](<https://pmc.ncbi.nlm.nih.gov/articles/PMC8854445>)",
also contained within the [Itch](<https://github.com/eriklarsen4/Itch>) repository.

Three separate publications involved the use of these functions (see [here](<https://pmc.ncbi.nlm.nih.gov/articles/PMC8854445>), [here](<https://pubmed.ncbi.nlm.nih.gov/37730546>), and [here](<https://pubmed.ncbi.nlm.nih.gov/39975166>)); these publications
have (or will have) their own repositories.

These functions were developed, in-house, years before any `Bioconductor` updates (APIs).
The goal was to enable direct analysis and plotting of data acquired from websites providing downstream GSEA analysis (visualizing [usegalaxy.org](<https://usegalaxy.org/>) DESeq2 results, [geneontology.org/PANTHER](<https://geneontology.org/>)/KEGG results, etc.).
I have not yet compared the functionality between `Bioconductor`/GO.db with this package.

Please see the [vignette](<https://github.com/eriklarsen4/TMEM/blob/main/vignettes/TMEM.md>) for more detail.
