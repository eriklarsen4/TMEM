## Package Install Options

### **Install the TMEM Package dev version from GitHub**

Install the package from `github`:

``` r 
remotes::install_github('https://github.com/eriklarsen4/TMEM') 
```

### **Install the TMEM Package from CRAN**

Alternatively, if the package is accepted by `CRAN`, install from `CRAN`

``` r
utils::install.packages('TMEM')
```

### **Attach TMEM Package**

Attach the package

``` r
library(TMEM)
```

We'll use the datasets included in the `TMEM` package to demonstrate its
functions

## Functions

### **get_GO_info**

After installing and attaching the `TMEM` package, one of its crucial functions
is the `get_GO_info` that gathers all the gene ontology information of a list of
gene or protein identifiers

This is called **gene set enrichment analysis (GSEA)**

The `GSEA` conducted on this dataset is included in full in
[a markdown of another repository](https://github.com/eriklarsen4/Itch/blob/main/Code/GSEA/ItchGSEA.md),
along with additional analyses for [that paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8854445/#_ad93_)

The real utility of the function is gathering all the relevant information in
one place without needing to generate multiple files and analyses in a browser--
it's all in-house

However, this `GSEA` utilizes `Bioconductor`'s **AnnotationDbi** package, thus,
there are discrepancies with, for example, `geneontology.org`'s annotations and its
statistics (test statistics and p-value), even when computed exactly as described in
`geneontology.org/PANTHER`'s [documentation](https://pantherdb.org/tips/tips_overrep.jsp)

Relatedly, `Bioconductor`'s gene ontologies do not employ `geneontology.org`'s gene ontology (nesting)
structure, which effects any statistics a user might try to extract

![]("C:\Users\Erik\Pictures\Screenshots\GOorgNesting.jpg")<!-- -->

I suggest using this function for exploratory work and plotting, while relying
on `geneontology.org`'s [web browser](https://geneontology.org/) directly for
publications-- I've noticed differences between `GO.db R package` results and
those from the browser.

The statistical discrepancies are noted in the examples section of this vignette

First, import the data in the package

``` r
data("aDRG_DEG_list")
```

Next, pass the variable into the function (data was derived from mouse)

``` r
results <- TMEM::get_GO_info(list_of_interest = aDRG_DEG_list,
                             species = 'mouse')
```

Inspect results by unpacking everything from the returned list

```r
  ## GO INFO about each Unique GO returned by the list
GO_info_by_term_df <- results |> 
  purrr::keep_at("GO_info_by_term_df") |> 
  as.data.frame() %>%
  dplyr::rename_with(.,
                    ~gsub(.x,
                          pattern = 'GO_info_by_term_df\\.',
                          replacement = ''))

  ## GO INFO about each gene in the list
gene_GO_info_df <- results |> 
  purrr::keep_at("gene_GO_info_df") |>
  as.data.frame() %>%
  dplyr::rename_with(.,
                    ~gsub(.x,
                          pattern = 'gene_GO_info_df\\.',
                          replacement = ''))

  ## every alias for each gene in the list
aliases <- results |> 
  purrr::keep_at("aliases")

  ## all aliases from the list without the genes of the list
list_of_interest_aliases <- results |> 
  purrr::keep_at("list_of_interest_aliases") |>
  unlist() |> 
  as.character()

  ## Unique GOs returned by the list
unique_GOs <- results |> 
  purrr::keep_at("unique_GOs")

  ## Unique GO IDs returned by the list
unique_GO_IDs <- results |> 
  purrr::keep_at("unique_GO_IDs")
```

### **get_orthologs_and_aliases**

Pass the gene list to the `get_orthologs_and_aliases`

Since the data is derived from `mouse`, `mouse` will serve as the `ref_species`

```r
ortholog_df <- TMEM::get_orthologs_and_aliases(ref_species = 'mouse', 
                                               list_of_interest = aDRG_DEG_list)
```

Print the `get_orthologs_and_aliases` function output to console

```r
head(ortholog_df)
```
```r
# A tibble: 6 × 7
  ref_species target_species input_gene ortholog_gene input_ensg   ortholog_ensg
  <chr>       <chr>          <chr>      <chr>         <chr>        <chr>        
1 mouse       human          Il31ra     IL31RA        ENSMUSG0000… ENSG00000164…
2 mouse       NA             Il31ra     NA            ENSMUSG0000… NA           
3 mouse       macaque        Il31ra     IL31RA        ENSMUSG0000… ENSMMUG00000…
4 mouse       human          Shtn1      SHTN1         ENSMUSG0000… ENSG00000187…
5 mouse       fly            Shtn1      Frl           ENSMUSG0000… FBgn0267795  
6 mouse       macaque        Shtn1      SHTN1         ENSMUSG0000… ENSMMUG00000…
# ℹ 1 more variable: description <chr>
```

It may be very useful to pass the list of genes with their aliases 
(c(`list_of_interest_aliases`, `aDRG_DEG_list`)) to the `get_GO_info` function

### **query_GO**

Another useful function with multiple uses is finding genes/proteins annotated 
to certain `GO Terms`

This can be useful for:
  + investigating molecular pathways
  + investigating complexes (for cross-checking or discovery)
  + finding multiple `GO Terms` that contain a regular expression string

```r
TMEM::query_GO(model_org = 'human',
               GO_db = GO.db::GO.db,
               string_terms = 'dense core vesicle|lysosome')
```
```r
[1] "CTSB" "CTSK" "CTSL" "CTSS" "LGMN"
```

Extracting the gene/protein IDs associated with the `GO Terms` derived from the
given regular expression(s) is non-trivial, depending on the function call.

Those gene/protein IDs associated with a given term can be extracted directly
using `dplyr` without storing all the function's results to an object in the
global environment: 

```r
TMEM::query_GO(model_org = 'human',
               GO_db = GO.db::GO.db,
               string_terms = 'dense core vesicle|lysosome') |> 
  purrr::keep_at("GO_df") |> 
  as.data.frame() %>%
  dplyr::rename_with(., ~gsub(.x, pattern = 'GO_df.', replacement = '')) %>%
  dplyr::rename_with(., ~gsub(.x, pattern = '\\.\\.', replacement = 'N.')) |>
  dplyr::filter(GO.Term == 'endolysosome membrane') |> 
  dplyr::select(Genes) |> 
  unlist() |> 
  as.character() %>%
  stringr::str_split_1(., pattern = ';')
```
```r
[1] "CTSB" "CTSK" "CTSL" "CTSS" "LGMN"
```

Alternatively, all of the results can be stored to an object in the global
environment from which (in this example) the gene/protein IDs can be extracted:

```r
query_GO_results <- TMEM::query_GO(model_org = 'human',
                                   GO_db = GO.db::GO.db,
                                   string_terms = 'dense core vesicle|lysosome')

query_GO_results$GO_df |> 
  dplyr::filter(`GO Term` == 'endolysosome membrane') |> 
  dplyr::select(Genes) |> 
  unlist() |> 
  as.character() %>%
  stringr::str_split_1(., pattern = ';')
```
```r
[1] "CTSB" "CTSK" "CTSL" "CTSS" "LGMN"
```

### **find_row_Z**

This function computes a Z-score for a provided matrix, currently agnostic to 
grouping

This is highly useful for generating transcriptional heatmaps

```r
data("aDRG_TPM")
head(TMEM::find_row_Z(Expression_Profile = aDRG_TPM))
```
```r
          GeneID          WT1         WT2        WT3        WT4       Mut1
1      -343C11.2  0.499032872 -0.06255327 -0.7815054 -0.5794017 -1.2862768
2 00R_AC107638.2 -1.173563905  0.73833452  0.5446953 -1.1735639 -0.2224270
3      00R_Pgap2 -0.823887096 -0.01608605  1.4992682 -0.8238871  0.4620882
4  0610009O20Rik  0.002055188 -0.59281334 -0.4189258 -0.3383803  1.0251002
5          1-Mar -0.848854933  0.81609116 -0.1023590 -0.4994568  1.2086874
6          1-Sep -0.374463026  0.05814900 -0.5354068 -0.6127807 -0.6432380
       Mut2       Mut3       Mut4
1 0.6287173  1.9089070 -0.3269201
2 1.7951670 -0.2210958 -0.2875463
3 1.3502780 -0.8238871 -0.8238871
4 0.4602527 -1.6633984  1.5261097
5 1.3367426 -0.5812425 -1.3296080
6 2.2909830  0.4158662 -0.5991097
```

## Example Uses

### **Finding genes related by a generic GO Term relevant to a list of interest**

Interested in finding out what has been documented to be related to a certain
cellular phenomenon?

When there are multiple terms involving a certain phrase (i.e. organelles or
processes), simply searching GO terms using regular expressions is an effective
approach

Obviously, this provides larger scope than investigating a specific term

```r
GO_info_by_term_df |> 
  dplyr::filter(grepl(GO_Term, pattern = 'dense core vesicle')) |> 
  dplyr::select(Gene_IDs_from_List) |> 
  unlist() |> 
  as.character() %>%
  paste0(., collapse = ';') %>%
  stringr::str_split_1(., pattern = ';') |> 
  unique()
```
```r
[1] "Sst"    "Syt4"   "Scg2"   "Syt5"   "Unc13c" "Syt6"
```

Similarly, the results returned by the `get_GO_info` function include all the information
on all the `GO Term`s-- not just including the genes/proteins from the list of
interest that are associated with a `GO Term` (these IDs are called
"`Overlap`"); it provides *all* the IDs associated with the `GO Term`

```r
GO_info_by_term_df |> 
  dplyr::filter(grepl(GO_Term, pattern = 'dense core vesicle')) |> 
  dplyr::select(Gene_IDs) |> 
  unlist() |> 
  as.character() %>%
  paste0(., collapse = ';') %>%
  stringr::str_split_1(., pattern = ';') |> 
  unique()
```
```r
 [1] "Adrb1"    "Adrb2"    "App"      "Avp"      "Bdnf"     "Cacna2d1"
 [7] "Calca"    "Chga"     "Dvl1"     "Fzd8"     "Gnai2"    "Htr1d"   
[13] "Igf1"     "Kif1a"    "Npy1r"    "Oprd1"    "Pdyn"     "Penk"    
[19] "Plat"     "Plcb2"    "Scg2"     "Sst"      "Syt4"     "Cadps"   
[25] "Syt5"     "Calcrl"   "Npff"     "Stxbp5"   "Npy"      "Stxbp5l" 
[31] "Slc18a2"  "Grp"      "P2rx2"    "Dmxl2"    "Vps13a"   "P2rx7"   
[37] "Rab3a"    "Snap25"   "Stxbp1"   "Unc13b"   "Syt6"     "Syt10"   
[43] "Rims1"    "Unc13c"   "Unc13a"
```

Lastly, to find those IDs across generic GO Terms **NOT** also in the provided
list of IDs:

```r
"%notin%" <- Negate("%in%")

c(GO_info_by_term_df |> 
    dplyr::filter(grepl(GO_Term, pattern = 'dense core vesicle')) |> 
    dplyr::select(Gene_IDs) |>
    unlist() |> 
    as.character() %>%
    paste0(., collapse = ';') %>%
    stringr::str_split_1(., pattern = ';') |> 
    unique() # All the genes in all the GO Terms
  )[
    which(
  c(GO_info_by_term_df |> 
      dplyr::filter(grepl(GO_Term, pattern = 'dense core vesicle')) |> 
      dplyr::select(Gene_IDs) |>
      unlist() |> 
      as.character() %>%
      paste0(., collapse = ';') %>%
      stringr::str_split_1(., pattern = ';') |> 
      unique()
  ) %notin%
  c(GO_info_by_term_df |>
      dplyr::filter(grepl(GO_Term, pattern = 'dense core vesicle')) |> 
      dplyr::select(Gene_IDs_from_List) |> 
      unlist() |> 
      as.character() %>%
      paste0(., collapse = ';') %>% 
      stringr::str_split_1(., pattern = ';') |> 
      unique() ) == T
)
] # All the genes in all the GO Terms indexed by those that are NOT in the list 
    # of interest
```
```r
 [1] "Adrb1"    "Adrb2"    "App"      "Avp"      "Bdnf"     "Cacna2d1"
 [7] "Calca"    "Chga"     "Dvl1"     "Fzd8"     "Gnai2"    "Htr1d"   
[13] "Igf1"     "Kif1a"    "Npy1r"    "Oprd1"    "Pdyn"     "Penk"    
[19] "Plat"     "Plcb2"    "Cadps"    "Calcrl"   "Npff"     "Stxbp5"  
[25] "Npy"      "Stxbp5l"  "Slc18a2"  "Grp"      "P2rx2"    "Dmxl2"   
[31] "Vps13a"   "P2rx7"    "Rab3a"    "Snap25"   "Stxbp1"   "Unc13b"  
[37] "Syt10"    "Rims1"    "Unc13a"
```

### **Finding genes related by multiple precise GO Terms**

Interested in which genes from your list of interest overlap across multiple
`GO Term`s elicited by your query/list of interest?

```r
c(GO_info_by_term_df |> 
  dplyr::filter(GO_Term == 'neuronal dense core vesicle' |
                  GO_Term == 'growth') |> 
  dplyr::select(Gene_IDs_from_List) |> 
  unlist() |> 
  as.character() %>%
  paste0(., collapse = ';')%>%
  stringr::str_split_1(., pattern = ';'))[which(GO_info_by_term_df |> 
  dplyr::filter(GO_Term == 'neuronal dense core vesicle' |
                  GO_Term == 'growth'
                ) |> 
  dplyr::select(Gene_IDs_from_List) |> 
  unlist() |>
  as.character() %>%
  paste0(., collapse = ';') %>%
  stringr::str_split_1(., pattern = ';') |> 
  duplicated() == T
)]
```
```r
[1] "Syt4"
```

### **Finding GO terms of any number of genes**

When a user is interested in all of the `GO Terms` that are shared by multiple
gene/protein IDs from the user's list of gene/protein IDs:

```r
head(gene_GO_info_df |> 
  dplyr::filter(GeneID == 'Lpar3' |
                  GeneID == 'Lpar5') |>
  dplyr::select(GO_Terms) |> 
  unlist() |> 
  as.character() %>%
  paste0(., collapse = ';') %>%
  stringr::str_split_1(., pattern = ';') |>
  unique()
  )
```
```r
[1] "MAPK cascade"                             
[2] "cell morphogenesis"                       
[3] "regulation of cell growth"                
[4] "G-protein alpha-subunit binding"          
[5] "molecular_function"                       
[6] "transmembrane signaling receptor activity"
```

Similarly, repeating the above for just one gene/protein is simpler-- remove the
`|` and subsequent `GeneID`s:

```r
head(gene_GO_info_df |>
  dplyr::filter(GeneID == 'Il31ra') |> 
  dplyr::select(GO_Terms) |> 
  unlist() |>
  as.character() %>%
  stringr::str_split_1(., pattern = ';')
  )
```
```r
[1] "cytokine production"                              
[2] "columnar/cuboidal epithelial cell differentiation"
[3] "glandular epithelial cell differentiation"        
[4] "adaptive immune response"                         
[5] "immune effector process"                          
[6] "cytokine production involved in immune response"
```

### **Finding genes of a specific GO Term**

```r
head(GO_info_by_term_df |> 
  dplyr::filter(GO_Term == 'neuronal dense core vesicle') |> 
  dplyr::select(Gene_IDs) |> 
  unlist() |> 
  as.character() %>%
  paste0(., collapse = ';') %>%
  stringr::str_split_1(., pattern = ';')
)
```
```r
[1] "Adrb1"    "Adrb2"    "App"      "Avp"      "Bdnf"     "Cacna2d1"
```

### **Finding enriched GOs from get_GO_info**

Determining which `GO`s (themes) are worth looking into further is the entire
point of `GSEA`

The `get_GO_info` function returns all the `GO Terms` relevant to a user's provided list,
and does not include all annotated `GO Terms`

The function also relies on `Bioconductor`'s species-specific databases, as well
as its own version of the `GO.db`'s annotations (which may, *themselves* differ
from annotations on the `geneontology.org` website)

Thus, there are discrepancies in which genes map to which terms, and in the
adaptation of [PANTHER](https://pantherdb.org/tools/)'s
(via [geneontology.org](https://geneontology.org/)) statistical methods

In brief, `PANTHER` permits the computation of a binomial test statistic or, via
Fisher's exact test, the probability of (in this case) observing at least the
number of gene/protein IDs (or more) from the submitted list in a given `GO Term`

`PANTHER` provides corresponding adjusted p-values via the
`Benjamini-Hochberg Method` for each statistic (Fisher's p-value, binomial test
statistic's p-value)

See below examples of this

```r
Genes_in_Human_Genome <- 23481
Genes_in_Mouse_Genome <- 54879

# note that there are 21836 uniquely mapped genes

GO_info_by_term_df2 <- GO_info_by_term_df |> 
  dplyr::arrange(desc(Overlap)) |> 
  dplyr::mutate(Term_Freq = GO_Term_Size/Genes_in_Mouse_Genome,
                Expected = length(aDRG_DEG_list) * Term_Freq,
                FE = Overlap/Expected) |> 
  dplyr::filter(!is.na(GO_Term))
```

Compute the p-values according to `PANTHER`'s documentation

```r
for (i in 1:nrow(GO_info_by_term_df2)) {
  
  # compute p-value for binomial test statistic according to PANTHER's 
  # documentation
  
    if (GO_info_by_term_df2$Overlap[i] > GO_info_by_term_df2$Expected[i]) { # for over-representation
    GO_INFO_by_TERM_df2$pval[i] <- sum(
      (GO_info_by_term_df2$Term_Freq[i]^(seq(GO_info_by_term_df2$Overlap[i],
                                             length(aDRG_DEG_list))
                                         ))*(1 - GO_info_by_term_df2$Term_Freq[i])^(seq(length(aDRG_DEG_list)-GO_info_by_term_df2$Overlap[i],
                                                                                      0
                                                                                     )
                                                                                 )
    )
    } else { # for under-representation
    GO_info_by_term_df2$pval[i] <- sum(
      ((GO_info_by_term_df2$Term_Freq[i])^(seq(0,GO_info_by_term_df2$Overlap[i])
                                         ))*(1 - GO_info_by_term_df2$Term_Freq[i])^(seq(length(aDRG_DEG_list),
                                                                                      GO_info_by_term_df2$Overlap[i]
                                                                                     
                                                                                     )
                                                                                 )
    )
    }
}
```

Compute p-values using `rstatix`'s binomial test for both `df`s

```r

GO_info_by_term_df2$binom_pval = NA_real_

for (i in 1:nrow(GO_info_by_term_df2)) {
  GO_info_by_term_df2$binom_pval[i] = rstatix::binom_test(
    x = c(GO_info_by_term_df2$Overlap[i],
          GO_info_by_term_df2$GO_Term_Size[i]-GO_info_by_term_df2$Overlap[i]),
    p = (GO_info_by_term_df2$GO_Term_Size[i]/Genes_in_Mouse_Genome),
    alternative = 'two.sided',
    conf.level = 0.95,
    detailed = TRUE) |> 
    dplyr::select(p) |> unlist() |> as.character() |> as.numeric()
}

```

Add Adjusted P-values using BHM

Note these, and all p-values do not match those from all the same results from
the `geneontology.org` website

```r

GO_info_by_term_df2$binom_adjp = NA_real_

for (i in 1:nrow(GO_info_by_term_df2)) {
  GO_info_by_term_df2$binom_adjp[i] = rstatix::binom_test(
    x = c(GO_info_by_term_df2$Overlap[i],
          GO_info_by_term_df2$GO_Term_Size[i]-GO_info_by_term_df2$Overlap[i]),
    p = (GO_info_by_term_df2$GO_Term_Size[i]/Genes_in_Mouse_Genome),
    alternative = 'two.sided',
    conf.level = 0.95,
    detailed = TRUE) |> 
    dplyr::select(p) |> unlist() |> as.character() |> as.numeric()
}

GO_info_by_term_df2 <- GO_info_by_term_df2 |> 
  dplyr::mutate(adjp = p.adjust(p = GO_info_by_term_df2$pval,
                                method = 'BH'), .after= pval)

GO_info_by_term_df2$binom_adjp = p.adjust(p = GO_info_by_term_df2$binom_pval,
                                          method = 'BH')
```
