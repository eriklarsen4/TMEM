#'
#' @docType package
#' @name TMEM
#' @format This package includes three datasets to ensure function examples are testable:
#' \describe{
#'    \item{aDRG_DEA_results}{a dataframe of the DESeq2 Wald test output of 2 groups (Wild-Type and Mutant) of 4 samples of bulk RNA derived from mouse dorsal root ganglion neurons}
#'    \item{aDRG_TPM}{a dataframe of DESeq2's computed normalized counts across 2 groups (Wild-Type and Mutant) of 4 samples of bulk RNA derived from mouse dorsal root ganglion neurons}
#'    \item{aDRG_DEG_list}{a character vector containing the mm9 Entrez GeneIDs of the genes DESeq2 identified as differentially expressed}
#' }
#' @source [TMEM184B is necessary for IL-31-induced itch](https://pmc.ncbi.nlm.nih.gov/articles/PMC8854445/#SD5)
#'
#' @title get_GO_info
#'
#' @description
#' \strong{`get_GO_info`} is a function that extracts gene ontology (GO)
#' information about a list (character vector) of human/mouse/drosophila gene
#' symbols
#'
#' @param list_of_interest a character vector of the gene symbols of a
#' list of genes or proteins
#' @param species string of the species of the provided character vector
#'
#' @returns
#' \strong{`GO_info_list`}: a list containing all the objects described below
#'
#' \strong{`gene_GO_info_df`}: a dataframe containing all the genes from the
#' provided list and all the gene ontology terms with which each gene/protein
#' has been associated with, imperically and/or theoretically
#' \strong{`GO_info_by_term_df`}: a dataframe containing all the GO information for
#' each GO term (i.e. how many and which genes are in each term) that is associated
#' with at least one gene symbol from the provided character vector
#' \strong{`unique_GOs`}: a character vector containing all the unique GO terms
#' associated with at least one gene symbol from the provided character vector
#' \strong{`unique_GO_IDs`}: a character vector containing all the unique GO terms'
#' IDs for the terms associated with at least one gene symbol from the provided
#' character vector
#' \strong{`list_of_interest_aliases`}: a character vector of the same-species
#' aliases associated with any gene symbol from the provided character vector AND
#' the provided gene symbol character vector
#' \strong{`aliases`}: a character vector of only the same-species aliases associated
#' with any gene symbol from the provided character vector
#'
#' @details
#' Queries the Bioconductor db to extract the GO info for each gene/protein in a
#' provided character vector.
#'
#' \strong{`species`} can be one of 'human', 'HS', 'homo sapiens', 'mouse', 'Mm', 'mus musculus', 'drosophila', 'fly', 'fruit fly', 'DM' (all case insensitive)
#'
#' @examples
#' data(aDRG_DEG_list)
#' \donttest{
#' get_GO_info(list_of_interest = aDRG_DEG_list, species = 'mouse')
#' }
#'
#' @import AnnotationDbi
#' @import stringr
#' @import dplyr
#' @import GO.db
#' @importFrom rlang .data
#'
#' @rdname get_go_info
#' @export
get_GO_info <- function(list_of_interest, species) {
  # First, define the species of interest ----
  if ( grepl(species, pattern = 'human|HS|homo sapiens', ignore.case = T) ) {
    abbrev_species_name <- org.Hs.eg.db
  } else if ( grepl(species, pattern = 'mouse|Mm|MM|Mus musculus', ignore.case = T) ) {
    abbrev_species_name <- org.Mm.eg.db
  } else if ( grepl(species,
                    pattern = 'drosophila melanogaster|drosophila|fly|DM|fruit fly', ignore.case = T) ) {
    abbrev_species_name <- org.Dm.eg.db
  }
  # Next, acquire all aliases within species ----
  aliases <- list()
  for (i in 1:length(list_of_interest)) {

    if (list_of_interest[i] %in% AnnotationDbi::keys(abbrev_species_name, keytype = 'SYMBOL')) {

      aliases[[i]] <- AnnotationDbi::mapIds(abbrev_species_name,
                                            keys = c(list_of_interest[i]),
                                            keytype = 'SYMBOL',
                                            column = 'ALIAS',
                                            multiVals = 'list')

    } else {

      next()

    }
  }

  ## make a new list of interest that includes the original IDs AND aliases
  list_of_interest_and_aliases <- aliases |> unlist() |> as.character()

  # Find Entrez (Ensembl) IDs from the input list ----
  Ensemble_IDs <- NULL
  Ensembl_IDs <- as.integer(
    AnnotationDbi::mapIds(abbrev_species_name,
                          keys = as.character(
                            c(list_of_interest[])
                          ),
                          keytype = 'SYMBOL',
                          column = 'ENTREZID')
  )

  # Extract GO Terms for all genes/proteins in the input list ----
  GO_IDs <- NULL
  GO_IDs <- c(
    AnnotationDbi::mapIds(abbrev_species_name,
                          keys = as.character(
                            c(list_of_interest[])
                          ),
                          keytype = 'SYMBOL',
                          column = 'GOALL',
                          multiVals = 'list')
  )

  # Re-query the GO IDs having filtered the list ----
  GO_IDs <- c(
      AnnotationDbi::mapIds(abbrev_species_name,
                            keys = as.character(
                              c(list_of_interest[])
                            ),
                            keytype = 'SYMBOL',
                            column = 'GOALL',
                            multiVals = 'list')
  )

  # GO Term dataframe storage for the genes in the list ----
  ## Create a matrix to house each gene/protein's GO information
  GENE_GO_INFO <- matrix(nrow = length(c(list_of_interest)), ncol = 4)
  # rownames(GENE_GO_INFO) = c(list_of_interest)
  colnames(GENE_GO_INFO) <- c("GeneID", "Num_GO_Terms", "GO_Term_IDs", "GO_Terms")

  ## Fill the GO Term ID column with the concatenated GO Term IDs associated with each gene/protein
  ## skip over the genes/proteins without any
  for (i in 1:length(GO_IDs)) {
    if ( !is.na(GO_IDs[[i]]) |> any() ) {

      GENE_GO_INFO[i,3] <- paste0(names(AnnotationDbi::mapIds(GO.db, GO_IDs[[i]][], 'TERM', 'GOID')), collapse = ';')

    } else {

      GENE_GO_INFO[i,3] <- 'None'

    }
  }

  for (i in 1:length(list_of_interest)) {

    GENE_GO_INFO[i,1] <- list_of_interest[i]

  }

  ## Extract the GO IDs and remove the redundant ones;
  ## Fill the 3rd column with the numbers of GO Terms associated with each gene/protein
  ## skip over the genes/proteins without any
  for (i in 1:length(GO_IDs)) {

    if ( !is.na(GO_IDs[[i]]) |> any() ) {

      x <- stringr::str_split(as.vector(GENE_GO_INFO[i,3]), pattern = ';', simplify = TRUE)
      x <- c(as.character(x[which(x %in% x == TRUE )]))
      GENE_GO_INFO[i,3] <- paste0(unique(x), collapse = ';')
      GENE_GO_INFO[i,2] <- paste0(as.numeric(length(unique(x))))

    } else {

      GENE_GO_INFO[i,3] <- 'None'
      GENE_GO_INFO[i,2] <- 0
    }


  }

  ## Fill the 4th column with the concatenated GO Terms associated with each gene/protein
  ## skip over the genes/proteins without any
  for (i in 1:length(GO_IDs)) {
    if ( !is.na(GO_IDs[[i]]) |> any() ) {

      x <- str_split(as.vector(GENE_GO_INFO[i,3]), pattern = ';', simplify = TRUE)
      GENE_GO_INFO[i,4] <- paste0(as.character(AnnotationDbi::mapIds(GO.db, keys = x, keytype = 'GOID', 'TERM')), collapse = ';')
    } else {
      GENE_GO_INFO[i,4] <- 'None'
    }
  }


  ## Extract all of the terms/IDs into one vector to find unique GO terms
  Unique_GOs <- c(unique(as.character(stringr::str_split(as.vector(GENE_GO_INFO[,4]), pattern = ';', simplify = T))))
  Unique_GO_IDs <- c(unique(as.character(stringr::str_split(as.vector(GENE_GO_INFO[,3]), pattern = ';', simplify = T))))

  if (any(Unique_GOs == '')) {

    Unique_GOs <- Unique_GOs[-which(Unique_GOs == '')]
    Unique_GO_IDs <- Unique_GO_IDs[-which(Unique_GO_IDs == '')]
  }

  # GO storage for all GOs ----
  # Find all the genes/proteins within the list of terms
  GO.list <- list()
  GO.list <- AnnotationDbi::mapIds(abbrev_species_name,
                                  keys = Unique_GO_IDs[],
                                  keytype = 'GOALL',
                                  column = 'SYMBOL',
                                  multiVals = 'list')

  GO_INFO <- matrix(nrow = length(Unique_GOs), ncol = 6)
  colnames(GO_INFO) <- c("GO_Term_ID", "GO_Term", "GO_Term_Size", "Gene_IDs", "Overlap", "Gene_IDs_from_List")

  ## Initialize columns
  GO_INFO[ ,c(1,2,4,6)] <- ""
  GO_INFO[ , c(3,5)] <- 0

  for (i in 1:length(Unique_GOs)) {
    ## Store the IDs
    GO_INFO[i,1] <- Unique_GO_IDs[i]
    GO_INFO[i,2] <- Unique_GOs[i]

    ## Number of genes in each term
    GO_INFO[i,3] <- as.numeric(length(unique(GO.list[[i]])))
    ## Genes in each term
    GO_INFO[i,4] <- paste0(unique(c(GO.list[[i]])), collapse = ';')
    ## Number of genes in each term also in the provided list
    GO_INFO[i,5] <- as.numeric(length(which(c(list_of_interest) %in% unique(GO.list[[i]]) == T)))
    ## Genes in each term also in the provided list
    GO_INFO[i,6] <- paste0(c(list_of_interest[which(c(list_of_interest) %in% unique(GO.list[[i]]) == T)]), collapse = ';')
  }
  ## Convert to df
  GO_INFO <- as.data.frame(GO_INFO) |> dplyr::mutate(GO_Term_Size = as.numeric(.data$GO_Term_Size), Overlap = as.numeric(.data$Overlap))
  GENE_GO_INFO_df <- as.data.frame(GENE_GO_INFO) |>
    dplyr::mutate(Num_GO_Terms = as.numeric(.data$Num_GO_Terms)) |>
    dplyr::filter(!is.na(.data$GeneID))


  # Export ----
  GO_info_list <- list("gene_GO_info_df" = GENE_GO_INFO_df,
                       "aliases" = aliases,
                       "list_of_interest_aliases" = list_of_interest_and_aliases,
                       "unique_GOs" = Unique_GOs,
                       "unique_GO_IDs" = Unique_GO_IDs,
                       "GO_info_by_term_df" = GO_INFO)
  # GENE_GO_INFO_df <<- GENE_GO_INFO_df
  # GO_INFO_by_TERM_df <<- GO_INFO
  # Unique_GOs <<- Unique_GOs
  # Unique_GO_IDs <<- Unique_GO_IDs
  # list_of_interest_aliases <<- list_of_interest_and_aliases
  # # all_unique_genes <<- all_unique_genes
  # aliases <<- aliases
  return(GO_info_list)
}

#' @title get_orthologs_and_aliases
#'
#' @description
#' \strong{`get_orthologs_and_aliases`} is a function that finds all the aliases and
#' orthologs of the other four main model organsims from a character vector of
#' gene symbols of another
#'
#' @param ref_species a string containing the nickname of a model species
#' @param list_of_interest a character vector containing gene symbols
#'
#' @returns
#' \strong{`result`}: a dataframe housing:
#'  + the reference species
#'  + the reference gene symbols
#'  + the reference gene symbols' Ensembl IDs
#'  + remaining (target) species names
#'  + target species' orthologous gene symbols
#'  + target species' orthologous gene Ensembl IDs
#'
#' @details
#' Queries the Bioconductor db to extract the aliases of the provided gene symbols,
#' then uses the orthogene package to extract orthologous gene information from
#' the 5 main biological model species, and corrals into a practical and intuitive dataframe
#'
#' \strong{`ref_species`} can be one of 'human', 'mouse', or 'fly'
#'
#' @examples
#' \donttest{data("aDRG_DEG_list")
#' get_orthologs_and_aliases(ref_species = 'mouse',
#'                          list_of_interest = aDRG_DEG_list)
#' }
#'
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Dm.eg.db
#' @import AnnotationDbi
#' @import dplyr
#' @import purrr
#' @import orthogene
#' @import tidyr
#' @importFrom rlang .data
#'
#' @references [orthogene](https://github.com/neurogenomics/orthogene)
#'
#' @export
get_orthologs_and_aliases <- function(ref_species, list_of_interest) {
  # initialize variables ----
  ## list for aliases in the reference species
  listy <- list()

  ## the list for the results
  result <- list()

  "%notin%" <- Negate("%in%")

  ## all the model organisms
  species <- c('human', 'mouse', 'fly', 'macaque', 'zebrafish')
  species <- species[which(species %notin% ref_species == T)]

  ## the reference species db object
  species_nickname_db_object <- vector()
  if (grepl(ref_species, pattern = 'human')) {
    species_nickname_db_object <- org.Hs.eg.db
  } else if (grepl(ref_species, pattern = 'mouse')) {
    species_nickname_db_object <- org.Mm.eg.db
  } else if (grepl(ref_species, pattern = 'fly')) {
    species_nickname_db_object <- org.Dm.eg.db
  }

  # extract the list of aliases ----
  for (i in 1:length(list_of_interest)) {
    if ( list_of_interest[i] %in% AnnotationDbi::keys(species_nickname_db_object, keytype = 'SYMBOL') ) {

      listy[[i]] <- AnnotationDbi::mapIds(species_nickname_db_object,
                                          keys = c(list_of_interest[i]),
                                          keytype = 'SYMBOL',
                                          column = 'ALIAS',
                                          multiVals = 'list')

    } else {

      next()

    }
  }

  ## convert the aliases list to a dataframe
  ALIAS <- purrr::list_flatten(listy) |>
    purrr::map_df(.f = as.data.frame) |>
    dplyr::rename(Aliases = .data$`.x[[i]]`) |>
    dplyr::mutate(SYMBOL = '', .before = 1)

  ## fill the "Symbol" column with the genes from the list of interest
  ## aliases of these genes are the 'Aliases' column
  for (i in 1:length(listy)) {

    ALIAS$SYMBOL[which(ALIAS$Aliases %in% listy[[i]][[1]] == T)] = listy[[i]] |> names()

  }

  ## convert "N/A" to NA
  ALIAS <- ALIAS |>
    dplyr::group_by(.data$SYMBOL) |>
    dplyr::mutate(Aliases = case_when(length(.data$Aliases) == 1 ~ NA, .data$SYMBOL == .data$Aliases ~ '', TRUE ~ .data$Aliases)) |>
    dplyr::ungroup() |>
    dplyr::filter(is.na(.data$Aliases) | .data$Aliases != '')

  # extract ortholog information on all genes (genes in the input list and aliases) ----
  for (i in 1:length(species) ) {

    result[[i]] <- orthogene::map_orthologs(genes = c(ALIAS |>
                                                        unlist() |>
                                                        as.character() |>
                                                        unique())[which(
                                                          !is.na(c(ALIAS |>
                                                                     unlist() |>
                                                                     as.character() |>
                                                                     unique())) == T)],
                                            input_species = ref_species,
                                            output_species = species[i]) |>
      dplyr::mutate(ref_species = as.character(ref_species), .before = 1)

  }

  result <- purrr::list_merge(result) %>%
    purrr::map_df(.f = as.data.frame) |>
    dplyr::mutate(input_number = as.character(.data$input_number)) |>
    tidyr::pivot_longer(cols = c(.data$input_number,
                                 .data$input_ensg,
                                 .data$ortholog_gene,
                                 .data$ortholog_ensg,
                                 .data$description)) %>%
    dplyr::mutate(across(c(1:5), ~case_when(.x == 'N/A' ~ NA, TRUE ~ .x))) |>
    tidyr::pivot_wider(id_cols = c(.data$ref_species,
                                   .data$input_gene),
                       names_from = .data$name,
                       values_from = .data$value,
                       values_fn = list) |>
    tidyr::unnest(cols = everything()) |>
    dplyr::distinct() |>
    dplyr::mutate(target_species = case_when(grepl(.data$ortholog_ensg,
                                                   pattern = 'FB') ~ 'fly',
                                             grepl(.data$ortholog_ensg,
                                                   pattern = 'ENSMM') ~ 'macaque',
                                             grepl(.data$ortholog_ensg,
                                                   pattern = 'ENSG') ~ 'human',
                                             grepl(.data$ortholog_ensg,
                                                   pattern = 'ENSD') ~ 'zebrafish',
                                             grepl(.data$ortholog_ensg,
                                                   pattern = 'ENSMU') ~ 'mouse'),
                  .after = 1) |>
    dplyr::select(-.data$input_number) |>
    dplyr::relocate(.data$ortholog_gene, .after = .data$input_gene)

  # output ----
  # result <<- result
  return(result)
}

#' @title query_GO
#' @description
#' \strong{`query_GO`} is a function that queries the GO.db with a string
#' term of interest for a given model organism
#'
#' @param model_org a string comprising one of: (for human)
#'  'human', 'HS', 'homo sapiens', (for mouse) 'mouse', 'MM', 'mus musculus',
#'  (for fly) drosophila', DM', 'fly'
#'
#' @param string_terms any GO Term string
#'
#' @returns
#' \strong{`query_GO_list`}: a list containing the following objects:
#' \strong{`GO.Terms`}: a character vector of all the GO Terms associated with the organism and string of interest
#' \strong{`GO.IDs`}: a character vector of all the GO IDs of the GO Terms associated with the organism and string of interest
#' \strong{`all_unique_genes`}: a character vector of all the gene symbols associated with provided string
#' \strong{`GO.list`}: a character vector of gene symbols associated with the GO string that is queried
#' \strong{`GO_df`}: a dataframe housing the GO Term, ID, GO size, and its associated gene symbols
#' \strong{`aliases`}: a character vector of the gene symbol aliases associated with the provided GO string
#'
#' @details
#' Queries the GO.db within a species across (potentially) multiple GO Terms to
#' aggregate all the terms and associated genes that overlap with a provided string.
#' Designed for exploration of gene symbols' associations with multiple GO Terms.
#'
#' For example, parsing which genes are associated with multiple MTOR terms (complexes, etc.)
#'
#' @examples
#' \donttest{
#' query_GO(
#'          model_org = 'human',
#'          string_terms = 'dense core vesicle'
#'          )
#' }
#'
#'
#' @import GO.db
#' @import AnnotationDbi
#' @import stringr
#'
#' @export
query_GO <- function(model_org, string_terms) {

  # GO_db<- eval(parse(text = GO.db::GO.db))
  abbrev_species_name <- vector()

  if (grepl(model_org, pattern = 'human|HS|homo sapiens', ignore.case = T)) {

    abbrev_species_name <- org.Hs.eg.db

  } else if (grepl(model_org, pattern = 'mouse|MM|mus musculus', ignore.case = T)) {

    abbrev_species_name <- org.Mm.eg.db

  } else if (grepl(model_org, pattern = 'drosophila|DM|fly|fruit fly', ignore.case = T)) {

    abbrev_species_name <- org.Dm.eg.db

  }

  ## first check that there are GO Terms associated with a provided string
  keys_var <- c(AnnotationDbi::keys(GO.db::GO.db, keytype = 'TERM')[
    which(
      grepl(AnnotationDbi::keys(GO.db::GO.db, keytype = 'TERM'), pattern = string_terms) == TRUE
    )
  ])

  if (length(keys_var) > 0 ) {

    GO.Terms = c(AnnotationDbi::keys(GO.db::GO.db, keytype = 'TERM')[
      which(
        grepl(AnnotationDbi::keys(GO.db::GO.db, keytype = 'TERM'), pattern = string_terms) == TRUE
      )
    ])

    GO.IDs = c(AnnotationDbi::keys(GO.db::GO.db, keytype = 'GOID')[
      which(
        grepl(AnnotationDbi::keys(GO.db::GO.db, keytype = 'TERM'), pattern = string_terms) == TRUE
      )
    ])

    GO.list = c(AnnotationDbi::mapIds(abbrev_species_name,
                       keys = c(as.character(
                         AnnotationDbi::mapIds(GO.db::GO.db,
                                keys = c(
                                  AnnotationDbi::keys(GO.db::GO.db, keytype = 'TERM')[
                                    which(
                                      grepl(AnnotationDbi::keys(GO.db::GO.db, keytype = 'TERM'), pattern = string_terms) == TRUE
                                    )
                                  ]
                                ),
                                keytype = 'TERM',
                                column = 'GOID',
                                multiVals = 'list')
                       )),
                       keytype = 'GOALL',
                       column = 'SYMBOL',
                       multiVals = 'list'))

    rm_idx = c(which(as.character(is.na(GO.list[][])) == "TRUE"))

    if (length(rm_idx) > 0) {

      GO.list = GO.list[-rm_idx]
      GO.terms = GO.Terms[-rm_idx]
      GO.IDs = GO.IDs[-rm_idx]

    }

    df = matrix(nrow = length(c(GO.Terms)), ncol = 4)
    colnames(df) = c("GO Term ID", "GO Term", "# of Genes", "Genes")

    df[ , c(1,2,4)] = ''
    df[ , 3] = 0

    for (i in 1:length(GO.Terms)) {

      df[i,1] = GO.IDs[i]
      df[i,2] = GO.Terms[i]
      x = c(as.character(unlist(as.vector(GO.list[i][]))))
      df[i,3] = paste(as.numeric(length(unique(x))))
      df[i,4] = paste(unique(x), collapse = ';')

    }
    df <- df |> as.data.frame()

    all_unique_genes =
      unique(
        as.character(
          unlist(
            AnnotationDbi::mapIds(abbrev_species_name,
                   keys = c(as.character(
                     AnnotationDbi::mapIds(GO.db::GO.db,
                            keys = c(AnnotationDbi::keys(GO.db::GO.db, keytype = 'TERM'))[
                              which(grepl(AnnotationDbi::keys(GO.db::GO.db, keytype = 'TERM'), pattern = string_terms) == TRUE)
                            ],
                            keytype = 'TERM',
                            column = 'GOID',
                            multiVals = 'list'))),
                     keytype = 'GOALL',
                     column = 'SYMBOL',
                     multiVals = 'list')
            )
          )
        )
    aliases = unique(
      as.character(
        unlist(
          AnnotationDbi::mapIds(abbrev_species_name,
                 keys = c(all_unique_genes),
                 keytype = 'SYMBOL',
                 column = 'ALIAS',
                 multiVals = 'list')
        )
      )
    )
    # print(head(all_unique_genes))
    # print(c("GO Terms = 'GO.Terms' in Global Env.",
    #         "List of genes in all GO Terms = 'all_unique_genes' in Global Env.",
    #         "Dataframe housing all info = 'GO_df' in Global Env.",
    #         "List of list housing raw search results = 'GO.list' in Global Env.",
    #         "aliases = 'aliases' in Global Env."))
    Query_GO_list <- list("GO.Terms" = GO.Terms,
                          "GO.IDs" = GO.IDs,
                          "all_unique_genes" = all_unique_genes[!is.na(all_unique_genes)],
                          "GO.list" = GO.list,
                          "GO_df" = df,
                          "aliases" = aliases)
    # GO.Terms <<- GO.Terms
    # GO.IDs <<- GO.IDs
    # all_unique_genes = all_unique_genes[!is.na(all_unique_genes)]
    # all_unique_genes <<- all_unique_genes
    # GO.list <<- GO.list
    # GO_df <<- df
    # aliases <<- aliases
    return(Query_GO_list)

  } else {

    GO.Terms = 'none'
    GO.IDs = 'none'
    all_unique_genes = 'none'
    GO.list = 'none'
    aliases = 'none'

    if (grepl(string_terms, pattern = '\\|') == TRUE ) {

      string_terms = c(str_split(as.vector(string_terms), pattern = '\\|', simplify = T))

      df = matrix(nrow = length(string_terms), ncol = 4)
      colnames(df) = c("GO Term ID", "GO Term", "# of Genes", "Genes")

      for (i in 1:length(string_terms)) {

        df[i,2] = string_terms[i]

      }
      df[,3] = 0
      df[,c(1,4)] = ""
      df <- df %>% as.data.frame()

    } else {

      df = matrix(nrow = 1, ncol = 4)
      colnames(df) = c("GO Term ID", "GO Term", "# of Genes", "Genes")
      df <- df |> as.data.frame()
      df$`GO Term` = string_terms
      df$`GO Term ID` = ''
      df$`# of Genes` = 0
      df$Genes = ''

    }

    # print(c("GO Terms = 'GO.Terms' in Global Env.",
    #         "List of genes in all GO Terms = 'all_unique_genes' in Global Env.",
    #         "Dataframe housing all info = 'GO_df' in Global Env.",
    #         "List of list housing raw search results = 'GO.list' in Global Env.",
    #         "aliases = 'aliases' in Global Env."))
    query_GO_list <- list("GO.Terms" = GO.Terms,
                          "GO.IDs" = GO.IDs,
                          "all_unique_genes" = all_unique_genes,
                          "GO.list" = GO.list,
                          "GO_df" = df,
                          "aliases" = aliases)
    # GO.Terms <<- GO.Terms
    # GO.IDs <<- GO.IDs
    # all_unique_genes <<- all_unique_genes
    # GO.list <<- GO.list
    # GO_df <<- GO_df
    # aliases <<- aliases
    return(Query_GO_list)
  }

}

#' @title find_row_Z
#' @description
#' \strong{`find_row_Z`} is a function that determines Z-scores of an expression matrix
#'
#' @param Expression_Profile a dataframe containing a list of gene symbols (first column);
#' remaining columns comprising samples (iterations) of expression values
#'
#' @returns
#' \strong{`GeneZ`}: a dataframe with Z-scored expression values
#'
#' @details
#' Z-scores are \strong{`row-wise`} across all columns containing samples (column 2 through the end of the dataframe)
#'
#'
#' @examples
#' \donttest{
#' data(aDRG_TPM)
#' find_row_Z(Expression_Profile = aDRG_TPM)
#' }
#'
#' @import stats
#'
#' @export
find_row_Z <- function(Expression_Profile){
  ## Create dataframes to calculate Z-scores across replicates of each gene.
  ## These will hold the original data frame values until filled in later commands
  ## "MeansAndSDs" will store the means and sds of each gene
  MeansAndSDs = Expression_Profile
  ## "GeneMeans" will store only the means by genes for all replicates
  GeneMeans = Expression_Profile
  ## "GeneSDs" will store only the sds by gene for all replicates
  GeneSDs = Expression_Profile
  ## "GeneZ" will store only the Z-scores, which will be used directly to create the heatmaps
  GeneZ = Expression_Profile

  ## Create a new column to fill with the means and standard deviations of each gene symbol's expression value mean&sd
  MeansAndSDs$mean = 0
  MeansAndSDs$sd = 0

  ## Loop through the dataframe and fill the "mean" and "sd" columns with their appropriate values
  for (i in 1:nrow(MeansAndSDs)){
    MeansAndSDs$mean[i] <- (sum(Expression_Profile[i,c(2:ncol(Expression_Profile))])/ncol(Expression_Profile[,c(2:ncol(Expression_Profile))]))
  }
  for (j in 1:nrow(MeansAndSDs)){
    MeansAndSDs$sd[j] <- sd(Expression_Profile[j,c(2:ncol(Expression_Profile))], na.rm = TRUE)
  }
  ## Create a dataframe storing the gene-specific mean normalized TPM in all columns/replicates for Z-score calculating
  GeneMeans[,c(2:ncol(Expression_Profile))] <- MeansAndSDs$mean
  ## Create a dataframe storing the gene-specific normalized TPM standard deviationsin all columns/replicates for Z-score calculating
  GeneSDs[,c(2:ncol(Expression_Profile))] <- MeansAndSDs$sd
  ## Create the Z-score dataframe
  GeneZ[,c(2:ncol(Expression_Profile))] <- (Expression_Profile[,c(2:ncol(Expression_Profile))] - GeneMeans[,c(2:ncol(Expression_Profile))])/GeneSDs[,c(2:ncol(Expression_Profile))]
  ## Remove the genes that were detected to have 0 TPMs across all samples
  ## Create a function that evaluates a vector/dataframe's (x's) numerical values. Returns equal length vector with T/F bools.
  ## Subset the dataframe that filters those genes.
  row_has_na = apply(GeneZ, 1, function(x){any(is.na(x))})
  GeneZ = GeneZ[!row_has_na,]
  return(GeneZ)
}
