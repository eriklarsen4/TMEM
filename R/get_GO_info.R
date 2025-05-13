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
#' @import assertthat
#' @importFrom AnnotationDbi, keys
#' @importFrom AnnotationDbi, mapIds
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Dm.eg.db
#' @import stringr
#' @import dplyr
#' @import GO.db
#' @importFrom rlang .data
#'
#' @rdname get_GO_info
#' @export
get_GO_info <- function(list_of_interest, species) {
  assertthat::assert_that(!is.null(list_of_interest))
  assertthat::assert_that(!is.null(species))
  assertthat::is.string(species)
  assertthat::is.string(list_of_interest)
  assertthat::not_empty(species)
  assertthat::not_empty(list_of_interest)
  species <- match.arg(species, c("human", "Hs", "HS", "homo sapiens",
                                  "mouse", "Mm", "MM", "Mus musculus",
                                  "fly", "fruit fly", "Dm", "DM", "drosophila", "drosophila melanogaster"))
  # First, define the species of interest ----
  if ( grepl(species, pattern = 'human|HS|homo sapiens', ignore.case = T) ) {
    abbrev_species_name <- org.Hs.eg.db
  } else if ( grepl(species, pattern = 'mouse|Mm|MM|Mus musculus', ignore.case = T) ) {
    abbrev_species_name <- org.Mm.eg.db
  } else if ( grepl(species,
                    pattern = 'drosophila melanogaster|drosophila|fly|DM|fruit fly', ignore.case = T) ) {
    abbrev_species_name <- org.Dm.eg.db
  } else {
    message("invalid 'species' input!")
  }
  # Next, acquire all aliases within species ----
  aliases <- list()
  for (i in 1:length(list_of_interest)) {

    if (list_of_interest[i] %in% AnnotationDbi::keys(abbrev_species_name, keytype = 'SYMBOL')) {

      aliases[[i]] <- suppressMessages(AnnotationDbi::mapIds(abbrev_species_name,
                                                             keys = c(list_of_interest[i]),
                                                             keytype = 'SYMBOL',
                                                             column = 'ALIAS',
                                                             multiVals = 'list')
      )

    } else {

      next()

    }
  }
  ## list returns all input gene/protein IDs AND alias gene/protein IDs
    ## if none returned, input gene/protein IDs are invalid, remaining function fails
  assertthat::see_if(length(aliases |>
                              unlist() |>
                              as.character()) > 0,
                     msg = 'No valid gene symbols!')

  ## make a new list of interest that includes the original IDs AND aliases
  list_of_interest_and_aliases <- aliases |> unlist() |> as.character()

  # Find Entrez (Ensembl) IDs from the input list ----
  Ensemble_IDs <- NULL
  Ensembl_IDs <- as.integer(
    suppressMessages(AnnotationDbi::mapIds(abbrev_species_name,
                                           keys = as.character(
                                             c(list_of_interest[])
                                           ),
                                           keytype = 'SYMBOL',
                                           column = 'ENTREZID')
    )
  )

  # Extract GO Terms for all genes/proteins in the input list ----
  GO_IDs <- NULL
  GO_IDs <- c(
    suppressMessages(AnnotationDbi::mapIds(abbrev_species_name,
                                           keys = as.character(
                                             c(list_of_interest[])
                                           ),
                                           keytype = 'SYMBOL',
                                           column = 'GOALL',
                                           multiVals = 'list')
    )
  )

  # Re-query the GO IDs having filtered the list ----
  GO_IDs <- c(
    suppressMessages(AnnotationDbi::mapIds(abbrev_species_name,
                                           keys = as.character(
                                             c(list_of_interest[])
                                           ),
                                           keytype = 'SYMBOL',
                                           column = 'GOALL',
                                           multiVals = 'list')
    )
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

      GENE_GO_INFO[i,3] <- paste0(names(suppressMessages(AnnotationDbi::mapIds(GO.db::GO.db, GO_IDs[[i]][], 'TERM', 'GOID'))), collapse = ';')

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

      x <- stringr::str_split(as.vector(GENE_GO_INFO[i,3]), pattern = ';', simplify = TRUE)
      GENE_GO_INFO[i,4] <- paste0(as.character(suppressMessages(AnnotationDbi::mapIds(GO.db::GO.db, keys = x, keytype = 'GOID', 'TERM'))), collapse = ';')
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
  GO.list <- suppressMessages(AnnotationDbi::mapIds(abbrev_species_name,
                                                    keys = Unique_GO_IDs[],
                                                    keytype = 'GOALL',
                                                    column = 'SYMBOL',
                                                    multiVals = 'list')
  )

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
