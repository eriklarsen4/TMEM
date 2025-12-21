#' @title get_orthologs_and_aliases
#'
#' @description
#' \strong{get_orthologs_and_aliases} is a function that finds all the aliases
#'    and orthologs of the other four main model organisms from a character
#'    vector of gene symbols of another
#'
#' @param ref_species a string containing the nickname of a model species
#' @param list_of_interest a character vector containing gene symbols
#'
#' @returns
#' \strong{result}: a dataframe housing:
#' \itemize{
#'    \item the reference species
#'    \item the reference gene symbols
#'    \item the reference gene symbols' Ensembl IDs
#'    \item remaining (target) species names
#'    \item target species' orthologous gene symbols
#'    \item target species' orthologous gene Ensembl IDs
#' }
#'
#' @details
#' Queries the Bioconductor db to extract the aliases of the provided gene
#'    symbols, then uses the \strong{orthogene} package to extract orthologous
#'    gene information from the 5 main biological model species, and corrals
#'    into a practical and intuitive dataframe
#'
#' \strong{ref_species} can be one of \strong{"human"}, \strong{"mouse"}, or
#'    \strong{"fly"}
#'
#' @examples
#' \donttest{
#' aDRG_DEG_list <- TMEM::aDRG_DEG_list
#' get_orthologs_and_aliases(ref_species = 'mouse',
#'                          list_of_interest = aDRG_DEG_list[c(1:5)])
#' }
#'
#' @import dplyr
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Dm.eg.db
#' @importFrom AnnotationDbi keys
#' @importFrom AnnotationDbi mapIds
#' @import purrr
#' @import orthogene
#' @import tidyr
#' @import assertthat
#'
#' @references [orthogene](https://github.com/neurogenomics/orthogene)
#'
#' @export
get_orthologs_and_aliases <- function(ref_species, list_of_interest) {
  assertthat::assert_that(!is.null(list_of_interest))
  assertthat::assert_that(!is.null(ref_species))
  assertthat::is.string(ref_species)
  assertthat::is.string(list_of_interest)
  assertthat::not_empty(ref_species)
  assertthat::not_empty(list_of_interest)
  assertthat::assert_that(grepl(ref_species, pattern = '(human)|(mouse)|(fly)'),
                          msg = "ref_species must be exactly one of 'human', 'mouse', or 'fly'!")
  ref_species <- match.arg(ref_species, c("human", "mouse", "fly"))
  # initialize variables ----
  ## list for aliases in the reference species
  listy <- list()

  ## the list for the results
  result <- list()

  "%notin%" <- Negate("%in%")

  ## all the model organisms
  species <- c('human', 'mouse', 'fly', 'macaque', 'zebrafish')
  species <- species[which(species %notin% ref_species == TRUE)]

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
  for (i in seq_len(length(list_of_interest)) ) {
    if ( list_of_interest[i] %in% AnnotationDbi::keys(species_nickname_db_object, keytype = 'SYMBOL') ) {

      listy[[i]] <- suppressMessages(AnnotationDbi::mapIds(species_nickname_db_object, keys = c(list_of_interest[i]), keytype = 'SYMBOL',
                                                           column = 'ALIAS', multiVals = 'list')
      )

    } else {
      next()
    }
  }

  ## convert the aliases list to a dataframe
  ALIAS <- purrr::list_flatten(listy) |>
    purrr::map_df(.f = as.data.frame) |>
    dplyr::rename(Aliases = 1) |>
    dplyr::mutate(SYMBOL = '', .before = 1)

  ## fill the "Symbol" column with the genes from the list of interest
  ## aliases of these genes are the 'Aliases' column
  for (i in seq_len(length(listy)) ) {

    ALIAS$SYMBOL[which(ALIAS$Aliases %in% listy[[i]][[1]] == T)] = listy[[i]] |> names()

  }

  ## convert "N/A" to NA
  ALIAS <- ALIAS |>
    dplyr::group_by(.data$SYMBOL) |>
    dplyr::mutate(Aliases = dplyr::case_when(length(.data$Aliases) == 1 ~ NA, .data$SYMBOL == .data$Aliases ~ '', TRUE ~ .data$Aliases)) |>
    dplyr::ungroup() |>
    dplyr::filter(is.na(.data$Aliases) | .data$Aliases != '')

  # extract ortholog information on all genes (genes in the input list and aliases) ----
  for (i in seq_len(length(species)) ) {

    result[[i]] <- suppressMessages(orthogene::map_orthologs(genes = c(ALIAS |>
                                                                         unlist() |>
                                                                         as.character() |>
                                                                         unique())[which(
                                                                           !is.na(c(ALIAS |>
                                                                                      unlist() |>
                                                                                      as.character() |>
                                                                                      unique())) == TRUE)],
                                                             input_species = ref_species,
                                                             output_species = species[i]) |>
                                      dplyr::mutate(ref_species = as.character(ref_species), .before = 1)
    )
  }

  result <- purrr::list_merge(result) %>%
    purrr::map_df(.f = as.data.frame) |>
    dplyr::mutate(input_number = as.character(.data$input_number)) |>
    tidyr::pivot_longer(cols = c(.data$input_number,
                                 .data$input_ensg,
                                 .data$ortholog_gene,
                                 .data$ortholog_ensg,
                                 .data$description)) %>%
    dplyr::mutate(dplyr::across(c(1:5), ~dplyr::case_when(.x == 'N/A' ~ NA, TRUE ~ .x))) |>
    tidyr::pivot_wider(id_cols = c(ref_species,
                                   .data$input_gene),
                       names_from = .data$name,
                       values_from = .data$value,
                       values_fn = list) |>
    tidyr::unnest(cols = everything()) |>
    dplyr::distinct() |>
    dplyr::mutate(target_species = dplyr::case_when(grepl(.data$ortholog_ensg,
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
