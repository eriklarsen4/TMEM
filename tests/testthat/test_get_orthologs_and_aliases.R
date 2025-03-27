library(testthat)
library(assertthat)
library(TMEM)
# define function
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
  assertthat::assert_that(ref_species %notin% species)

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

  assertthat::assert_that(colnames(result) %in% c('ref_species', 'target_species',
                                                'input_gene', 'ortholog_gene',
                                                'input_ensg', 'ortholog_ensg',
                                                'description') |> all())
  # output ----
  # result <<- result
  return(result)
}
data("aDRG_DEG_list")

# test whether species input parameter is valid ----
testthat::test_that("get_orthologs_and_aliases ref_species is correct string option", {
  testthat::expect_error(
    TMEM::get_orthologs_and_aliases(ref_species = 'cow',
                                    list_of_interest = aDRG_DEG_list))
})

# test whether list of interest contains valid parameters ----
testthat::test_that("get_orthologs_and_aliases input parameters are not omitted", {
  testthat::expect_error(
    TMEM::get_orthologs_and_aliases(ref_species = NULL,
                                    list_of_interest = aDRG_DEG_list)
  )
  testthat::expect_error(
    TMEM::get_orthologs_and_aliases(ref_species = 'human',
                                    list_of_interest = NULL)
  )
  testthat::expect_error(
    TMEM::get_orthologs_and_aliases(ref_species = NULL,
                                    list_of_interest = NULL)
  )
})

# test that the output is a non-empty dataframe ----
testthat::test_that("get_orthologs_and_aliases returns a tibble with 7 columns",{
  result <- TMEM::get_orthologs_and_aliases(ref_species = 'mouse', list_of_interest = 'Mtor')
  testthat::expect_message(regexp = "'select()' returned 1:many mapping between keys and columns")
  testthat::expect_message(regexp = "'select()' returned 1:1 mapping between keys and columns")
  testthat::expect_is(result, class = c(base::class(result)))
  testthat::expect_length(nrow(result), n = 7)
  testthat::expect_gt(nrow(result), expected = 1)

})
