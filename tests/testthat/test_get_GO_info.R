library(testthat)
library(TMEM)
# define function
get_GO_info <- function(list_of_interest, species) {
  # First, define the species of interest ----
  if ( grepl(species, pattern = 'human|HS|homo sapiens', ignore.case = T) ) {
    abbrev_species_name <- org.Hs.eg.db
  } else if ( grepl(species, pattern = 'mouse|Mm|MM|Mus musculus', ignore.case = T) ) {
    abbrev_species_name <- org.Mm.eg.db
  } else if ( grepl(species,
                    pattern = 'drosophila melanogaster|drosophila|fly|DM|fruit fly', ignore.case = T) ) {
    abbrev_species_name <- org.Dm.eg.db
  } else {
    message("Invalid species name! Please see documentation for options!")
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
  # Find the genes/proteins in the input list that don't have Ensembl/Entrez IDs ----
  indeces_without_Ensembl_IDs <- c(
    which(is.na(Ensembl_IDs[]) == TRUE )
  )

  # Find the genes/proteins that don't have GO IDs ----
  indeces_without_GO_IDs <- as.integer(c(
    which(is.na(GO_IDs) == TRUE)
  ))

  # Remove genes/proteins that don't have IDs ----
  if ( length(indeces_without_Ensembl_IDs) >= 1 &
       length(indeces_without_GO_IDs) == 0 ) {

    list_of_interest <- list_of_interest[-indeces_without_Ensembl_IDs]

  } else if (

    length(indeces_without_Ensembl_IDs) == 0 &
    length(indeces_without_GO_IDs) >= 1

  ) {

    list_of_interest <- list_of_interest[-indeces_without_GO_IDs]

  } else if (

    length(indeces_without_Ensembl_IDs) >= 1 &
    length(indeces_without_GO_IDs) >= 1

  ) {

    list_of_interest <- list_of_interest[-c(unique(indeces_without_GO_IDs, indeces_without_Ensembl_IDs))]

  } else {

    list_of_interest <- list_of_interest

  }

  # Re-query the GO IDs having filtered the list ----
  GO_IDs <- c(
    unique(
      AnnotationDbi::mapIds(abbrev_species_name,
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
  colnames(GENE_GO_INFO) <- c("GeneID", "GO_Term_Size", "GO_Term_IDs", "GO_Terms")

  ## Fill the GO Term ID column with the concatenated GO Term IDs associated with each gene/protein
  ## skip over the genes/proteins without any
  for (i in 1:length(GO_IDs)) {

    GENE_GO_INFO[i,1] <- list_of_interest[i]

    ## Initialize numeric columns with 0s
    GENE_GO_INFO[i,2] <- 0
    ## Initialize character columns with "None"s
    GENE_GO_INFO[i, c(3,4)] <- "None"

    GENE_GO_INFO[i,3] <- paste0(names(AnnotationDbi::mapIds(GO.db, GO_IDs[[i]][], 'TERM', 'GOID')), collapse = ';')

  }

  ## Extract the GO IDs and remove the redundant ones;
  ## Fill the 3rd column with the numbers of GO Terms associated with each gene/protein
  ## skip over the genes/proteins without any
  for (i in 1:length(GO_IDs)) {

    x <- stringr::str_split(as.vector(GENE_GO_INFO[i,3]), pattern = ';', simplify = TRUE)
    x <- c(as.character(x[which(x %in% x == TRUE )]))
    GENE_GO_INFO[i,3] <- paste0(unique(x), collapse = ';')
    GENE_GO_INFO[i,2] <- paste0(as.numeric(length(unique(x))))
  }

  ## Fill the 4th column with the concatenated GO Terms associated with each gene/protein
  ## skip over the genes/proteins without any
  for (i in 1:length(GO_IDs)) {

    x <- str_split(as.vector(GENE_GO_INFO[i,3]), pattern = ';', simplify = TRUE)
    GENE_GO_INFO[i,4] <- paste0(as.character(AnnotationDbi::mapIds(GO.db, keys = x, keytype = 'GOID', 'TERM')), collapse = ';')
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
  GENE_GO_INFO_df <- as.data.frame(GENE_GO_INFO) |> dplyr::mutate(Num_GO_Term_Size = as.numeric(.data$GO_Term_Size))


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
"%notin%" <- Negate("%in%")
data("aDRG_DEG_list")

# test whether the input string for species is valid ----
testthat::test_that("get_GO_info species input is correct string option", {
  ok_species <- c('human', 'HUMAN', 'Human', 'hs', 'HS', 'homo sapiens',
               'Homo sapiens', 'Homo Sapiens', 'HOMO SAPIENS',
               'mouse', 'MOUSE', 'Mouse', 'mm', 'MM', 'mus musculus',
               'Mus musculus', 'Mus Musculus', 'MUS MUSCULUS',
               'fly', 'FLY', 'Fly', 'fruit fly', 'FRUIT FLY', 'Fruit fly',
               'Fruit Fly', 'dm', 'DM', 'Dm', 'drosophila melanogaster',
               'DROSOPHILA MELANOGASTER', 'Drosophila Melanogaster',
               'Drosophila melanogaster')

  testthat::expect_error(TMEM::get_GO_info(list_of_interest = aDRG_DEG_list, species = 'macaque'))

})

# test whether the input list has at least one valid string ----
testthat::test_that("get_GO_info list_of_interest input has 1+ valid gene", {
  testthat::expect_error(TMEM::get_GO_info(list_of_interest = 'Z', species = 'hs'))
})

# test whether required input parameters were specified ----
testthat::test_that("get_GO_info both parameters specified", {
  testthat::expect_error(TMEM::get_GO_info(list_of_interest = NULL, species = 'hs'))
  testthat::expect_error(TMEM::get_GO_info(list_of_interest = aDRG_DEG_list, species = NULL))
})


# test that the returned object is a list ----
testthat::test_that("get_GO_info output is list of 6 objects", {
  result <- TMEM::get_GO_info(list_of_interest = aDRG_DEG_list, species = 'mm')
  testthat::expect_is(result, 'list')
  testthat::expect_length(result, n = 6)
})

# test the returned object contents are not empty ----
testthat::test_that("get_GO_info output does not contain empty elements", {
  result <- TMEM::get_GO_info(list_of_interest = aDRG_DEG_list, species = 'mm')
  testthat::expect_gt(result$list_of_interest_aliases, expected = 1)
  testthat::expect_gt(result$unique_GOs, expected = 1)
  testthat::expect_gt(result$unique_GO_IDs, expected = 1)
  testthat::expect_gt(nrow(result$GO_info_by_term_df), expected = 1)
  testthat::expect_gt(result$gene_GO_info_df, expected = 1)

  testthat::expect_error(TMEM::get_GO_info(list_of_interest = 'Z', species = 'mm'))
})
