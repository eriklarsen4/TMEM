library(testthat)
library(assertthat)
library(TMEM)
# define function
query_GO <- function(model_org, string_terms) {
  assertthat::assert_that(grepl(model_org,
                                pattern = 'human|HS|homo sapiens|mouse|mm|mus musculus|fly|fruit fly|dm|drosophila|drosophila melanogaster'))
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

# test that model org is valid string option ----
testthat::test_that("query_GO has valid model_org user input", {
  testthat::expect_error(TMEM::query_GO(model_org = 'mouse', string_terms = 'dense core vesicle'))
})

# test that string terms are valid ----
testthat::test_that("query_GO has valid string_terms input", {
  testthat::expect_error(TMEM::query_GO(model_org = 'mouse', string_terms = ''))
})

# test that no input parameters are omitted ----
testthat::test_that("query_GO has no omitted parameters", {
  testthat::expect_error(TMEM::query_GO(model_org = NULL, string_terms = 'vesicle'))
  testthat::expect_error(TMEM::query_GO(model_org = 'human', string_terms = NULL))
  testthat::expect_error(TMEM::query_GO(model_org = NULL, string_terms = NULL))
})

# test that output contains appropriate elements ----
testthat::test_that("query_GO returns a list of character vectors, lists, and dataframes", {
  result <- TMEM::query_GO(model_org = 'mouse', string_terms = 'vesicle')
  testthat::expect_message(regexp = "'select()' returned 1:many mapping between keys and columns")
  testthat::expect_message(regexp = "'select()' returned 1:1 mapping between keys and columns")
  testthat::expect_gte(nrow(result$GO_df), expected = 1)
  testthat::expect_gte(length(result$GO.list |> unlist() |> as.character()), expected = 1)
  testthat::expect_gte(length(result$GO.IDs), expected = 1)
  testthat::expect_gte(length(result$GO.Terms), expected = 1)
  testthat::expect_gte(length(result$all_unique_genes), expected = 1)
  testthat::expect_gte(length(result$aliases), expected = 1)
})
