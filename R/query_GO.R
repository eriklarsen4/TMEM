#' @title query_GO
#' @description
#' \strong{`query_GO`} is a function that queries the GO.db with a string
#' term of interest for a given model organism
#'
#' @param model_org a string comprising one of: (for human)
#'  'human', 'HS', 'homo sapiens', (for mouse) 'mouse', 'MM', 'mus musculus',
#'  (for fly) 'drosophila', 'DM', 'fly'
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
#' @importFrom AnnotationDbi, keys
#' @importFrom AnnotationDbi, mapIds
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Dm.eg.db
#' @import stringr
#' @import assertthat
#'
#' @export
query_GO <- function(model_org, string_terms) {
  assertthat::assert_that(!is.null(string_terms))
  assertthat::assert_that(!is.null(model_org))
  assertthat::is.string(model_org)
  assertthat::is.string(string_terms)
  assertthat::not_empty(model_org)
  assertthat::not_empty(string_terms)

  model_org <- match.arg(model_org, c("human", "Hs", "HS", "homo sapiens",
                                      "mouse", "Mm", "MM", "mus musculus",
                                      "fly", "fruit fly", "Dm", "DM", "drosophila melanogaster"))

  # GO_db<- eval(parse(text = GO.db::GO.db))
  abbrev_species_name <- vector()

  if (grepl(model_org, pattern = 'human|HS|homo sapiens', ignore.case = T)) {
    abbrev_species_name <- org.Hs.eg.db
  } else if (grepl(model_org, pattern = 'mouse|MM|mus musculus', ignore.case = T)) {
    abbrev_species_name <- org.Mm.eg.db
  } else if (grepl(model_org, pattern = 'drosophila|DM|fly|fruit fly', ignore.case = T)) {
    abbrev_species_name <- org.Dm.eg.db
  } else {
    message("invalid 'model_org' input!")
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

    GO.list = c(suppressMessages(AnnotationDbi::mapIds(abbrev_species_name,
                                                       keys = c(as.character(
                                                         suppressMessages(AnnotationDbi::mapIds(GO.db::GO.db,
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
                                                         )
                                                       )),
                                                       keytype = 'GOALL',
                                                       column = 'SYMBOL',
                                                       multiVals = 'list')
    )
    )

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
            suppressMessages(AnnotationDbi::mapIds(abbrev_species_name,
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
      )
    aliases = unique(
      as.character(
        unlist(
          suppressMessages(AnnotationDbi::mapIds(abbrev_species_name,
                                                 keys = c(all_unique_genes),
                                                 keytype = 'SYMBOL',
                                                 column = 'ALIAS',
                                                 multiVals = 'list')
          )
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

      string_terms = c(stringr::str_split(as.vector(string_terms), pattern = '\\|', simplify = T))

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
