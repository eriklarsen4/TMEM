#'
#' aDRG_DEG_list dataset
#'
#' List of differentially expressed genes as determined by DESeq2 Wald test on bulk RNA extracted from mouse dorsal root ganglion neurons of Tmem184b-mutant and Wild-Type (C57BL/6) mice
#'
#' @docType data
#' @format one of 3 datasets for the package to test and execute package functions:
#' \describe{
#'    \item{aDRG_DEG_list}{a character vector of differentially expressed Entrez gene identifiers for mus musculus build 9, according to the DESeq2 Wald test output of 2 groups (Wild-Type and Mutant) of 4 samples of bulk RNA derived from mouse dorsal root ganglion neurons}
#'    }
#'
#' @details
#' Additional manipulations on the data were made to exclude pseudogenes, unannotated genes, ribosomal genes, and genes not containing an adjusted P-value.
#' Genes were then filtered by those containing an adjusted P-value less than 0.05.
#' The code used to derive these genes was as follows:
#' \code{
#' library(readxl)
#' aDRG_DEG_list <- read_xlsx('path_to_file_from_source/NIHMS1770888-supplement-SDC_1_-_Adult_DRG_RNAseq.xlsx') |>
#'  dplyr::rename(Base.Mean = `Base Mean`,
#'                log2.FC = `log2(FC)`,
#'                Wald.stats = `Wald stats`,
#'                `\%WT` = `\% WT`) |>
#'  dplyr::filter(!is.na(AdjP)) |>
#'  dplyr::filter(!grepl(GeneID,
#'                        pattern = 'Rps.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|Gm.+.?$|^[A-Z+]+[A-Z].+.?$|^[0-9]+.+.?$|mt.+.?$')) |>
#'  dplyr::filter(is.na(`...9`)) |> # removes Nppb which we noted as having been excluded by DESeq2 because of persistently low expression across mutant samples
#'  dplyr::select(c(1:8)) |>
#'  dplyr::mutate(AdjP = as.numeric(AdjP)) |>
#'  dplyr::filter(AdjP <= 0.05) |>
#'  dplyr::filter(GeneID != "Tmem184b") |>
#'  dplyr::select(GeneID) |>
#'  unlist() |>
#'  as.character()
#' }
#'
#' @source [TMEM184B is necessary for IL-31-induced itch](https://pmc.ncbi.nlm.nih.gov/articles/PMC8854445/#SD5)
#' @references [DESeq2 Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
#'
"aDRG_DEG_list"
