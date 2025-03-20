#'
#' aDRG_TPM dataset
#'
#' DESeq2 normalized counts dataframe for bulk RNA extracted from mouse dorsal root ganglion neurons of Tmem184b-mutant and Wild-Type (C57BL/6) mice
#'
#' @docType data
#' @format one of 3 datasets for the package to test and execute package functions:
#' \describe{
#'    \item{aDRG_TPM}{a dataframe of the DESeq2 estimations of normalized counts of 2 groups (Wild-Type and Mutant) of 4 samples of bulk RNA derived from mouse dorsal root ganglion neurons}
#'    \item{WT1}{the sample derived from Wild-Type mouse 1; numeric}
#'    \item{WT2}{the sample derived from Wild-Type mouse 2; numeric}
#'    \item{WT3}{the sample derived from Wild-Type mouse 3; numeric}
#'    \item{WT4}{the sample derived from Wild-Type mouse 4; numeric}
#'    \item{Mut1}{the sample derived from Tmem184b-mutant mouse 1; numeric}
#'    \item{Mut2}{the sample derived from Tmem184b-mutant mouse 2; numeric}
#'    \item{Mut3}{the sample derived from Tmem184b-mutant mouse 3; numeric}
#'    \item{Mut4}{the sample derived from Tmem184b-mutant mouse 4; numeric}
#'    }
#'
#' @details
#' Data was not included in paper submission, and was used to generate heatmaps, thus, this dataset is raw/in-house.
#' Additional manipulations on the data were made for formatting purposes.
#' These manipulations were as follows:
#' \code{
#' aDRG_TPM <- read.csv('path_to_file/RNASeqRepResults.csv')
#` colnames(aDRG_TPM) <-  c("GeneID",
#`                         "WT1", "WT2", "WT3", "WT4",
#`                         "Mut1", "Mut2", "Mut3", "Mut4")
#' }
#'
#'
#'
#' @source [TMEM184B is necessary for IL-31-induced itch](https://pmc.ncbi.nlm.nih.gov/articles/PMC8854445/#SD5)
#' @references [DESeq2 Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
#'
"aDRG_TPM"
