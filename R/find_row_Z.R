#' @title find_row_Z
#' @description
#' \strong{find_row_Z} is a function that determines Z-scores of an expression
#'  matrix
#'
#' @param Expression_Profile a dataframe containing a list of gene symbols
#'  (first column);
#'  remaining columns comprising samples (iterations) of expression values
#'
#' @returns
#' \strong{GeneZ}: a dataframe with Z-scored expression values
#'
#' @details
#' Z-scores are \strong{row-wise} across all columns containing samples
#' (column 2 through the end of the dataframe)
#'
#'
#' @examples
#' \donttest{
#' data(aDRG_TPM)
#' find_row_Z(Expression_Profile = aDRG_TPM)
#' }
#'
#' @import assertthat
#' @import dplyr
#' @import purrr
#'
#' @export
find_row_Z <- function(Expression_Profile){
  assertthat::assert_that(is.data.frame(Expression_Profile), msg = "Expression_Profile must be a dataframe with 2+ numeric columns!")
  assertthat::assert_that(!purrr::is_empty(Expression_Profile), msg = 'Expression_Profile must contain values to compute!')
  assertthat::assert_that(ncol(Expression_Profile) > 2, msg = 'Expression_Profile must be a dataframe with 2+ numeric columns!')
  assertthat::assert_that(!purrr::is_null(Expression_Profile), msg = 'Expression_Profile must not contain any NAs!')
  assertthat::assert_that(grepl(lapply(Expression_Profile, class), pattern = 'numeric|double|float') |> any(), msg = 'Expression_Profile must contain numeric values to compute Z values!')

  message('find_row_Z expects the first column of input dataframe, `Expression_Profile`, to contain character IDs.\nComputations are conducted across c(2:ncol(`Expression_Profile`))!')

  # create an index of column numbers (samples) with expression values
  idx <- which(vapply(Expression_Profile, is.numeric) == TRUE) |> as.numeric()

  Z <- Expression_Profile |>
    dplyr::rowwise() %>%
    dplyr::mutate(row.mean = base::mean(dplyr::c_across(idx))) |> # compute gene mean
    dplyr::mutate(row.sd = stats::sd(dplyr::c_across(idx))) |> # comput gene sd
    dplyr::mutate(
      dplyr::across(idx,~((.x - row.mean) / row.sd))) |> # convert to Z-scores
    dplyr::ungroup() |>
    dplyr::filter(dplyr::if_any(c(idx), ~ !is.na(.))) # remove genes with NaNs

  return(Z)
}
