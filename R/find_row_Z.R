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
#' @import assertthat
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
    MeansAndSDs$mean[i] <- (sum(Expression_Profile[i,c(2:ncol(Expression_Profile))])/ncol(Expression_Profile[,c(2:ncol(Expression_Profile))])) # row mean
  }
  # for (j in 1:nrow(MeansAndSDs)){
  #   MeansAndSDs$sd[j] <- sd(Expression_Profile[j,c(2:ncol(Expression_Profile))], na.rm = TRUE) # row sd
  # }

  for (i in 1:nrow(MeansAndSDs)){
    for (j in 2:ncol(MeansAndSDs)){
      MeansAndSDs$sd[i] <- sum(
        (Expression_Profile[i,j] - (sum(Expression_Profile[i,c(2:ncol(Expression_Profile))])/ncol(Expression_Profile[,c(2:ncol(Expression_Profile))]))
        )^2)/ncol(Expression_Profile[,c(2:ncol(Expression_Profile))])
    }
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
