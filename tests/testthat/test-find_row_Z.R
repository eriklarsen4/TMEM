testthat::test_that("find_row_Z requires at least two numeric columns in input dataframe", {
  testthat::expect_error(find_row_Z(Expression_Profile = TMEM::aDRG_TPM[c(1:5),c(1,2)]))
  testthat::expect_gt(length(which(grepl(lapply(TMEM::aDRG_TPM[c(1:5),], class), pattern = 'numeric|double|float') == T)), 2)
  testthat::expect_no_error(find_row_Z(Expression_Profile = TMEM::aDRG_TPM[c(1:5),c(2:9)]))
})

testthat::test_that("find_row_Z expects Expression_Profile's first column to be IDs", {
  testthat::expect_no_error(find_row_Z(Expression_Profile = TMEM::aDRG_TPM[c(1:5),c(2:9)]))
  testthat::expect_error(find_row_Z(Expression_Profile = c(TMEM::aDRG_TPM |> dplyr::mutate(col1 = 1, .before = 1) |> dplyr::slice(c(1:5)))))
})
