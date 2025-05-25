testthat::test_that("get_GO_info requires valid species strings", {
  testthat::expect_error(get_GO_info(species = 'Rm'))
})

testthat::test_that("get_GO_info requires at least one valid gene ID", {
  testthat::expect_no_error(get_GO_info(list_of_interest = c("TMEM184B", "ERIKLARSEN"), species = 'HS'))
  testthat::expect_error(get_GO_info(list_of_interest = c("MARTHAB", "ERIKLARSEN"), species = "HS"))
})
