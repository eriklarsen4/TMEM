testthat::test_that("query_GO requires valid species strings", {
  testthat::expect_error(query_GO(model_org = 'Rm', string_terms = 'mTORC1'))
})

testthat::test_that("query_GO requires at least one valid GO ID term", {
  testthat::expect_no_error(query_GO(model_org = 'HS', string_terms = 'vacuolar acidification'))
  testthat::expect_error(query_GO(model_org = 'HS', string_terms = 'Itch'))
})
