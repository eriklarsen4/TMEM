testthat::test_that("get_orthologs_and_aliases requires valid ref_species string", {
  testthat::expect_error(get_orthologs_and_aliases(ref_species = 'monkey', list_of_interest = 'MTOR'))
  testthat::expect_no_error(get_orthologs_and_aliases(ref_species = 'human', list_of_interest = 'MTOR'))
})

testthat::test_that("get_orthologs_and_aliases requires at least one valid gene/protein Entrez ID symbol", {
  testthat::expect_no_error(get_orthologs_and_aliases(ref_species = 'human', list_of_interest = 'MTOR'))
  testthat::expect_error(get_orthologs_and_aliases(ref_species = 'human', list_of_interest = 'ERIKLARSEN'))
  testthat::expect_error(get_orthologs_and_aliases(ref_species = 'human', list_of_interest = ''))
})
