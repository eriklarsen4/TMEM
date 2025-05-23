
devtools::document(pkg = 'C:/Users/Erik/Desktop/Programming/R/Bio/TMEM')
devtools::build(pkg = 'C:/Users/Erik/Desktop/Programming/R/Bio/TMEM', path = 'C:/Users/Erik/Desktop/Programming/R/Bio/TMEM', binary = T, vignettes = F)

devtools::check(pkg = 'C:/Users/Erik/Desktop/Programming/R/Bio/TMEM', document = T, vignettes = F, cran = T)
devtools::load_all(export_all = T)

devtools::install('C:/Users/Erik/Desktop/Programming/R/Bio/TMEM', reload = T, build = T, dependencies = T)
devtools::install('C:/Users/Erik/Desktop/Programming/R/Bio/TMEM', reload = T, build = F, dependencies = T)

# pkgbuild::check_build_tools(debug = TRUE)

usethis::use_news_md()

library(testthat)
usethis::use_github_action("check-standard")
install.packages("covr")
library(covr)

usethis::use_news_md()
usethis::use_vignette("TMEM")
testthat::test_that("TMEM")

roxygen2::roxygenize(clean = TRUE, package.dir = 'C:/Users/Erik/Desktop/Programming/R/Bio/TMEM')
roxygen2::
# data load ----
library(tidyverse)
aDRG_DEA_results <- read.csv('C:/Users/Erik/Downloads/DESeq2 Expression Results.csv') |>
  dplyr::rename(`%WT` = `X..WT`) |>
  dplyr::filter(!is.na(AdjP)) |> # remove any NAs in p-value columns
  dplyr::filter(!grepl(GeneID, # remove rRNAs, mitochondrial tRNAs, pseudogenes
                       pattern = 'Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$|mt.+.?$'))
usethis::use_data(aDRG_DEA_results)

aDRG_TPM <- read.csv('C:/Users/Erik/Downloads/RNASeqRepResults.csv')
colnames(aDRG_TPM) <-  c("GeneID",
                          "WT1", "WT2", "WT3", "WT4",
                          "Mut1", "Mut2", "Mut3", "Mut4")
usethis::use_data(aDRG_TPM, overwrite = TRUE)
library(readxl)
aDRG_DEA_results <- read_xlsx('C:/Users/Erik/Documents/NIHMS1770888-supplement-SDC_1_-_Adult_DRG_RNAseq.xlsx') |> dplyr::rename(Base.Mean = `Base Mean`, log2.FC = `log2(FC)`, Wald.stats = `Wald stats`, `%WT` = `% WT`) |>
  dplyr::filter(!is.na(AdjP)) |>
  dplyr::filter(!grepl(GeneID,
                       pattern = 'Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$|mt.+.?$')) |> dplyr::filter(is.na(`...9`)) |>
  dplyr::select(c(1:8))


aDRG_DEG_list <- read_xlsx('C:/Users/Erik/Documents/NIHMS1770888-supplement-SDC_1_-_Adult_DRG_RNAseq.xlsx') |> dplyr::rename(Base.Mean = `Base Mean`, log2.FC = `log2(FC)`, Wald.stats = `Wald stats`, `%WT` = `% WT`) |>
  dplyr::filter(!is.na(AdjP)) |>
  dplyr::filter(!grepl(GeneID,
                       pattern = 'Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$|mt.+.?$')) |>
  dplyr::select(c(1:8)) |>
  dplyr::mutate(AdjP = as.double(AdjP)) |>
  dplyr::filter(AdjP <= 0.05) |>
  dplyr::filter(GeneID != "Tmem184b") |>
  dplyr::select(GeneID) |>
  unlist() |>
  as.character()
usethis::use_data(aDRG_DEG_list)
# testing ----
usethis::use_testthat()
usethis::use_test()

results <- usethis::use_cran_comments()

usethis::use_travis()

usethis::use_github_action()

#install rhub
install.packages("rhub")

# check and release ----
#sign in
rhub::validate_email()
# check
devtools::check_rhub()

devtools::check_win_devel()

cran_checks <- rhub::check_for_cran()

#check package viability on all OSs
devtools::check_win_devel()

results$cran_summary()

devtools::test_coverage()

#release to cran
devtools::release()
devtools::spell_check()

##
gitcreds::gitcreds_set(url = 'https://www.github.com/eriklarsen4/TMEM')
rhub::rhub_check(gh_url = 'https://www.github.com/eriklarsen4/TMEM')
usethis::create_github_token()

usethis::use_spell_check(vignettes = TRUE, lang = 'en-US', error = FALSE)

# misc ----
install.packages("annotationDbi", lib = 'C:/Users/Erik/AppData/Local/R/win-library/4.3')
BiocManager::install("AnnotationDbi")
BiocManager::install("AnnotationDbi", lib = 'C:/Users/Erik/AppData/Local/R/win-library/4.3')

# ----
