library(scITD)
library(testthat)
library(Matrix)

test_that("colMeanVars() functionality", {
  donor_by_gene <- rbind(c(9,2,1,5),c(3,3,1,2))
  donor_by_gene <- Matrix(donor_by_gene, sparse = TRUE)
  result <- colMeanVars(donor_by_gene,rowSel = NULL)
  expected_result <- cbind(c(6,2.5,1,3.5),c(9,.25,0,2.25),c(2,2,2,2))
  expected_result <- as.data.frame(expected_result)
  colnames(expected_result) <- c('m','v','nobs')
  expect_equal(result, expected_result)
})

test_that("get_ctype_data() functionality", {
  expected_result <- test_container$scMinimal_ctype[['CD4+ T']]$data_sparse
  test_container$scMinimal_ctype <- NULL
  test_container <- get_ctype_data(test_container,donor_min_cells = 0)
  result <- test_container$scMinimal_ctype[['CD4+ T']]$data_sparse
  expect_equal(result, expected_result)
})

test_that("run_tucker_ica() functionality", {
  expected_result <- test_container$tucker_results
  test_container$tucker_results <- NULL
  test_container <- run_tucker_ica(test_container, ranks=c(2,4,2),
                                   shuffle=FALSE)
  result <- test_container$tucker_results
  expect_equal(result, expected_result)
})


