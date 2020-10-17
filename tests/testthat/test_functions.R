library(scITD)
library(testthat)

test_that("colMeanVars() functionality", {
  donor_by_gene <- rbind(c(9,2,1,5),c(3,3,1,2))
  donor_by_gene <- methods::as(donor_by_gene,'sparseMatrix')
  result <- colMeanVars(donor_by_gene,rowSel = NULL)
  expected_result <- cbind(c(6,2.5,1,3.5),c(9,.25,0,2.25),c(2,2,2,2))
  expected_result <- as.data.frame(expected_result)
  colnames(expected_result) <- c('m','v','nobs')
  expect_equal(result, expected_result)
})
