
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

test_that("form_tensor() functionality", {
  expected_result <- test_container$scMinimal_ctype[['CD4+ T']]$pseudobulk
  test_container$pseudobulk <- NULL
  test_container <- form_tensor(test_container, donor_min_cells=0, gene_min_cells=0,
                                norm_method='trim', scale_factor=10000,
                                vargenes_method='norm_var', vargenes_thresh=500,
                                scale_var = TRUE, var_scale_power = 1.5)
  result <- test_container$scMinimal_ctype[['CD4+ T']]$pseudobulk
  expect_equal(result, expected_result)
})


test_that("tucker() functionality", {

  expected_result <- test_container$tucker_decomp
  test_container$tucker_decomp <- NULL

  # get the tensor
  tnsr <- test_container$tensor_data[[4]]

  # run just tucker
  tucker_decomp <- rTensor::tucker(rTensor::as.tensor(tnsr), ranks=c(2,4,2))

  expect_equal(tucker_decomp, expected_result)
})


test_that("run_tucker_ica() functionality", {
  expected_result <- test_container$tucker_results
  test_container$tucker_results <- NULL
  test_container <- run_tucker_ica(test_container, ranks=c(2,4,2),
                                   tucker_type = 'regular', rotation_type = 'ica')
  result <- test_container$tucker_results
  expect_equal(result, expected_result)
})






