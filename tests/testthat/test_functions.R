
library(scITD)
library(testthat)
library(Matrix)


# test_that("colMeanVars() functionality", {
#   donor_by_gene <- rbind(c(9,2,1,5),c(3,3,1,2))
#   donor_by_gene <- Matrix(donor_by_gene, sparse = TRUE)
#   result <- colMeanVars(donor_by_gene,rowSel = NULL)
#   expected_result <- cbind(c(6,2.5,1,3.5),c(9,.25,0,2.25),c(2,2,2,2))
#   expected_result <- as.data.frame(expected_result)
#   colnames(expected_result) <- c('m','v','nobs')
#   expect_equal(result, expected_result)
# })

# test_that("form_tensor() functionality", {
#   expected_result <- test_container$scMinimal_ctype[['CD4+ T']]$pseudobulk
#   test_container$scMinimal_ctype[['CD4+ T']]$pseudobulk <- NULL
#   test_container <- form_tensor(test_container, donor_min_cells=0, gene_min_cells=0,
#                                 norm_method='trim', scale_factor=10000,
#                                 vargenes_method='norm_var', vargenes_thresh=500,
#                                 scale_var = TRUE, var_scale_power = 1.5)
#   result <- test_container$scMinimal_ctype[['CD4+ T']]$pseudobulk
#   expect_equal(result, expected_result)
# })


# test_that("tucker() functionality", {
#
#   expected_result <- test_container$tucker_decomp
#
#   # get the tensor
#   tnsr <- test_container$tensor_data[[4]]
#
#   # run just tucker
#   result <- rTensor::tucker(rTensor::as.tensor(tnsr), ranks=c(2,4,2))
#
#   expect_equal(result, expected_result)
# })

# test_that("tucker_ica_helper() functionality", {
#   expected_result <- test_container$pre_tucker_results
#   test_container$pre_tucker_results <- NULL
#
#   tensor_data <- test_container$tensor_data
#   tucker_res <- tucker_ica_helper(tensor_data, ranks=c(2,4,2), tucker_type='regular',
#                                   rotation_type='ica')
#   expect_equal(tucker_res, expected_result)
# })

# test_that("get_factor_exp_var() functionality", {
#   expected_result <- test_container$exp_var
#   result <- c(get_factor_exp_var(test_container,1),get_factor_exp_var(test_container,2))
#   expect_equal(result, expected_result)
# })

test_that("icafast() functionality", {
  expected_result <- test_df[[2]]
  X_dat <- test_df[[1]]
  result <- icafast2(X_dat,2,center=FALSE,alg='def')$S
  expect_equal(result, expected_result)
})

# test_that("icafast() functionality", {
#   expected_result <- test_container$donor_mat_rot
#   dm <- test_container$tucker_decomp$U[[1]]
#   result <- icafast2(dm,2,center=FALSE,alg='def')$S
#   expect_equal(result, expected_result)
# })

# test_that("ica_p1() functionality", {
#   expected_result <- test_container$new_X
#   dm <- test_container$tucker_decomp$U[[1]]
#   result <- ica_p1(dm,2,center=FALSE,alg='def')
#   expect_equal(result, expected_result)
# })

# test_that("ica_p2() functionality", {
#   expected_result <- test_container$my_eigen
#   result <- ica_p2(test_container$dmat,2,center=FALSE,alg='def')
#   expect_equal(result, expected_result)
# })

# test_that("ica_p2() functionality", {
#   expected_result <- test_df[[2]]
#   dm <- test_df[[1]]
#   result <- ica_p2(dm,2,center=FALSE,alg='def')
#   expect_equal(result, expected_result)
# })

# test_that("ica in parts functionality", {
#   expected_result <- test_container$my_eigen
#   dm <- test_container$tucker_decomp$U[[1]]
#   result1 <- ica_p1(dm,2,center=FALSE,alg='def')
#   result <- eigen(result1/5,symmetric=TRUE)$vectors
#   expect_equal(result, expected_result)
# })

# test_that("eigen() functionality", {
#   expected_result <- test_container$my_eigen
#   result <- eigen(test_container$new_X/5,symmetric=TRUE)$vectors
#   expect_equal(result, expected_result)
# })

# test_that("scale functionality", {
#   expected_result <- colMeans(test_container$donor_mat_rot)
#   result <- c(0,0)
#   expect_equal(result, expected_result)
# })

# test_that("kronecker() functionality", {
#   expected_result <- test_container$kron_prod_test
#
#   tensor_data <- test_container$tensor_data
#   gene_nm  <- tensor_data[[2]]
#   ctype_nm  <- tensor_data[[3]]
#   gene_by_factors <- test_container$tucker_decomp$U[[2]]
#   rownames(gene_by_factors) <- gene_nm
#   ctype_by_factors <- test_container$tucker_decomp$U[[3]]
#   rownames(ctype_by_factors) <- ctype_nm
#
#   result <- kronecker(ctype_by_factors,gene_by_factors,make.dimnames = TRUE)
#   expect_equal(result, expected_result)
# })


# test_that("run_tucker_ica() functionality", {
#   expected_result <- test_container$tucker_results
#   test_container <- run_tucker_ica(test_container, ranks=c(2,4,2),
#                                    tucker_type = 'regular', rotation_type = 'ica')
#   result <- test_container$tucker_results
#   expect_equal(result, expected_result)
# })



