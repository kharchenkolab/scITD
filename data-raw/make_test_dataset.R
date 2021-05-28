
pbmc_counts <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_counts_v2.rds')

pbmc_meta <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_meta_v2.rds')

feature.names <- readRDS('/home/jmitchel/data/van_der_wijst/genes.rds')

# need to do this one donor at a time!!!
small_meta1 <- pbmc_meta[pbmc_meta$donors == 's5',]
small_meta2 <- pbmc_meta[pbmc_meta$donors == 's40',]
small_meta3 <- pbmc_meta[pbmc_meta$donors == 's12',]
small_meta4 <- pbmc_meta[pbmc_meta$donors == 's20',]
small_meta5 <- pbmc_meta[pbmc_meta$donors == 's22',]

ct1 <- 'CD4+ T'
ct2 <- 'CD8+ T'

cells_use_1_ct1 <- sample(rownames(small_meta1[small_meta1$ctypes==ct1,]),5)
cells_use_1_ct2 <- sample(rownames(small_meta1[small_meta1$ctypes==ct2,]),5)
cells_use_2_ct1 <- sample(rownames(small_meta2[small_meta2$ctypes==ct1,]),5)
cells_use_2_ct2 <- sample(rownames(small_meta2[small_meta2$ctypes==ct2,]),5)
cells_use_3_ct1 <- sample(rownames(small_meta3[small_meta3$ctypes==ct1,]),5)
cells_use_3_ct2 <- sample(rownames(small_meta3[small_meta3$ctypes==ct2,]),5)
cells_use_4_ct1 <- sample(rownames(small_meta4[small_meta4$ctypes==ct1,]),5)
cells_use_4_ct2 <- sample(rownames(small_meta4[small_meta4$ctypes==ct2,]),5)
cells_use_5_ct1 <- sample(rownames(small_meta5[small_meta5$ctypes==ct1,]),5)
cells_use_5_ct2 <- sample(rownames(small_meta5[small_meta5$ctypes==ct2,]),5)

cells_use <- c(cells_use_1_ct1,cells_use_1_ct2,cells_use_2_ct1,cells_use_2_ct2,
               cells_use_3_ct1,cells_use_3_ct2,cells_use_4_ct1,cells_use_4_ct2,
               cells_use_5_ct1,cells_use_5_ct2)

baby_counts <- pbmc_counts[1:30,cells_use]
baby_meta <- pbmc_meta[cells_use,]

# counts are too sparse for tests so just using random matrix of ints...
baby_counts_rand <- matrix(sample.int(15, size = ncol(baby_counts)*nrow(baby_counts),
                                 replace = TRUE), nrow = nrow(baby_counts),
                      ncol = ncol(baby_counts))
colnames(baby_counts_rand) <- colnames(baby_counts)
rownames(baby_counts_rand) <- rownames(baby_counts)


baby_counts_rand <- Matrix(baby_counts_rand, sparse = TRUE)


param_list <- initialize_params(ctypes_use = c("CD4+ T", "CD8+ T"),
                                ncores = 2, rand_seed = 10)

test_container <- make_new_container(count_data=baby_counts_rand, meta_data=baby_meta,
                                     params=param_list)

test_container <- form_tensor(test_container, donor_min_cells=0, gene_min_cells=0,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var', vargenes_thresh=500,
                              scale_var = TRUE, var_scale_power = 1.5)

test_container <- run_tucker_ica(test_container, ranks=c(2,4,2),
                                 tucker_type = 'regular', rotation_type = 'ica')

test_container <- run_jackstraw(test_container, ranks=c(2,4,2),
                                n_fibers=10, n_iter=500, tucker_type='regular',
                                rotation_type='ica')

save(test_container,file='/home/jmitchel/scITD/data/test_container.RData',compress = "xz")


## now saving the raw tucker results for testing purposes
# get the tensor
tnsr <- test_container$tensor_data[[4]]

# run tucker
tucker_decomp <- rTensor::tucker(rTensor::as.tensor(tnsr), ranks=c(2,4,2))
test_container$tucker_decomp <- tucker_decomp
test_container$dmat <- tucker_decomp$U[[1]]
my_eigen <- ica_p2(test_container$dmat,2,center=FALSE,alg='def')
test_container$my_eigen <- my_eigen
save(test_container,file='/home/jmitchel/scITD/data/test_container.RData',compress = "xz")


## now saving the raw tucker helper results
# run tucker
tensor_data <- test_container$tensor_data
tucker_res <- tucker_ica_helper(tensor_data, ranks=c(2,4,2), tucker_type='regular',
                                rotation_type='ica')

test_container$pre_tucker_results <- tucker_res


# adding rotation results as a test
dm <- test_container$tucker_decomp$U[[1]]
donor_mat_rot <- ica::icafast(dm,2,center=FALSE,alg='def')$S
test_container$donor_mat_rot <- donor_mat_rot

save(test_container,file='/home/jmitchel/scITD/data/test_container.RData',compress = "xz")


# calculating first and second parts of ica fn for testing
dm <- test_container$tucker_decomp$U[[1]]
new_X <- ica_p1(dm,2,center=FALSE,alg='def')
my_eigen <- ica_p2(dm,2,center=FALSE,alg='def')
test_container$new_X <- new_X
test_container$my_eigen <- my_eigen
save(test_container,file='/home/jmitchel/scITD/data/test_container.RData',compress = "xz")


# # store result for testing kronecker product
# tensor_data <- test_container$tensor_data
# gene_nm  <- tensor_data[[2]]
# ctype_nm  <- tensor_data[[3]]
# gene_by_factors <- test_container$tucker_decomp$U[[2]]
# rownames(gene_by_factors) <- gene_nm
# ctype_by_factors <- test_container$tucker_decomp$U[[3]]
# rownames(ctype_by_factors) <- ctype_nm
#
# kron_prod <- kronecker(ctype_by_factors,gene_by_factors,make.dimnames = TRUE)
# test_container$kron_prod_test <- kron_prod

# save(test_container,file='/home/jmitchel/scITD/data/test_container.RData',compress = "xz")



# making separate df to test ica

set.seed(123)
nobs <- 50
Amat <- cbind(ica::icasamp("a","rnd",nobs),ica::icasamp("b","rnd",nobs))
Bmat <- matrix(2*runif(4),2,2)
X_dat <- tcrossprod(Amat,Bmat)
# X_dat <- matrix( rnorm(10,mean=0,sd=1), 5, 2)
test_res <- ica::icafast(X_dat,2,center=FALSE,alg='def')$S
test_df <- list(X_dat,test_res)
save(test_df,file='/home/jmitchel/scITD/data/test_df.RData',compress = "xz")

X_dat <- test_container$tucker_decomp$U[[1]]
test_res <- ica::icafast(X_dat,2,center=FALSE,alg='def')$S
test_df <- list(X_dat,test_res)
save(test_df,file='/home/jmitchel/scITD/data/test_df.RData',compress = "xz")



## making new test_container object
pbmc_counts <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_counts_v2.rds')
pbmc_meta <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_meta_v2.rds')

param_list <- initialize_params(ctypes_use = c("CD4+ T", "CD8+ T"),
                                ncores = 2, rand_seed = 10)

test_container2 <- make_new_container(count_data=pbmc_counts, meta_data=pbmc_meta,
                                     params=param_list)

test_container2 <- form_tensor(test_container2, donor_min_cells=0, gene_min_cells=0,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var', vargenes_thresh=500,
                              scale_var = TRUE, var_scale_power = 1.5)

tnsr <- test_container2$tensor_data[[4]]
tucker_decomp <- rTensor::tucker(rTensor::as.tensor(tnsr), ranks=c(2,4,2))
X_dat <- tucker_decomp$U[[1]]
set.seed(123)
test_res <- ica::icafast(X_dat,2,center=FALSE,alg='def',maxit = 200,tol = 1e-15)$S
# test_res <- fastICA(X_dat, 2, alg.typ = "deflation",
#         fun = "logcosh", alpha = 1.0, method = 'R',
#         row.norm = FALSE, maxit = 1000, tol = 1e-10, verbose = FALSE,
#         w.init = NULL)$S
test_df <- list(X_dat,test_res)
save(test_df,file='/home/jmitchel/scITD/data/test_df.RData',compress = "xz")


## seeing if icafast and fastica give the same results
library(fastICA)
r1 <- fastICA(X_dat, 2, alg.typ = "deflation",
        fun = "logcosh", alpha = 1.0, method = 'R',
        row.norm = FALSE, maxit = 100, tol = 1e-06, verbose = FALSE,
        w.init = NULL)
r2 <- ica::icafast(X_dat,2,center=FALSE,alg='def')$S
head(r1$S)
head(r2)









# ## making new test_container object
# pbmc_counts <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_counts_v2.rds')
# pbmc_meta <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_meta_v2.rds')
#
#
# param_list <- initialize_params(ctypes_use = c("CD4+ T", "CD8+ T"),
#                                 ncores = 2, rand_seed = 10)
#
# test_container2 <- make_new_container(count_data=pbmc_counts, meta_data=pbmc_meta,
#                                      params=param_list)
#
# test_container2 <- form_tensor(test_container2, donor_min_cells=0, gene_min_cells=0,
#                               norm_method='trim', scale_factor=10000,
#                               vargenes_method='norm_var', vargenes_thresh=50,
#                               scale_var = TRUE, var_scale_power = 1.5)
#
# test_container2 <- run_tucker_ica(test_container2, ranks=c(2,4,2),
#                                  tucker_type = 'regular', rotation_type = 'ica')
#
# tensor_data <- test_container2$tensor_data
# tucker_results <- test_container2$tucker_results
# test_container2 <- list(tensor_data=tensor_data,tucker_results=tucker_results)
#
# save(test_container2,file='/home/jmitchel/scITD/data/test_container2',compress = "xz")






