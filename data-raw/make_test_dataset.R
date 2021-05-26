
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


## now saving the raw tucker results for testing purposes
# get the tensor
tnsr <- test_container$tensor_data[[4]]

# run tucker
tucker_decomp <- rTensor::tucker(rTensor::as.tensor(tnsr), ranks=c(2,4,2))
test_container$tucker_decomp <- tucker_decomp


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

# set.seed(123)
# nobs <- 1000
# Amat <- cbind(icasamp("a","rnd",nobs),icasamp("b","rnd",nobs))
# Bmat <- matrix(2*runif(4),2,2)
# X_dat <- tcrossprod(Amat,Bmat)
# test_res <- ica::icafast(X_dat,2,center=FALSE,alg='def')$S
# test_df <- list(X_dat,test_res)
# save(test_df,file='/home/jmitchel/scITD/data/test_df.RData',compress = "xz")

X_dat <- test_container$tucker_decomp$U[[1]]
test_res <- ica::icafast(X_dat,2,center=FALSE,alg='def')$S
test_df <- list(X_dat,test_res)
save(test_df,file='/home/jmitchel/scITD/data/test_df.RData',compress = "xz")








