
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

test_container <- form_tensor(test_container, donor_min_cells=0,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var', vargenes_thresh=500,
                              scale_var = TRUE, var_scale_power = .5)

test_container <- run_tucker_ica(test_container, ranks=c(2,4),
                                 tucker_type = 'regular', rotation_type = 'hybrid')

test_container <- get_lm_pvals(test_container)

save(test_container,file='/home/jmitchel/scITD/data/test_container.RData',compress = "xz")










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
# test_df <- list(tensor_data=tensor_data,tucker_results=tucker_results)
#
# save(test_df,file='/home/jmitchel/scITD/data/test_df.RData',compress = "xz")






