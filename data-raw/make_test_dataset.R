
pbmc_counts <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_counts_v2.rds')

pbmc_meta <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_meta_v2.rds')

feature.names <- readRDS('/home/jmitchel/data/van_der_wijst/genes.rds')

# need to do this one donor at a time!!!
small_meta1 <- pbmc_meta[pbmc_meta$donors == 's5',]
small_meta2 <- pbmc_meta[pbmc_meta$donors == 's40',]
small_meta3 <- pbmc_meta[pbmc_meta$donors == 's12',]

ct1 <- 'CD4+ T'
ct2 <- 'CD8+ T'

cells_use_1_ct1 <- sample(rownames(small_meta1[small_meta1$ctypes==ct1,]),5)
cells_use_1_ct2 <- sample(rownames(small_meta1[small_meta1$ctypes==ct2,]),5)
cells_use_2_ct1 <- sample(rownames(small_meta2[small_meta2$ctypes==ct1,]),5)
cells_use_2_ct2 <- sample(rownames(small_meta2[small_meta2$ctypes==ct2,]),5)
cells_use_3_ct1 <- sample(rownames(small_meta3[small_meta3$ctypes==ct1,]),5)
cells_use_3_ct2 <- sample(rownames(small_meta3[small_meta3$ctypes==ct2,]),5)

cells_use <- c(cells_use_1_ct1,cells_use_1_ct2,cells_use_2_ct1,cells_use_2_ct2,
               cells_use_3_ct1,cells_use_3_ct2)

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


save(test_container,file='/home/jmitchel/scITD/data/test_container.RData',compress = "xz")




















