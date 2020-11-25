
# read in the processed data
gtex_tpm_sub_transform <- readRDS(file='/home/jmitchel/data/gtex/7_tissues_counts.rds')
gtex_meta <- readRDS(file='/home/jmitchel/data/gtex/7_tissues_meta.rds')
feature.names.final <- readRDS(file='/home/jmitchel/data/gtex/genes.rds')


# basic scITD pipeline
ctypes_use <- unique(gtex_meta$ctypes)
pbmc_scMinimal <- instantiate_scMinimal(data_sparse=gtex_tpm_sub_transform,
                                        meta_data=gtex_meta)
pbmc_container <- make_new_container(pbmc_scMinimal,
                                     ctypes_use = ctypes_use,
                                     gn_convert = feature.names.final,
                                     scale_var = TRUE,
                                     var_scale_power = 1, rotate_modes = 'donors',
                                     ncores = 30, rand_seed = 10)

pbmc_container <- get_ctype_data(pbmc_container,donor_min_cells=0)

pbmc_container <- get_ctype_vargenes(pbmc_container, method="norm_var", thresh=1000)

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(10,20,7), shuffle=FALSE)

pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex'),
                                    cluster_by_meta='sex')
pbmc_container$plots$donor_matrix

pbmc_container <- optimize_var_scale_power(pbmc_container,min_ranks_test=c(2,10,7),
                                           max_ranks_test=c(7,15,7),
                                           min_power_test=0.25,max_power_test=1)
pbmc_container$plots$var_scale_plot

pbmc_container <- set_experiment_params(pbmc_container, var_scale_power = 1)

pbmc_container <- determine_ranks_tucker(pbmc_container, max_ranks_test=c(12,20,7),
                                         method='svd', shuffle_level='tensor', num_iter=5)
pbmc_container$plots$rank_determination_plot

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(9,20,7), shuffle=FALSE)

pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex'),
                                    cluster_by_meta='sex', show_donor_ids=FALSE)
pbmc_container$plots$donor_matrix

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=9,
                                      use_sig_only=F, annot='none', display_genes=F)
pbmc_container$plots$single_lds_plot

pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=9, method="fgsea",
                                      thresh=0.005, db_use="GO", num_iter=10000)
pbmc_container$plots$gsea[['Factor2']][['up']]