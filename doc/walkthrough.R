## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(scITD)

## -----------------------------------------------------------------------------
# load raw counts matrix
invisible(pbmc_sub_counts)

# load the normalized, log-transformed counts matrix
invisible(pbmc_sub_transformed)

# load the metadata
invisible(pbmc_sub_meta)

## -----------------------------------------------------------------------------
feature.names = read.delim(file.path(find.package('scITD'),'data','genes.tsv'),
                           header = FALSE,
                           stringsAsFactors = FALSE)

## -----------------------------------------------------------------------------
pbmc_scMinimal <- instantiate_scMinimal(pbmc_sub_transformed, pbmc_sub_counts, pbmc_sub_meta)
pbmc_scMinimal <- identify_sex_metadata(pbmc_scMinimal)
pbmc_container <- make_new_container(pbmc_scMinimal,
                                     ctypes_use = c("CD4+ T", "cMonocyte", "CD56(dim) NK"),
                                     gn_convert = feature.names,
                                     scale_var = T,
                                     rotate_modes = 'donors',
                                     ncores = 30, rand_seed = 10)

## -----------------------------------------------------------------------------
pbmc_container <- get_ctype_vargenes(pbmc_container, method="empir", thresh=0.01)

## ---- fig.show = "hold", fig.width = 8.5, fig.height = 3.5--------------------
pbmc_container <- optimize_var_scale_power(pbmc_container, max_ranks_test=c(6,10,3))
pbmc_container$plots$var_scale_plot

## -----------------------------------------------------------------------------
pbmc_container <- set_experiment_params(pbmc_container, var_scale_power = 1.25)

## ---- fig.show = "hold", fig.width = 8.5, fig.height = 5----------------------
pbmc_container <- determine_ranks_tucker(pbmc_container, max_ranks_test=c(6,10,3),
                                 method='tucker', num_iter=5)
pbmc_container$plots$rank_determination_plot

## -----------------------------------------------------------------------------
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(3,5,3), shuffle=F)

## ---- fig.show = "hold", fig.width = 6, fig.height = 6------------------------
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex'),
                                    cluster_by_meta='sex')
pbmc_container$plots$donor_matrix

## ---- fig.show = "hold", fig.width = 8, fig.height = 5------------------------
pbmc_container <- get_all_lds_factor_plots(pbmc_container)
render_all_lds_plots(pbmc_container$plots$all_lds_plots)

## -----------------------------------------------------------------------------
pbmc_container <- run_jackstraw(pbmc_container, n_fibers=200, n_iter=1000)


## ---- fig.show = "hold", fig.width = 6, fig.height = 6------------------------
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=3,
                                        use_sig_only=T, annot='sig_genes',
                                        sig_thresh=0.05, display_genes=T)
pbmc_container$plots$single_lds_plot

## ---- fig.show = "hold", fig.width = 6, fig.height = 6------------------------
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1,
                                        use_sig_only=T, annot='sig_genes',
                                        sig_thresh=0.05, display_genes=F)
pbmc_container$plots$single_lds_plot

## -----------------------------------------------------------------------------
gsea_res <- run_fgsea(pbmc_container, factor_select=1,
                      "cMonocyte", db_use=c("GO","Reactome"),
                      num_iter=10000)

## -----------------------------------------------------------------------------
gsea_res <- run_fgsea(pbmc_container, factor_select=1,
                      "CD4+ T", db_use=c("GO","Reactome"),
                      num_iter=10000)

## ---- fig.show = "hold", fig.width = 6, fig.height = 6------------------------
gene_sets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
               "GO_DEFENSE_RESPONSE_TO_VIRUS",
               "GO_CELL_CHEMOTAXIS",
               "REACTOME_INTERLEUKIN_10_SIGNALING",
               "GO_MONOCYTE_CHEMOTAXIS",
               "GO_RNA_CATABOLIC_PROCESS")

# highlight genes in pathways - limit heatmap to significant genes only
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1,
                                      use_sig_only=T, annot='pathways',
                                      pathways=gene_sets, sig_thresh=0.05,
                                      display_genes=F)
pbmc_container$plots$single_lds_plot

