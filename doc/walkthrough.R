## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(scITD)

## -----------------------------------------------------------------------------
pbmc <- readRDS(file.path(find.package('scITD'),'data','PBMC_vanderwijst_clean.rds'))

## -----------------------------------------------------------------------------
feature.names = read.delim(file.path(find.package('scITD'),'data','genes.tsv'),
                           header = FALSE,
                           stringsAsFactors = FALSE)

## -----------------------------------------------------------------------------
pbmc_scMinimal <- seurat_to_scMinimal(pbmc)
pbmc_container <- make_new_container(pbmc_scMinimal,
                                     ctypes_use = c("CD4+ T", "cMonocyte", "B",
                                                  "CD8+ T", "CD56(dim) NK"),
                                     gn_convert = feature.names,
                                     scale_var = T,
                                     rotate_modes = 'donors',
                                     ncores = 30, rand_seed = 10)

## -----------------------------------------------------------------------------
pbmc_container <- get_ctype_vargenes(pbmc_container, method="empir", thresh=0.01)

## ---- fig.show = "hold", fig.width = 10, fig.height = 4-----------------------
pbmc_container <- optimize_var_scale_power(pbmc_container, max_ranks_test=c(7,15,5))
pbmc_container$plots$var_scale_plot

## -----------------------------------------------------------------------------
pbmc_container <- set_experiment_params(pbmc_container, var_scale_power = 1)

## ---- fig.show = "hold", fig.width = 8, fig.height = 8------------------------
pbmc_container <- determine_ranks_tucker(pbmc_container, max_ranks_test=c(7,15,5),
                                 method='svd', num_iter=10)
pbmc_container$plots$rank_determination_plot

## -----------------------------------------------------------------------------
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,8,5), shuffle=F)

## ---- fig.show = "hold", fig.width = 6, fig.height = 6------------------------
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex'),
                                    cluster_by_meta='sex')
pbmc_container$plots$donor_matrix

## ---- fig.show = "hold", fig.width = 10, fig.height = 6-----------------------
pbmc_container <- get_all_lds_factor_plots(pbmc_container)
render_all_lds_plots(pbmc_container$plots$all_lds_plots)

## -----------------------------------------------------------------------------
pbmc_container <- run_jackstraw(pbmc_container, n_fibers=200, n_iter=1000)


## ---- fig.show = "hold", fig.width = 6, fig.height = 6------------------------
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1,
                                        use_sig_only=T, annot='sig_genes',
                                        sig_thresh=0.05, display_genes=T)
pbmc_container$plots$single_lds_plot

## ---- fig.show = "hold", fig.width = 6, fig.height = 6------------------------
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=3,
                                        use_sig_only=T, annot='sig_genes',
                                        sig_thresh=0.05, display_genes=F)
pbmc_container$plots$single_lds_plot

## -----------------------------------------------------------------------------
gsea_res <- run_fgsea(pbmc_container, factor_select=3,
                      "cMonocyte", db_use=c("GO","Reactome"),
                      num_iter=10000)

## -----------------------------------------------------------------------------
gsea_res <- run_fgsea(pbmc_container, factor_select=3,
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
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=3,
                                      use_sig_only=T, annot='pathways',
                                      pathways=gene_sets, sig_thresh=0.05,
                                      display_genes=F)
pbmc_container$plots$single_lds_plot

