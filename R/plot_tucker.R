
#' Plot matrix of donor scores extracted from Tucker decomposition
#' @importFrom circlize colorRamp2
#' @import ComplexHeatmap
#' @importFrom grid gpar
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param meta_vars character Names of metadata variables to plot alongside
#' the donor scores. Can include more than one variable. (default=NULL)
#' @param cluster_by_meta character One metadata variable to cluster the heatmap
#' by. If NULL, donor clustering is done using donor scores. (default=NULL)
#' @param show_donor_ids logical Set to TRUE to show donor id as row name on the
#' heamap (default=FALSE)
#' @param add_meta_associations character Adds meta data associations with each
#' factor as top annotation. These should be generated first with
#' plot_meta_associations(). Set to 'pval' if used 'pval' in plot_meta_associations(),
#' otherwise set to 'rsq'. If NULL, no annotation is added. (default=NULL)
#' @param show_var_explained logical Set to TRUE to display the explained variance for
#' each factor (default=TRUE)
#' @param donors_sel character A vector of a subset of donors to include in the plot
#' (default=NULL)
#'
#' @return the project container with the plot in container$plots$donor_matrix
#' @export
plot_donor_matrix <- function(container, meta_vars=NULL, cluster_by_meta=NULL,
                              show_donor_ids=FALSE, add_meta_associations=NULL,
                              show_var_explained=TRUE, donors_sel=NULL) {

  # check that Tucker has been run
  if (is.null(container$tucker_results)) {
    stop("Need to run run_tucker_ica() first.")
  }

  donor_mat <- container$tucker_results[[1]]
  donor_mat <- as.data.frame(as.matrix(donor_mat))

  # rename columns of score matrix
  colnames(donor_mat) <- sapply(1:ncol(donor_mat),function(x){
    paste0("Factor ", x)
  })

  if (show_var_explained) {
    col_fun2 = circlize::colorRamp2(c(0, max(container$exp_var)), c("white", "black"))
    ba <- HeatmapAnnotation(exp_var=container$exp_var,col = list(exp_var = col_fun2),
                            border=TRUE, show_annotation_name=FALSE)
  } else {
    ba <- NULL
  }

  if (!is.null(add_meta_associations)) {
    if (add_meta_associations=='rsq') {
      col_fun_annot = colorRamp2(c(0, 1), c("white", "forest green"))
      ta <- HeatmapAnnotation(rsq=t(container$meta_associations),col = list(rsq = col_fun_annot),
                              border=TRUE,annotation_name_side = "right")
    } else {
      col_fun_annot = colorRamp2(c(0, -log10(.05), 5), c("white", "white", "forest green"))
      logpv <- -log10(container$meta_associations)
      ta <- HeatmapAnnotation('-log_10_pval'=t(logpv),col = list('-log_10_pval'=col_fun_annot),
                              border=TRUE,annotation_name_side="right")
    }

  } else {
    ta <- NULL
  }

  # make colormap for hmap
  color_lim <- max(abs(donor_mat))
  # col_fun = colorRamp2(c(-color_lim, 0, color_lim), c("blue", "white", "red"))

  nintieth_per <- stats::quantile(as.matrix(abs(donor_mat)), c(.95))
  if (color_lim > (2*nintieth_per)) {
    col_fun = colorRamp2(c(-nintieth_per, 0, nintieth_per), c("blue", "white", "red"))
  } else {
    col_fun = colorRamp2(c(-color_lim, 0, color_lim), c("blue", "white", "red"))
  }

  if (is.null(meta_vars)) {
    if (!is.null(donors_sel)) {
      donor_mat <- donor_mat[donors_sel,]
    }
    myhmap <- Heatmap(as.matrix(donor_mat), name = "score",
                      cluster_columns = FALSE,show_column_dend = FALSE,
                      cluster_rows = TRUE, show_row_dend = FALSE,
                      column_names_gp = gpar(fontsize = 10),
                      col = col_fun, row_title = "Donors",
                      row_title_gp = gpar(fontsize = 14),
                      show_row_names = show_donor_ids,
                      border = TRUE, top_annotation=ta,
                      bottom_annotation=ba)

  } else {
    meta <- container$scMinimal_full$metadata[,c('donors',meta_vars)]
    meta <- unique(meta)
    rownames(meta) <- meta$donors
    meta$donors <- NULL

    # make all columns of meta to be factors
    for (i in 1:ncol(meta)) {
      meta[,i] <- as.factor(unlist(meta[,i]))
    }

    # limit rows of meta to those of donor_mat
    meta <- meta[rownames(donor_mat),,drop=FALSE]

    # reorder meta rows by specified meta covariate
    if (!is.null(cluster_by_meta)) {
      meta <- meta[order(meta[,cluster_by_meta]),,drop=FALSE]

      # order rows of main matrix by metadata ordering
      donor_mat <- donor_mat[rownames(meta),]

      do_row_clust <- FALSE
    } else {
      do_row_clust <- TRUE
    }

    if (!is.null(donors_sel)) {
      donor_mat <- donor_mat[donors_sel,]
      meta <- meta[donors_sel,]
    }

    myhmap <- Heatmap(as.matrix(donor_mat), name = "score",cluster_columns = FALSE,
                      cluster_rows = do_row_clust,show_row_dend = FALSE,
                      column_names_gp = gpar(fontsize = 10),
                      col = col_fun, row_title = "Donors",
                      row_title_gp = gpar(fontsize = 14),
                      show_row_names = show_donor_ids,
                      border = TRUE, top_annotation=ta,
                      bottom_annotation=ba)

    for (j in 1:ncol(meta)) {
      if (colnames(meta)[j]=='sex') {
        mycol <- RColorBrewer::brewer.pal(n = 3, name = "Accent")
        names(mycol) <- unique(meta[,j])

        myhmap <- myhmap +
          Heatmap(as.matrix(meta[,j,drop=FALSE]), name = colnames(meta)[j], cluster_rows = FALSE,
                  cluster_columns = FALSE, show_column_names = FALSE,
                  show_row_names = FALSE, col = mycol, border = TRUE)

      } else {
        myhmap <- myhmap +
          Heatmap(as.matrix(meta[,j,drop=FALSE]), name = colnames(meta)[j], cluster_rows = FALSE,
                  cluster_columns = FALSE, show_column_names = FALSE,
                  show_row_names = FALSE, border = TRUE)
      }
    }

    if (show_donor_ids) {
      myhmap <- myhmap + rowAnnotation(rn = anno_text(rownames(donor_mat)))
    }
  }

  # save plot
  container$plots$donor_matrix <- myhmap

  return(container)
}


#' Plot the gene x cell type loadings for a factor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The factor to plot
#' @param use_sig_only logical If TRUE, includes only significant genes
#' from jackstraw in the heatmap. If FALSE, includes all the variable genes.
#' (default = FALSE)
#' @param nonsig_to_zero logical If TRUE, makes the loadings of all nonsignificant genes 0
#' (default=FALSE)
#' @param annot character If set to "pathways" then creates an adjacent heatmap
#' showing which genes are in which pathways. If set to "sig_genes" then creates
#' an adjacent heatmap showing which genes were significant from jackstraw. If
#' set to "none" no adjacent heatmap is plotted. (default="none")
#' @param pathways character Gene sets to plot if annot is set to "pathways"
#' (default=NULL)
#' @param sim_de_donor_group numeric To plot the ground truth significant genes from a
#' simulation next to the heatmap, put the number of the donor group that corresponds to
#' the factor being plotted (default=NULL)
#' @param sig_thresh numeric Pvalue significance threshold to use. If use_sig_only is
#' TRUE the threshold is used as a cutoff for genes to include. If annot is "sig_genes"
#' this value is used in the gene significance colormap as a minimum threshold. (default=0.05)
#' @param display_genes logical If TRUE, displays the names of gene names (default=FALSE)
#' @param gene_callouts logical If TRUE, then adds gene callout annotations to the heatmap
#' (default=FALSE)
#' @param callout_n_gene_per_ctype numeric To use if gene_callouts is TRUE. Sets the number
#' of largest magnitude significant genes from each cell type to include in gene callouts.
#' (default=5)
#' @param callout_ctypes character To use if gene_callouts is TRUE. Specifies which cell types
#' to get gene callouts for. If NULL, then gets gene callouts for largest magnitude significant
#' genes for all cell types. (default=NULL)
#' @param le_set_callouts character Pass a vector of gene set names to show leading edge genes
#' for a select set of gene sets (default=NULL)
#' @param le_set_colormap character A named vector with names as gene sets and values as colors.
#' If NULL, then selects first n colors of Set3 color palette. (default=NULL)
#' @param le_set_num_per numeric The number of leading edge genes to show for each gene set (default=5)
#' @param show_le_legend logical Set to TRUE to show the color map legend for leading edge genes (default=FALSE)
#' @param show_xlab logical If TRUE, displays the xlabel 'genes' (default=TRUE)
#' @param show_var_explained logical If TRUE then shows an anotation with the explained variance
#' for each cell type (default=TRUE)
#' @param reset_other_factor_plots logical Set to TRUE to set all other loadings plots to NULL.
#' Useful if run get_all_lds_factor_plots but then only want to show one or two plots. (default=FALSE)
#' @param draw_plot logical Set to TRUE to show the plot. Plot is stored regardless. (default=TRUE)
#'
#' @return container with the plot put in container$plots$all_lds_plots and the legend put in
#' container$plots$all_legends. Use draw(<hmap obj>,annotation_legend_list = <hmap legend obj>)
#' to re-render the plot with legend
#' @export
plot_loadings_annot <- function(container, factor_select, use_sig_only=FALSE, nonsig_to_zero=FALSE, annot='none',
                                pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL,
                                le_set_callouts=NULL, le_set_colormap=NULL, le_set_num_per=5, show_le_legend=FALSE,
                                show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE,
                                draw_plot=TRUE) {

  # check that Tucker has been run
  if (is.null(container$tucker_results)) {
    stop("Need to run run_tucker_ica() first.")
  }

  # remove other loadings plots if indicated
  if (reset_other_factor_plots) {
    container$plots$all_lds_plots <- NULL
    container$plots$all_legends <- NULL
  }

  ldngs <- container$tucker_results[[2]]

  # break down a factor from the loadings matrix
  genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
  ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})

  sr_col <- ldngs[factor_select,]

  tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)

  # limit loadings df to just sig genes if specified
  if (use_sig_only) {
    # make sure jackstraw has been run
    if (is.null(container$gene_score_associations)) {
      stop('Run jackstraw first to display significant genes')
    }

    sig_vectors <- get_significance_vectors(container,
                                            factor_select, colnames(tmp_casted_num))
    # convert list to df
    sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

    # order df same way as in tmp_casted_num
    sig_df <- sig_df[rownames(tmp_casted_num),colnames(tmp_casted_num)]

    # reduce tmp_casted_num to genes significant in at least one cell type
    tmp_casted_num <- tmp_casted_num[rowSums(sig_df < sig_thresh) > 0,]

    if (nonsig_to_zero) {
      tmp_casted_num[sig_df[rownames(tmp_casted_num),colnames(tmp_casted_num)] > sig_thresh] <- 0
    }
  }

  if (show_xlab) {
    rt <- "Genes"
  } else {
    rt <- ""
  }

  if (gene_callouts) {
    gene_callouts <- get_callouts_annot(container, tmp_casted_num, factor_select, sig_thresh,
                       top_n_per_ctype=callout_n_gene_per_ctype, ctypes=callout_ctypes)
  } else {
    gene_callouts <- NULL
  }

  if (!is.null(le_set_callouts)) {
    # get leading edge genes to plot
    le_genes <- get_leading_edge_genes(container, factor_select, gsets=le_set_callouts,
                           num_genes_per=le_set_num_per)

    # get colors for each gene by its gene set
    if (is.null(le_set_colormap)) { # need to pick random colors if not specified
      le_set_colormap <- RColorBrewer::brewer.pal(n = length(le_set_callouts), name = "Set3")
      names(le_set_colormap) <- le_set_callouts
    }
    le_colors <- c()
    for (i in 1:length(le_genes)) {
      gs <- le_genes[i]
      mycolor <- le_set_colormap[gs]
      le_colors[i] <- mycolor
    }

    # get indices for each gene
    le_ndx <- match(names(le_genes),rownames(tmp_casted_num))

    gene_callouts <- rowAnnotation(callouts = anno_mark(at = le_ndx, which='row',
                                                        labels = names(le_genes),
                                                        labels_gp = gpar(col = le_colors, fontsize = 11),
                                                        link_gp = gpar(lwd=1.25, col = le_colors),
                                                        padding = unit(.75, "mm")))

    # make legend if specified to do so
    if (show_le_legend) {
      le_legend <- Legend(labels = names(le_set_colormap),
                          legend_gp = gpar(fill = le_set_colormap), title = "gene sets",
             grid_height = unit(1, "mm"), grid_width = unit(3, "mm"))
    }


  } else {
    gene_callouts <- NULL
  }

  hm_legends <- list()

  if (show_var_explained) {

    explained_variances <- c()
    for (i in 1:ncol(tmp_casted_num)) {
      ct<- colnames(tmp_casted_num)[i]
      exp_var <- get_ctype_exp_var(container,factor_select,ct)
      explained_variances[i] <- exp_var
    }
    col_fun2 = circlize::colorRamp2(c(0, max(explained_variances)), c("white", "black"))
    var_annot <- ComplexHeatmap::HeatmapAnnotation(exp_var = explained_variances,col=list(exp_var=col_fun2),
                                                   show_annotation_name=FALSE, border=TRUE,
                                                   show_legend = FALSE)
    hm_legends[[2]] <- Legend(col_fun = col_fun2, title = "var exp",
                              grid_height = unit(1, "mm"), grid_width = unit(3, "mm"),
                              title_position = "leftcenter-rot")

  } else {
    var_annot <- NULL
  }

  # make colormap for hmap
  # color_lim <- max(abs(tmp_casted_num))
  # nintieth_per <- stats::quantile(as.matrix(abs(tmp_casted_num)), c(.95))
  # if (color_lim > (1.5*nintieth_per)) {
  #   col_fun = colorRamp2(c(-nintieth_per, 0, nintieth_per), c("blue", "white", "red"))
  # } else {
  #   col_fun = colorRamp2(c(-color_lim, 0, color_lim), c("blue", "white", "red"))
  # }

  # color_lim <- stats::quantile(as.matrix(abs(tmp_casted_num)), c(.9999))
  color_lim <- stats::quantile(as.matrix(abs(tmp_casted_num)), c(.99))
  col_fun = colorRamp2(c(-color_lim, 0, color_lim), c("blue", "white", "red"))


  hm_legends[[1]] <- Legend(col_fun = col_fun, title = "loading",
                            grid_height = unit(1, "mm"), grid_width = unit(3, "mm"),
                            title_position = "leftcenter-rot")

  # 'median' clustering method works well
  hm_list <- Heatmap(tmp_casted_num, show_row_dend = FALSE, show_column_dend = FALSE,
                     name = "loading", show_row_names = display_genes,
                     column_names_gp = gpar(fontsize = 12), cluster_columns = FALSE,
                     clustering_method_rows = "median",
                     row_names_side = "left", col=col_fun,
                     column_title = paste0('Factor ', factor_select),
                     column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                     row_title = rt, row_title_gp = gpar(fontsize = 14), border = TRUE,
                     row_labels = convert_gn(container,rownames(tmp_casted_num)),
                     right_annotation = gene_callouts, top_annotation=var_annot,
                     show_heatmap_legend = FALSE,
                     width = unit(10, "cm"),
                     height = unit(20, "cm")) #used to use w=8, h=14. or 10, 18

  # turn off heatmap message saying callouts require pdf view or zoom view
  ht_opt$message = FALSE

  if (annot == 'pathways') {
    if (is.null(container$gn_convert)) {
      stop('Gene symbols are not present in your data and no gene name conversion was provided')
    }

    ### fix this to be same way as sig_genes so can work with reduced df
    gene_set_vectors <- get_gene_set_vectors(container, pathways, tmp_casted_num)

    for (i in 1:length(gene_set_vectors)) {
      g_vec <- gene_set_vectors[[i]]
      hm_list <- hm_list +
        Heatmap(g_vec + 0, name = pathways[i],
                col = c("0" = "white", "1" = "black"),
                show_heatmap_legend = FALSE, width = unit(5, "mm"),
                column_names_gp = gpar(fontsize = 7), border = TRUE)
    }
  } else if (annot == 'sig_genes') {
    # make sure jackstraw has been run
    if (is.null(container$gene_score_associations)) {
      stop('Run jackstraw first to display significant genes')
    }

    sig_vectors <- get_significance_vectors(container,
                                            factor_select, colnames(tmp_casted_num))
    # convert list to df
    sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

    # limit to just the genes in tmp_casted_num
    sig_df <- sig_df[rownames(tmp_casted_num),colnames(tmp_casted_num)]

    # add gene significance heatmap to total hmap
    col_fun = colorRamp2(c(sig_thresh, 0), c("white", "green"))
    hm_list <- hm_list +
      Heatmap(sig_df,
              name = "adj p-value", cluster_columns = FALSE,
              col = col_fun, show_row_names = FALSE,
              show_heatmap_legend = TRUE, show_column_dend = FALSE,
              column_names_gp = gpar(fontsize = 12), border = TRUE)
  }

  if (!is.null(sim_de_donor_group)) {

    de_genes <- sim_de_donor_group[[factor_select]][rownames(tmp_casted_num),]
    color_lim <- stats::quantile(as.matrix(abs(de_genes)), c(.95))
    col_fun = colorRamp2(c(0,0.6,color_lim), c("white","white", "violet"))
    hm_list <- hm_list +
      Heatmap(as.matrix(de_genes), col=col_fun,
              name = "True DE Genes", cluster_columns = FALSE,
              show_row_names = FALSE,
              show_heatmap_legend = TRUE, show_column_dend = FALSE,
              column_names_gp = gpar(fontsize = 12), border = TRUE)


  }

  # save plot in the container
  container$plots$all_lds_plots[[as.character(factor_select)]] <- hm_list

  # store matrix that generated the plot
  container$plots$lds_plots_data[[as.character(factor_select)]] <- tmp_casted_num

  # pack and save legend in container
  if (show_var_explained) {
    pd <- packLegend(hm_legends[[1]], hm_legends[[2]], direction = "vertical")
  } else {
    pd <- hm_legends[[1]]
  }
  container$plots$all_legends[[as.character(factor_select)]] <- pd

  # optionally draw the plot
  if (draw_plot) {
    if (show_le_legend) {
      draw(hm_list,annotation_legend_list = pd,
           legend_grouping = "original",
           heatmap_legend_list = le_legend, heatmap_legend_side = "bottom",
           newpage=FALSE)
    } else {
      draw(hm_list,annotation_legend_list = pd,
           legend_grouping = "original",
           newpage=FALSE)
    }

  }

  return(container)
}

#' Reshape loadings for a factor from linearized to matrix form
#'
#' @param ldngs_row numeric A vector of loadings values for one factor
#' @param genes character The gene identifiers corresponding to each loading
#' @param ctypes character The cell type corresponding to each loading
#'
#' @return a loadings matrix with dimensions of genes by cell types
#' @export
reshape_loadings <- function(ldngs_row,genes,ctypes) {
  # create df with genes ctype and value
  tmp <- cbind(ctypes,genes,ldngs_row)

  # transform to gene by cell type matrix
  tmp_casted <- reshape2::dcast(as.data.frame(tmp),
                      genes~ctypes, value.var='ldngs_row')
  rownames(tmp_casted) <- tmp_casted$genes
  tmp_casted$genes <- NULL

  # remove any rows with NA
  tmp_casted[,"NA"] <- NULL
  tmp_casted <- tmp_casted[rowSums(is.na(tmp_casted)) == 0, ]

  # convert all columns to numeric
  tmp_casted_num <- sapply(tmp_casted, as.numeric)
  rownames(tmp_casted_num) <- rownames(tmp_casted)
  return(tmp_casted_num)
}

#' Get logical vectors indicating which genes are in which pathways
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param gene_sets character Vector of gene sets to extract genes for
#' @param tmp_casted_num matrix The gene by cell type loadings matrix
#'
#' @return list of the logical vectors for each pathway
#' @export
get_gene_set_vectors <- function(container,gene_sets,tmp_casted_num) {
  m_df = msigdbr::msigdbr(species = "Homo sapiens")
  my_pathways = split(m_df$gene_symbol, f = m_df$gs_name)

  gene_set_vectors <- lapply(gene_sets,function(x) {
    set_x <- my_pathways[[x]]
    gs_indices <- convert_gn(container, rownames(tmp_casted_num)) %in% set_x
  })

  return(gene_set_vectors)
}

#' Get vectors indicating which genes are significant in which cell types
#' for a factor of interest
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The factor to query
#' @param ctypes character The cell types used in all the analysis ordered
#' as they appear in the loadings matrix
#'
#' @return a list of pvalues for each gene in each cell type
#' @export
get_significance_vectors <- function(container, factor_select, ctypes) {
  # parse the gene significance results to get only gene_ctype combos for factor of interest
  padj <- container$gene_score_associations
  padj_factors <- sapply(names(padj),function(x) {
    tmp <- strsplit(x,split = '.', fixed = TRUE)[[1]]
    return(tmp[[length(tmp)]])
  })
  padj_use <- padj[which(padj_factors == as.character(factor_select))]

  # get out pvals for gene_ctype combos of specific ctypes
  padj_all_ctypes <- list()
  for (ct in ctypes) {
    padj_ct <- sapply(names(padj_use),function(x) {
      tmp <- strsplit(x,split = '.', fixed = TRUE)[[1]]
      return(tmp[[length(tmp)-1]])
    })
    padj_ct <- padj_use[which(padj_ct == ct)]

    names(padj_ct) <- sapply(names(padj_ct),function(x) {
      tmp <- strsplit(x,split = '.',fixed = TRUE)[[1]]

      if (length(tmp)>3){
        return(paste0(tmp[[1]],".",tmp[[2]]))
      } else {
        return(tmp[[1]])
      }
    })

    padj_all_ctypes[[ct]] <- padj_ct
  }
  return(padj_all_ctypes)
}


#' Get gene callouts annotation for a loadings heatmap
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param tmp_casted_num matrix The gene by cell type loadings matrix
#' @param factor_select numeric The factor to investigate
#' @param sig_thresh numeric Pvalue cutoff for significant genes
#' @param top_n_per_ctype numeric The number of significant, largest magnitude
#' genes from each cell type to generate callouts for (default=5)
#' @param ctypes character The cell types for which to get the top genes to make
#' callouts for. If NULL then uses all cell types. (default=NULL)
#'
#' @return HeatmapAnnotation for the gene callouts
get_callouts_annot <- function(container, tmp_casted_num, factor_select, sig_thresh, top_n_per_ctype=5, ctypes=NULL) {

  # extract the genes to show
  if (is.null(ctypes)) {
    ctypes <- container$experiment_params$ctypes_use
  }
  sig_vecs <- get_significance_vectors(container,factor_select,ctypes)
  genes_plot <- c()
  for (ct in ctypes) {
    # get significant genes for the ctype
    ct_sig_genes <- sig_vecs[[ct]]
    ct_sig_genes <- ct_sig_genes[ct_sig_genes < sig_thresh]

    # get top loading genes of the significant ones
    ct_sig_loadings <- tmp_casted_num[names(ct_sig_genes),ct]

    ct_sig_loadings <- ct_sig_loadings[order(abs(ct_sig_loadings),decreasing=TRUE)]
    ct_top_genes <- names(ct_sig_loadings)[1:top_n_per_ctype]
    genes_plot <- c(genes_plot,ct_top_genes)
  }

  gene_callouts <- unique(genes_plot)

  ndx <- match(gene_callouts,rownames(tmp_casted_num))
  callouts <- list()
  callouts[[1]] <- ndx
  callouts[[2]] <- convert_gn(container, gene_callouts)

  myannot <- rowAnnotation(callouts = anno_mark(at = callouts[[1]], which='row',
                                                labels = callouts[[2]]))
  return(myannot)
}


#' Generate loadings heatmaps for all factors
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param use_sig_only logical If TRUE, includes only significant genes
#' from jackstraw in the heatmap. If FALSE, includes all the variable genes.
#' (default = FALSE)
#' @param nonsig_to_zero logical If TRUE, makes the loadings of all nonsignificant genes 0
#' (default=FALSE)
#' @param annot character If set to "pathways" then creates an adjacent heatmap
#' showing which genes are in which pathways. If set to "sig_genes" then creates
#' an adjacent heatmap showing which genes were significant from jackstraw. If
#' set to "none" no adjacent heatmap is plotted. (default="none")
#' @param pathways_list list A list of sets of pathways for each factor. List index
#' should be the number corresponding to the factor. (default=NULL)
#' @param sim_de_donor_group numeric To plot the ground truth significant genes from a
#' simulation next to the heatmap, put the number of the donor group that corresponds to
#' the factor being plotted. Here it should be a vector corresponding to the factors.
#' (default=NULL)
#' @param sig_thresh numeric Pvalue significance threshold to use. If use_sig_only is
#' TRUE the threshold is used as a cutoff for genes to include. If annot is "sig_genes"
#' this value is used in the gene significance colormap as a minimum threshold. (default=0.05)
#' @param display_genes logical If TRUE, displays the names of gene names (default=FALSE)
#' @param gene_callouts logical If TRUE, then adds gene callout annotations to the heatmap
#' (default=FALSE)
#' @param callout_n_gene_per_ctype numeric To use if gene_callouts is TRUE. Sets the number
#' of largest magnitude significant genes from each cell type to include in gene callouts.
#' (default=5)
#' @param callout_ctypes list To use if gene_callouts is TRUE. Specifies which cell types
#' to get gene callouts for. Each entry of the list should be a character vector of ctypes for
#' the respective factor. If NULL, then gets gene callouts for largest magnitude significant
#' genes for all cell types. (default=NULL)
#' @param show_var_explained logical If TRUE then shows an anottation with the explained variance
#' for each cell type (default=TRUE)
#'
#' @return the project container with the list of plots placed in container$plots$all_lds_plots
#' @export
get_all_lds_factor_plots <- function(container, use_sig_only=FALSE, nonsig_to_zero=FALSE, annot='none',
                                     pathways_list=NULL, sim_de_donor_group=NULL,
                                     sig_thresh=0.05, display_genes=FALSE,
                                     gene_callouts=FALSE, callout_n_gene_per_ctype=5,
                                     callout_ctypes=NULL,
                                     show_var_explained=TRUE) {

  num_fact <- nrow(container$tucker_results[[2]])
  for (i in 1:num_fact) {
    container <- plot_loadings_annot(container, factor_select=i,
                                     use_sig_only=use_sig_only,
                                     nonsig_to_zero=nonsig_to_zero,
                                     annot=annot, pathways=pathways_list[[i]],
                                     sig_thresh=sig_thresh,
                                     sim_de_donor_group=sim_de_donor_group[i],
                                     display_genes=display_genes,
                                     gene_callouts=gene_callouts,
                                     callout_n_gene_per_ctype=callout_n_gene_per_ctype,
                                     callout_ctypes=callout_ctypes[[i]],
                                     show_xlab=TRUE,
                                     show_var_explained=show_var_explained,
                                     reset_other_factor_plots=FALSE,
                                     draw_plot=FALSE)

  }

  return(container)
}

#' Creates a figure of all loadings plots arranged
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param data_type character Can be either "loadings", "gsea", or "dgenes". This
#' determines which list of heatmaps to organize into the figure.
#' @param max_cols numeric The max number of columns to plot. Can only either be 2
#' or 3 since these are large plots. (default=3)
#'
#' @return the multi-plot figure
#' @export
render_multi_plots <- function(container,data_type,max_cols=3) {
  # if (!(max_cols == 2 || max_cols == 3)) {
  #   stop('max_cols can only be set to 2 or 3')
  # }

  if (data_type == "loadings") {
    hm_list <- container$plots$all_lds_plots
    hm_legends <- container$plots$all_legends
  } else if (data_type == "gsea") {
    hm_list <- container$plots$gsea
  } else if (data_type == "dgenes") {
    hm_list <- container$plots$donor_sig_genes
  }

  # order the list of heatmaps by factor number
  hm_order <- order(as.numeric(names(hm_list)),decreasing=FALSE)
  hm_list <- hm_list[hm_order]
  if (data_type=='ldngs') {
    hm_legends <- hm_legends[hm_order]
  }

  grob_lst <- list()
  for (i in 1:length(hm_list)) {
    if (data_type=='loadings') {
      gb <- grid::grid.grabExpr(draw(hm_list[[i]],annotation_legend_list = hm_legends[[i]],
                                     legend_grouping = "original",
                                     newpage=FALSE))
    } else {
      gb <- grid::grid.grabExpr(draw(hm_list[[i]], newpage=FALSE))
    }

    grob_lst[[i]] <- gb
  }

  num_plots <- length(grob_lst)

  if (num_plots > max_cols) {
    num_rows <- floor(num_plots/max_cols)
    num_bottom <- num_plots %% max_cols

    top_rows <- grob_lst[1:(num_plots-num_bottom)]
    top_rows <- cowplot::plot_grid(plotlist=top_rows,ncol=max_cols,align = "v")

    if (num_bottom==1) {
      bottom_row <- list(NULL,grob_lst[[num_plots]],NULL)
      bottom_row <- cowplot::plot_grid(plotlist=bottom_row, ncol=max_cols,rel_widths=c((1/2.825),(1/3),(1/3)))
      fig <- cowplot::plot_grid(top_rows, bottom_row, ncol=1, rel_heights=c(num_rows,1),align = "v")
    } else if (num_bottom==2) {
      bottom_row <- list(NULL,grob_lst[[num_plots-1]],grob_lst[[num_plots]],NULL)
      bottom_row <- cowplot::plot_grid(plotlist=bottom_row, ncol=4, rel_widths=c((1/3)/2,(1/3),(1/3),(1/3)/2))
      fig <- cowplot::plot_grid(top_rows, bottom_row, ncol=1, rel_heights=c(num_rows,1),align = "v")
    } else {
      fig <- top_rows
    }
  } else {
    fig <- cowplot::plot_grid(plotlist=grob_lst, ncol=num_plots)
  }


  return(fig)
}


#' Generate heatmap showing top genes in each cell type significantly associated
#' with a given factor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The factor to query
#' @param top_n_per_ctype numeric Vector of the number of top genes from each cell type
#' to plot
#' @param ctypes_use character The cell types for which to get the top genes to make
#' callouts for. If NULL then uses all cell types. (default=NULL)
#' @param show_donor_labels logical Set to TRUE to display donor labels (default=FALSE)
#' @param additional_meta character Another meta variable to plot (default=NULL)
#' @param add_genes character Additional genes to plot for all ctypes (default=NULL)
#'
#' @return the project container with the plot in the slot
#' container$plots$donor_sig_genes$Factor#
#' @export
plot_donor_sig_genes <- function(container, factor_select, top_n_per_ctype,
                                 ctypes_use=NULL, show_donor_labels=FALSE,
                                 additional_meta=NULL, add_genes=NULL) {
  ## add catch in case they havent run jackstraw yet...

  # extract tensor information
  tensor_data <- container$tensor_data
  donor_nm <- tensor_data[[1]]
  gene_nm  <- tensor_data[[2]]
  ctype_nm  <- tensor_data[[3]]
  tnsr <- tensor_data[[4]]

  # get the loadings matrix
  ldngs <- container$tucker_results[[2]]

  # break down a factor from the loadings matrix
  genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
  ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})

  sr_col <- ldngs[factor_select,]

  tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)

  # extract the genes to show
  if (is.null(ctypes_use)) {
    ctypes <- container$experiment_params$ctypes_use
  } else {
    ctypes <- ctypes_use
  }
  sig_vecs <- get_significance_vectors(container,factor_select,ctypes)
  genes_plot <- c()
  ct_in_hmap <- c()
  for (i in 1:length(ctypes)) {
    ct <- ctypes[i]
    if (length(top_n_per_ctype)==1) {
      top_n <- top_n_per_ctype
    } else {
      top_n <- top_n_per_ctype[i]
    }

    # get significant genes for the ctype
    ct_sig_genes <- sig_vecs[[ct]]
    ct_sig_genes <- ct_sig_genes[ct_sig_genes<0.05]

    # get top loading genes of the significant ones
    ct_sig_loadings <- tmp_casted_num[names(ct_sig_genes),ct]

    ct_sig_loadings <- ct_sig_loadings[order(abs(ct_sig_loadings),decreasing=TRUE)]
    ct_top_genes <- names(ct_sig_loadings)[1:top_n]
    if (!is.null(add_genes)) {
      ct_top_genes <- unique(c(ct_top_genes,add_genes))
    }
    ct_top_genes <- sapply(ct_top_genes,function(x) {paste0(x,"_",ct)})
    genes_plot <- c(genes_plot,ct_top_genes)
    ct_in_hmap <- c(ct_in_hmap, rep(ct,top_n))
  }

  ct_in_hmap <- factor(ct_in_hmap)

  # unfold tensor along donor mode
  donor_unfold <- rTensor::k_unfold(rTensor::as.tensor(tnsr),1)@data

  gn_ctype_cnames <- c()
  for (ct in ctype_nm) {
    for (gn in gene_nm) {
      gn_ctype_cnames <- c(gn_ctype_cnames,paste0(gn,"_",ct))
    }
  }

  colnames(donor_unfold) <- gn_ctype_cnames
  rownames(donor_unfold) <- donor_nm

  ## testing out scaling the data to unit variance
  donor_unfold <- scale(donor_unfold)

  # subset data to just genes to plot
  donor_unfold_sub <- donor_unfold[,genes_plot]
  donor_unfold_sub <- t(donor_unfold_sub)

  # reorder donors by their score for the factor
  donor_scores <- container$tucker_results[[1]]
  donor_scores <- donor_scores[,factor_select]
  donor_unfold_sub <- donor_unfold_sub[,order(donor_scores)]
  donor_scores <- donor_scores[order(donor_scores)]

  donor_scores <- unlist(donor_scores)
  col_fun2 = circlize::colorRamp2(c(min(donor_scores), 0, max(donor_scores)), c("purple", "white", "green"))
  ha <- ComplexHeatmap::HeatmapAnnotation(score = donor_scores,col=list(score=col_fun2),
                                          show_annotation_name=FALSE)


  if (!is.null(additional_meta)) {
    meta <- container$scMinimal_full$metadata[,c('donors',additional_meta)]
    meta <- unique(meta)
    rownames(meta) <- meta$donors
    meta$donors <- NULL
    meta <- meta[colnames(donor_unfold_sub),,drop=FALSE]

    # make all columns of meta to be factors
    for (i in 1:ncol(meta)) {
      meta[,i] <- factor(unlist(meta[,i]),levels=unique(unlist(meta[,i]))[order(unique(unlist(meta[,i])))])
    }

    set.seed(30)
    if (length(levels(meta)) < 3) {
      mycol <- RColorBrewer::brewer.pal(n = 3, name = "Paired")
    } else {
      mycol <- RColorBrewer::brewer.pal(n = length(levels(meta)), name = "Paired")
    }
    names(mycol) <- levels(meta)
    ta <- ComplexHeatmap::HeatmapAnnotation(df = meta, show_annotation_name=TRUE,
                                            col = list(df = mycol))
  } else {
    ta <- NULL
  }



  # rename genes
  rownames(donor_unfold_sub) <- sapply(rownames(donor_unfold_sub),function(x) {
    gn <- strsplit(x,split="_")[[1]][1]
    ct <- strsplit(x,split="_")[[1]][2]
    gn <- convert_gn(container,gn)
    return(paste0(gn,"_",ct))
  })

  rn_show <- sapply(rownames(donor_unfold_sub),function(x){
    strsplit(x,split="_")[[1]][[1]]
  })
  ct_show <- sapply(rownames(donor_unfold_sub),function(x){
    strsplit(x,split="_")[[1]][[2]]
  })
  ct_show <- factor(ct_show,levels=ctypes_use)

  set.seed(10)
  mycol <- RColorBrewer::brewer.pal(n = length(ctypes_use), name = "Accent")
  names(mycol) <- ctypes_use

  ct_annot <- ComplexHeatmap::rowAnnotation(cell_types=anno_simple(ct_show),
                                            show_annotation_name=FALSE,
                                            col = list(cell_types = mycol))

  # create the hmap
  col_fun = colorRamp2(c(min(donor_unfold_sub), 0, max(donor_unfold_sub)), c("blue", "white", "red"))

  myhmap <- Heatmap(donor_unfold_sub, name = "expr",
                    cluster_columns = FALSE,
                    cluster_rows = TRUE,
                    cluster_row_slices=FALSE,
                    column_names_gp = gpar(fontsize = 8),
                    row_names_gp = gpar(fontsize = 10),
                    col = col_fun, bottom_annotation=ha, row_split = ct_show,
                    row_labels=rn_show,border=TRUE, show_column_names=show_donor_labels,
                    left_annotation=ct_annot, show_row_dend = FALSE,
                    column_title = paste0('Factor ',as.character(factor_select)),
                    column_title_gp = gpar(fontsize = 20),
                    column_title_side = "top",
                    top_annotation=ta)

  container$plots$donor_sig_genes[[as.character(factor_select)]] <- myhmap
  return(container)
}


#' Generate heatmap showing top genes in each cell type significantly associated
#' with a given factor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The factor to query
#' @param top_n_per_ctype numeric Vector of the number of top genes from each cell type
#' to plot
#' @param ctypes_use character The cell types for which to get the top genes to make
#' callouts for. If NULL then uses all cell types. (default=NULL)
#' @param show_donor_labels logical Set to TRUE to display donor labels (default=FALSE)
#' @param additional_meta character Another meta variable to plot (default=NULL)
#'
#' @return the project container with the plot in the slot
#' container$plots$donor_sig_genes$Factor#
#' @export
plot_donor_sig_genes_v2 <- function(container, factor_select, top_n_per_ctype,
                                 ctypes_use=NULL, show_donor_labels=FALSE,
                                 additional_meta=NULL) {
  ## add catch in case they havent run jackstraw yet...

  # extract tensor information
  tensor_data <- container$tensor_data
  donor_nm <- tensor_data[[1]]
  gene_nm  <- tensor_data[[2]]
  ctype_nm  <- tensor_data[[3]]
  tnsr <- tensor_data[[4]]

  # get the loadings matrix
  ldngs <- container$tucker_results[[2]]

  # break down a factor from the loadings matrix
  genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
  ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})

  sr_col <- ldngs[factor_select,]

  tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)

  # extract the genes to show
  if (is.null(ctypes_use)) {
    ctypes <- container$experiment_params$ctypes_use
  } else {
    ctypes <- ctypes_use
  }
  sig_vecs <- get_significance_vectors(container,factor_select,ctypes)
  genes_plot <- c()
  ct_in_hmap <- c()
  for (i in 1:length(ctypes)) {
    ct <- ctypes[i]
    if (length(top_n_per_ctype)==1) {
      top_n <- top_n_per_ctype
    } else {
      top_n <- top_n_per_ctype[i]
    }

    # get significant genes for the ctype
    ct_sig_genes <- sig_vecs[[ct]]
    ct_sig_genes <- ct_sig_genes[ct_sig_genes<0.05]

    # get top loading genes of the significant ones
    ct_sig_loadings <- tmp_casted_num[names(ct_sig_genes),ct]

    # ## testing not using gene significance
    # ct_sig_loadings <- tmp_casted_num[,ct]

    ct_sig_loadings <- ct_sig_loadings[order(abs(ct_sig_loadings),decreasing=TRUE)]
    ct_top_genes <- names(ct_sig_loadings)[1:top_n]
    ct_top_genes <- sapply(ct_top_genes,function(x) {paste0(x,"_",ct)})
    genes_plot <- c(genes_plot,ct_top_genes)
    ct_in_hmap <- c(ct_in_hmap, rep(ct,top_n))
  }

  ct_in_hmap <- factor(ct_in_hmap)

  # unfold tensor along donor mode
  donor_unfold <- rTensor::k_unfold(rTensor::as.tensor(tnsr),1)@data

  gn_ctype_cnames <- c()
  for (ct in ctype_nm) {
    for (gn in gene_nm) {
      gn_ctype_cnames <- c(gn_ctype_cnames,paste0(gn,"_",ct))
    }
  }

  colnames(donor_unfold) <- gn_ctype_cnames
  rownames(donor_unfold) <- donor_nm

  ## testing out scaling the data to unit variance
  donor_unfold <- scale(donor_unfold)

  # subset data to just genes to plot
  donor_unfold_sub <- donor_unfold[,genes_plot]
  donor_unfold_sub <- t(donor_unfold_sub)

  # reorder donors by their score for the factor
  donor_scores <- container$tucker_results[[1]]
  donor_scores_fact <- donor_scores[,factor_select]
  dsc_ord <- order(donor_scores_fact)
  donor_unfold_sub <- donor_unfold_sub[,dsc_ord]
  donor_scores <- donor_scores[dsc_ord,]

  colnames(donor_scores) <- sapply(1:ncol(donor_scores), function(x) {
    paste0('Factor',as.character(x))
  })
  col_fun2 = circlize::colorRamp2(c(min(donor_scores), 0, max(donor_scores)), c("purple", "white", "green"))
  dscores_hmap <- Heatmap(as.matrix(t(donor_scores)),name='dscores',
                          cluster_columns = FALSE,
                          cluster_rows = FALSE,
                          show_row_names = TRUE,
                          show_column_names = FALSE,
                          col = col_fun2, border = TRUE,
                          height = unit(3.5, "cm"),
                          column_title = 'Donors',
                          column_title_side = "bottom")


  if (!is.null(additional_meta)) {
    meta <- container$scMinimal_full$metadata[,c('donors',additional_meta)]
    meta <- unique(meta)
    rownames(meta) <- meta$donors
    meta$donors <- NULL
    meta <- meta[colnames(donor_unfold_sub),,drop=FALSE]

    # make all columns of meta to be factors
    for (i in 1:ncol(meta)) {
      meta[,i] <- as.factor(unlist(meta[,i]))
    }

    set.seed(30)
    if (length(levels(meta)) < 3) {
      mycol <- RColorBrewer::brewer.pal(n = 3, name = "Paired")
    } else {
      mycol <- RColorBrewer::brewer.pal(n = length(levels(meta)), name = "Paired")
    }
    names(mycol) <- levels(meta)
    ta <- ComplexHeatmap::HeatmapAnnotation(df = meta, show_annotation_name=TRUE,
                                            col = list(df = mycol))

  } else {
    ta <- NULL
  }



  # rename genes
  rownames(donor_unfold_sub) <- sapply(rownames(donor_unfold_sub),function(x) {
    gn <- strsplit(x,split="_")[[1]][1]
    ct <- strsplit(x,split="_")[[1]][2]
    gn <- convert_gn(container,gn)
    return(paste0(gn,"_",ct))
  })

  rn_show <- sapply(rownames(donor_unfold_sub),function(x){
    strsplit(x,split="_")[[1]][[1]]
  })
  ct_show <- sapply(rownames(donor_unfold_sub),function(x){
    strsplit(x,split="_")[[1]][[2]]
  })
  ct_show <- factor(ct_show,levels=ctypes_use)

  set.seed(10)
  mycol <- RColorBrewer::brewer.pal(n = length(ctypes_use), name = "Accent")
  names(mycol) <- ctypes_use

  ct_annot <- ComplexHeatmap::rowAnnotation(cell_types=anno_simple(ct_show),
                                            show_annotation_name=FALSE,
                                            col = list(cell_types = mycol))

  # create the hmap
  col_fun = colorRamp2(c(min(donor_unfold_sub), 0, max(donor_unfold_sub)), c("blue", "white", "red"))

  myhmap <- Heatmap(donor_unfold_sub, name = "expr",
                    cluster_columns = FALSE,
                    cluster_rows = TRUE,
                    column_names_gp = gpar(fontsize = 8),
                    row_names_gp = gpar(fontsize = 10),
                    col = col_fun, row_split = ct_show,
                    row_labels=rn_show,border=TRUE, show_column_names=show_donor_labels,
                    left_annotation=ct_annot,show_row_dend = FALSE,
                    top_annotation=ta,
                    column_title = paste0('Factor ',as.character(factor_select)),
                    column_title_gp = gpar(fontsize = 20),
                    column_title_side = "top")

  myhmap <- myhmap %v% dscores_hmap

  container$plots$donor_sig_genes[[as.character(factor_select)]] <- myhmap
  return(container)
}



#' Pairwise comparison of factors from two separate decompositions
#'
#' @param tucker_res1 list The container$tucker_res from first decomposition
#' @param tucker_res2 list The container$tucker_res from first decomposition
#' @param decomp_names character Names of the two decompositions that will go
#' on the axes of the heatmap
#' @param meta_anno1 matrix The result of calling get_meta_associations()
#' corresponding to the first decomposition, which is stored in
#' container$meta_associations
#' @param meta_anno2 matrix The result of calling get_meta_associations()
#' corresponding to the second decomposition, which is stored in
#' container$meta_associations
#' @param use_text logical If TRUE, then displays correlation coefficients in cells
#' (default=TRUE)
#'
#' @export
compare_decompositions <- function(tucker_res1,tucker_res2,decomp_names,meta_anno1,meta_anno2,use_text=TRUE) {
  ## first get donor scores comparison
  # ensure donors in same order
  tr1 <- tucker_res1[[1]]
  tr2 <- tucker_res2[[1]]

  # get donors present in both decompositions
  donors_use <- intersect(rownames(tr1),rownames(tr2))
  tr1 <- tr1[donors_use,]
  tr2 <- tr2[donors_use,]

  res_cor <- cor(tr1,tr2)
  rownames(res_cor) <- sapply(1:ncol(tr1),function(x){paste0('Factor',as.character(x))})
  colnames(res_cor) <- sapply(1:ncol(tr2),function(x){paste0('Factor',as.character(x))})

  res_orig <- res_cor

  # order max vals along the diagonal
  mx_dimension <- which(dim(res_cor)==max(dim(res_cor)))
  if (length(mx_dimension)>1) {
    mx_dimension <- 1
  }

  if (mx_dimension==1) {
    # order columns by max value in column
    col_maxes <- apply(abs(res_cor), 2, function(x) max(x, na.rm = TRUE))
    res_cor <- res_cor[,order(col_maxes,decreasing=TRUE)]
    # loop through columns, rearranging rows to make max on diagonal
    for (j in 1:ncol(res_cor)) {
      new_row_order <- order(abs(res_cor[j:nrow(res_cor),j]),decreasing=TRUE)
      res_tmp <- res_cor[j:nrow(res_cor),,drop=FALSE]
      res_tmp <- res_tmp[new_row_order,,drop=FALSE]
      res_cor[j:nrow(res_cor),] <- res_tmp
      rownames(res_cor)[j:nrow(res_cor)] <- rownames(res_tmp)
    }
  } else if (mx_dimension==2) {
    # order rows by max value in row
    row_maxes <- apply(abs(res_cor), 1, function(x) max(x, na.rm = TRUE))
    res_cor <- res_cor[order(row_maxes,decreasing=TRUE),]
    # loop through rows, rearranging columns to make max on diagonal
    for (j in 1:nrow(res_cor)) {
      new_col_order <- order(abs(res_cor[j,j:ncol(res_cor)]),decreasing=TRUE)
      res_tmp <- res_cor[,j:ncol(res_cor),drop=FALSE]
      res_tmp <- res_tmp[,new_col_order,drop=FALSE]
      res_cor[,j:ncol(res_cor)] <- res_tmp
      colnames(res_cor)[j:ncol(res_cor)] <- colnames(res_tmp)
    }
  }

  col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  new_row_order <- match(rownames(res_cor),rownames(res_orig))
  new_col_order <- match(colnames(res_cor),colnames(res_orig))
  col_fun_annot = colorRamp2(c(0, 1), c("white", "forest green"))
  la <- rowAnnotation(rsq=t(meta_anno1),col = list(rsq = col_fun_annot),
                      border=TRUE,annotation_name_side = "top")
  ba <- HeatmapAnnotation(rsq=t(meta_anno2),col = list(rsq = col_fun_annot),
                          border=TRUE,annotation_name_side = "left",show_legend=FALSE)
  dscores_hmap <- Heatmap(res_orig, name = "Pearson r",
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    column_names_gp = gpar(fontsize = 8),
                    row_names_gp = gpar(fontsize = 10),
                    col = col_fun,border=TRUE, show_column_names=TRUE,
                    show_row_names=TRUE,show_row_dend = FALSE,
                    show_column_dend = FALSE,
                    row_title = decomp_names[1],
                    row_title_gp = gpar(fontsize = 20),
                    column_title = decomp_names[2],
                    column_title_gp = gpar(fontsize = 20),
                    column_title_side = "bottom",
                    bottom_annotation=ba,
                    left_annotation=la,
                    row_names_side = "left",
                    row_order = new_row_order,
                    column_order = new_col_order,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      if (use_text) {
                        grid::grid.text(sprintf("%.2f", res_orig[i, j]), x, y, gp = gpar(fontsize = 10))
                      }
                    })

  ## now to get loadings comparison with same factor ordering
  # ensure genes in same order
  tr1 <- tucker_res1[[2]]
  tr2 <- tucker_res2[[2]]

  # get gene_ctype combos present in both decompositions
  gc_use <- intersect(colnames(tr1),colnames(tr2))
  tr1 <- t(tr1[,gc_use])
  tr2 <- t(tr2[,gc_use])

  res_cor <- cor(tr1,tr2)
  rownames(res_cor) <- sapply(1:ncol(tr1),function(x){paste0('Factor',as.character(x))})
  colnames(res_cor) <- sapply(1:ncol(tr2),function(x){paste0('Factor',as.character(x))})

  loadings_hmap <- Heatmap(res_cor, name = "Loadings Pearson r",
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    column_names_gp = gpar(fontsize = 8),
                    row_names_gp = gpar(fontsize = 10),
                    col = col_fun,border=TRUE, show_column_names=TRUE,
                    show_row_names=TRUE,show_row_dend = FALSE,
                    show_column_dend = FALSE,
                    row_title = decomp_names[1],
                    row_title_gp = gpar(fontsize = 20),
                    column_title = decomp_names[2],
                    column_title_gp = gpar(fontsize = 20),
                    column_title_side = "bottom",
                    bottom_annotation=ba,
                    left_annotation=la,
                    row_names_side = "left",
                    row_order = new_row_order,
                    column_order = new_col_order,
                    show_heatmap_legend = FALSE,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      if (use_text) {
                        grid::grid.text(sprintf("%.2f", res_cor[i, j]), x, y, gp = gpar(fontsize = 10))
                      }
                    })

  hmlist <- list(dscores_hmap,loadings_hmap)
  hmlist <- dscores_hmap + loadings_hmap

  draw(hmlist, padding = unit(c(2, 2, 10, 2), "mm")) # add space for titles
  decorate_heatmap_body("Pearson r", {
    grid::grid.text("Donor Scores Comparison", y = unit(1, "npc") + unit(2, "mm"), just = "bottom")
  })
  decorate_heatmap_body("Loadings Pearson r", {
    grid::grid.text("Loadings Comparison", y = unit(1, "npc") + unit(2, "mm"), just = "bottom")
  })

}


#' Plot dotplots for each factor to compare donor scores between meta groups
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param meta_var character The meta data variable to compare groups for
#'
#' @return a figure of comparison plots (one for each factor) placed in
#' container$plots$indv_meta_scores_associations
#' @export
plot_scores_by_meta <- function(container,meta_var) {
  dscores <- container[["tucker_results"]][[1]]

  meta <- container$scMinimal_full$metadata[,c('donors',meta_var)]
  meta <- unique(meta)
  rownames(meta) <- meta$donors
  meta$donors <- NULL
  meta_vals <- as.character(unique(meta[[meta_var]]))

  if (sum(is.na(meta_vals))>0) {
    meta_vals <- meta_vals[!is.na(meta_vals)]
  }

  # make all columns of meta to be factors
  for (i in 1:ncol(meta)) {
    meta[,i] <- as.factor(unlist(meta[,i]))
  }

  all_plots <- list()
  all_pvals <- data.frame(matrix(nrow=0,ncol=4))
  all_dat <- data.frame(matrix(nrow=0,ncol=3))

  for (j in 1:ncol(dscores)) {
    f <- dscores[,j]

    # limit rows of meta to those in dscores
    meta <- meta[names(f),,drop=FALSE]

    tmp <- as.data.frame(cbind(f,meta,rep(j,nrow(meta))))

    # get p-value by t-test for each group comparison (pairwise)
    g_compare <- utils::combn(meta_vals,2)
    for (i in 1:ncol(g_compare)) {
      g1 <- g_compare[1,i]
      g2 <- g_compare[2,i]
      t_res <- stats::t.test(tmp[tmp[[meta_var]]==g1,1], tmp[tmp[[meta_var]]==g2,1],
                      alternative = "two.sided",var.equal = FALSE)
      pval <- t_res$p.value
      all_pvals <- rbind(all_pvals,c(j,'dscore',g1,g2,pval))
      all_dat <- rbind(all_dat,tmp)
    }
  }

  colnames(all_pvals) <- c('myfactor','.y.','group1','group2','p.adj')
  all_pvals$myfactor <- as.factor(all_pvals$myfactor)
  all_pvals$group1 <- as.character(all_pvals$group1)
  all_pvals$group2 <- as.character(all_pvals$group2)


  colnames(all_dat) <- c('dscore','Status','myfactor')
  all_dat$myfactor <- as.factor(all_dat$myfactor)

  # apply fdr correction
  all_pvals$p.adj <- p.adjust(all_pvals$p.adj,method='fdr')

  # write p-vals as text
  all_pvals$p.adj <- sapply(all_pvals$p.adj,function(x) {
    paste0('p = ',as.character(round(x,digits=5)))
  })

  for (j in 1:ncol(dscores)) {
    tmp_dat <- all_dat[all_dat$myfactor==j,]
    tmp_pvals <- all_pvals[all_pvals$myfactor==j,,drop=FALSE]
    p <- ggplot(tmp_dat,aes(x=Status,y=dscore)) +
      geom_violin() +
      # geom_dotplot(binaxis = 'y', stackdir = 'center', method = 'histodot',
      #              dotsize = 2.5, binwidth = .005) +
      # geom_boxplot() +
      ggpubr::stat_pvalue_manual(
        tmp_pvals,
        y.position = max(tmp$f)+.2, step.increase = .1,
        label = "p.adj"
      ) +
      ylab('Score') +
      ggtitle(paste0('Factor ',as.character(j))) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_y_continuous(expand = expansion(mult = c(.15, 0.25)))

    all_plots[[j]] <- p
  }

  if (ncol(dscores) >= 5) {
    nc <- 5
    nr <- ceiling(ncol(dscores) / 5)
  } else {
    nc <- ncol(dscores)
    nr <- 1
  }
  p_final <- ggpubr::ggarrange(plotlist=all_plots, nrow = nr, ncol = nc)

  container$plots$indv_meta_scores_associations <- p_final

  return(container)
}



#' Create UMAP for donor distances
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param color_by_factor numeric Number of factor to color donors by. Can be a
#' vector of multiple factor numbers to make several plots. (default=NULL)
#' @param color_by_meta character Names of meta variables to color donors by.
#' Can be a vector of multiple names to make several plots. (default=NULL)
#' @param n_col numeric The number of columns to orde the figure into (default=1)
#'
#' @return the project container with the figure in container$plots$dscores_umap
#' @export
plot_donor_umap <- function(container, color_by_factor=NULL, color_by_meta=NULL, n_col=1) {

  dscores <- container$tucker_results[[1]]
  um <- as.data.frame(umap::umap(dscores)$layout)
  colnames(um) <- c('UMAP1','UMAP2')

  all_plots <- list()
  if (!is.null(color_by_factor)) {
    for (i in 1:length(color_by_factor)) {
      fact <- color_by_factor[i]
      score <- container$tucker_results[[1]][,fact]
      tmp <- as.data.frame(cbind(um,score))
      p <- ggplot(tmp,aes(x=UMAP1, y=UMAP2, color=score)) +
        geom_point() +
        scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                              high = "red", space = "Lab" ) +
        xlab('') +
        ylab('') +
        ggtitle(paste0("Factor ",as.character(fact))) +
        theme(plot.title = element_text(hjust = 0.5))

      all_plots[[i]] <- p
    }

    fig <- cowplot::plot_grid(plotlist=all_plots,ncol=n_col,scale = 0.95)

  } else if (!is.null(color_by_meta)) {
    for (i in 1:length(color_by_meta)) {
      mv <- color_by_meta[i]
      meta <- container$scMinimal_full$metadata[,c('donors',mv)]
      meta <- unique(meta)
      rownames(meta) <- meta$donors
      meta$donors <- NULL
      score <- meta[rownames(dscores),]
      tmp <- as.data.frame(cbind(um,score))
      p <- ggplot(tmp,aes(x=UMAP1, y=UMAP2, color=score)) +
        geom_point() +
        xlab('') +
        ylab('') +
        ggtitle(mv) +
        theme(plot.title = element_text(hjust = 0.5))

      all_plots[[i]] <- p
    }

    fig <- cowplot::plot_grid(plotlist=all_plots,ncol=n_col,scale = 0.95)

  } else {
    fig <- ggplot(tmp,aes(x=UMAP1, y=UMAP2)) +
      geom_point()
  }



  container$plots$dscores_umap <- fig

  return(container)

}




plot_dscore_enr <- function(container,factor_use,meta_var) {
  meta <- unique(container$scMinimal_full$metadata[,c('donors',meta_var)])
  rownames(meta) <- meta$donors

  meta_vals <- unlist(unique(as.character(meta[,meta_var])))
  mypaths <- list()
  for (i in 1:length(meta_vals)) {
    mypaths[[meta_vals[i]]] <- rownames(meta)[meta[,meta_var]==meta_vals[i]]
  }

  myranks <- container$tucker_results[[1]][,factor_use]

  fgseaRes <- fgsea(pathways = mypaths,
                    stats    = myranks,
                    minSize  = 0,
                    maxSize  = 5000)

  print(fgseaRes)

  plt_lst <- list()
  for (i in 1:length(meta_vals)) {
    plt <- plotEnrichment(mypaths[[meta_vals[i]]],
                              myranks) + labs(title=paste0(meta_vals[i],' - Factor ',as.character(factor_use)))
    plt <- plt +
      annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
               label=paste0('adj pval: ',
                            round(fgseaRes[fgseaRes$pathway==meta_vals[i],'padj'],digits=4)))

    plt_lst[[i]] <- plt
  }

  fig <- cowplot::plot_grid(plotlist=plt_lst,nrow=1)
  return(fig)
}



get_leading_edge_genes <- function(container,factor_select,gsets,num_genes_per=5) {
  factor_name <- paste0('Factor',as.character(factor_select))

  all_le <- list()
  for (gs in gsets) {
    # get ctypes where gs is significant
    ct_sig <- c()
    for (ct in container$experiment_params$ctypes_use) {
      gsea_res <- container[["gsea_res_full"]][[factor_name]][[ct]]
      padj <- gsea_res$padj[gsea_res$pathway==gs]
      if (padj < 0.05) {
        ct_sig <- c(ct_sig,ct)
      }
    }

    # loop through cell types where gs is significant and get gene counts
    g_counts <- list()
    for (ct in ct_sig) {
      gsea_res <- container[["gsea_res_full"]][[factor_name]][[ct]]
      le_genes <- gsea_res$leadingEdge[gsea_res$pathway==gs][[1]]
      for (g in le_genes) {
        if (g %in% names(g_counts)) {
          g_counts[[g]] <- g_counts[[g]] + 1
        } else {
          g_counts[[g]] <- 1
        }
      }
    }
    # unlist and order g_counts decreasing order
    g_counts <- unlist(g_counts)
    g_counts <- g_counts[order(g_counts,decreasing=TRUE)]

    # if (length(ct_sig)==1) {
    #   g_counts <- sample(g_counts)
    # }
    all_le[[gs]] <- g_counts
  }

  # ensure no genes are selected that are in multiple sets
  final_le <- list()
  for (i in 1:length(all_le)) {
    g_counts <- all_le[[i]]
    # g_counts <- sample(g_counts)
    # print(g_counts)
    track <- 0 #keeps track of number genes accepted as unique for the set
    ndx <- 1
    while (track < num_genes_per && ndx <= length(g_counts)) {
      mygene <- names(g_counts)[ndx]

      # test if gene in any other leading edge gene sets
      is_unique <- TRUE
      for (j in 1:length(all_le)) {
        if (j != i) {
          g_counts2 <- all_le[[j]]
          if (mygene %in% names(g_counts2)) {
            is_unique <- FALSE
            break
          }
        }
      }

      if (is_unique) {
        final_le[[mygene]] <- names(all_le)[i]
        track <- track + 1
      }
      ndx <- ndx + 1
    }
  }
  return(unlist(final_le))
}





