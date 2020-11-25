
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
#'
#' @return the project container with the plot in container$plots$donor_matrix
#' @export
plot_donor_matrix <- function(container, meta_vars=NULL, cluster_by_meta=NULL, show_donor_ids=FALSE) {

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
    myhmap <- Heatmap(as.matrix(donor_mat), name = "score",
                      cluster_columns = FALSE,show_column_dend = FALSE,
                      cluster_rows = TRUE, show_row_dend = FALSE,
                      column_names_gp = gpar(fontsize = 10),
                      col = col_fun, row_title = "Donors",
                      row_title_gp = gpar(fontsize = 14),
                      show_row_names = show_donor_ids,
                      border = TRUE)
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
    }

    myhmap <- Heatmap(as.matrix(donor_mat), name = "score",
                      cluster_columns = FALSE,
                      cluster_rows = FALSE,
                      column_names_gp = gpar(fontsize = 10),
                      col = col_fun, row_title = "Donors",
                      row_title_gp = gpar(fontsize = 14),
                      show_row_names = show_donor_ids,
                      border = TRUE)

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
#' the factor being plotted (default=NULL).
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
#' @param show_xlab logical If TRUE, displays the xlabel 'genes' (default=TRUE)
#' @param show_var_eplained logical If TRUE then shows an anottation with the explained variance
#' for each cell type (default=TRUE)
#' @param show_all_legends logical Set to TRUE in order to show legends (default=TRUE)
#'
#' @return container with the plot put in container$plots$single_lds_plot
#' @export
plot_loadings_annot <- function(container, factor_select, use_sig_only=FALSE, nonsig_to_zero=FALSE, annot='none',
                                pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE, 
                                gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, show_xlab=TRUE,
                                show_var_eplained=TRUE, show_all_legends=TRUE) {
  # check that Tucker has been run
  if (is.null(container$tucker_results)) {
    stop("Need to run run_tucker_ica() first.")
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
  
  hm_legends <- list()
  
  if (show_var_eplained) {
    ctype_var_exp <- get_explained_var(container, tmp_casted_num, factor_select)
    ctype_var_exp <- unlist(ctype_var_exp)
    ctype_var_exp <- ctype_var_exp[colnames(tmp_casted_num)]
    col_fun2 = circlize::colorRamp2(c(0, max(ctype_var_exp)), c("white", "black"))
    var_annot <- ComplexHeatmap::HeatmapAnnotation(var = ctype_var_exp,col=list(var=col_fun2),
                                                   show_annotation_name=FALSE, border=TRUE,
                                                   show_legend = show_all_legends)
    
    hm_legends[[2]] <- Legend(col_fun = col_fun2, title = "var exp",
                              grid_height = unit(1, "mm"), grid_width = unit(3, "mm"),
                              title_position = "leftcenter-rot")
    
  } else {
    var_annot <- NULL
  }
  
  # make colormap for hmap
  color_lim <- max(abs(tmp_casted_num))
  nintieth_per <- stats::quantile(as.matrix(abs(tmp_casted_num)), c(.95))
  if (color_lim > (1.5*nintieth_per)) {
    col_fun = colorRamp2(c(-nintieth_per, 0, nintieth_per), c("blue", "white", "red"))
  } else {
    col_fun = colorRamp2(c(-color_lim, 0, color_lim), c("blue", "white", "red"))
  }
  
  hm_legends[[1]] <- Legend(col_fun = col_fun, title = "loading",
                            grid_height = unit(1, "mm"), grid_width = unit(3, "mm"),
                            title_position = "leftcenter-rot")
  
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
                     show_heatmap_legend = show_all_legends)
  
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
              column_names_gp = gpar(fontsize = 20), border = TRUE)
  }

  if (!is.null(sim_de_donor_group)) {
    genes_df <- tmp_casted_num
    genes <- rownames(tmp_casted_num)
    for (j in 1:ncol(genes_df)) {
      for (g in genes) {
        groups <- strsplit(g,split = '_')
        cur_group <- paste0('Group',(ncol(genes_df)*(sim_de_donor_group-1))+j)
        if (cur_group %in% groups[[1]]) {
          genes_df[g,j] <- 1
        } else {
          genes_df[g,j] <- 0
        }
      }
    }


    # add gene significance heatmap to total hmap
    col_fun <- structure(c("white", "cyan"), names = c("0", "1"))
    hm_list <- hm_list +
      Heatmap(genes_df,
              name = "True DE Genes", cluster_columns = FALSE,
              col = col_fun, show_row_names = FALSE,
              show_heatmap_legend = TRUE, show_column_dend = FALSE,
              column_names_gp = gpar(fontsize = 20), border = TRUE)
  }

  # save plot in the container
  container$plots$single_lds_plot <- hm_list
  
  pd = packLegend(hm_legends[[1]], hm_legends[[2]], direction = "vertical")
  container$plots$ldng_legends <- pd
  
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
    strsplit(x,split = '.', fixed = TRUE)[[1]][[3]]
  })
  padj_use <- padj[which(padj_factors == as.character(factor_select))]

  # get out pvals for gene_ctype combos of specific ctypes
  padj_all_ctypes <- list()
  for (ct in ctypes) {
    padj_ct <- sapply(names(padj_use),function(x) {
      strsplit(x,split = '.', fixed = TRUE)[[1]][[2]]
    })
    padj_ct <- padj_use[which(padj_ct == ct)]

    names(padj_ct) <- sapply(names(padj_ct),function(x) {
      strsplit(x,split = '.',fixed = TRUE)[[1]][[1]]
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

#' Get explained variance for each cell type for one factor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param tmp_casted_num matrix The gene by cell type loadings matrix
#' @param factor_use numeric The factor to investigate
#'
#' @return explained variance for each cell type in a list
get_explained_var <- function(container, tmp_casted_num, factor_use) {
  rnk <- container$experiment_params$ranks
  tnsr <- container$tensor_data[[4]]
  donor_scores <- container$tucker_results[[1]]
  ctypes <- container$experiment_params$ctypes_use
  
  ctype_errors <- list()
  for (i in 1:length(ctypes)) {
    ct <- ctypes[i]
    
    # expression from data tensor (not loadings tensor)
    t_slice <- tnsr[,,i]
    
    # compute reconstruction of cell type gene expression
    ldngs <- tmp_casted_num[,ct,drop=FALSE]
    recon <-  as.matrix(donor_scores[,factor_use]) %*% as.matrix(t(ldngs))
    
    diff_mat <- t_slice
    colnames(diff_mat) <- container$tensor_data[[2]]
    rownames(diff_mat) <- container$tensor_data[[1]]
    row_ndx <- match(rownames(recon),rownames(diff_mat))
    col_ndx <- match(colnames(recon),colnames(diff_mat))
    diff_mat[row_ndx,col_ndx] <- diff_mat[row_ndx,col_ndx] - recon
      
    exp_var_rel <- 1-((norm(diff_mat,"F")**2)/(norm(t_slice,"F")**2))
    exp_var <- exp_var_rel * sum(apply(t_slice, 2, stats::var))
    
    ctype_errors[[ct]] <- exp_var
  }
  return(ctype_errors)
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
#' @param show_var_eplained logical If TRUE then shows an anottation with the explained variance
#' for each cell type (default=TRUE)
#'
#' @return the project container with the list of plots placed in container$plots$all_lds_plots
#' @export
get_all_lds_factor_plots <- function(container, use_sig_only=FALSE, nonsig_to_zero=FALSE, annot='none',
                                     pathways_list=NULL, sig_thresh=0.05, display_genes=FALSE,
                                     gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL,
                                     show_var_eplained=TRUE) {

  # save any plot previously in the single lds plot slot because will be overwrittern
  prev_lds_plot <- container$plots$single_lds_plot

  num_fact <- nrow(container$tucker_results[[2]])
  hm_list <- list()
  lgnd_list <- list()
  for (i in 1:num_fact) {
    if (i==1 && !gene_callouts) {
      show_xlab <- TRUE
    } else {
      show_xlab <- FALSE
    }
    container <- plot_loadings_annot(container, factor_select=i, 
                                     use_sig_only=use_sig_only,
                                     nonsig_to_zero=nonsig_to_zero,
                                     annot=annot, pathways=pathways_list[[i]],
                                     sig_thresh=sig_thresh,
                                     display_genes=display_genes,
                                     gene_callouts=gene_callouts,
                                     callout_n_gene_per_ctype=callout_n_gene_per_ctype,
                                     callout_ctypes=callout_ctypes[[i]],
                                     show_xlab=show_xlab,
                                     show_var_eplained=show_var_eplained,
                                     show_all_legends=FALSE)
    
    hm_list[[i]] <- container$plots$single_lds_plot
    lgnd_list[[i]] <- container$plots$ldng_legends
  }

  # rewrite the single plot back into its original slot
  container$plots$single_lds_plot <- prev_lds_plot

  # save list of all plots
  container$plots$all_lds_plots <- hm_list
  container$plots$all_legends <- lgnd_list
  
  return(container)
}

#' Creates a figure of all loadings plots side by side
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param n_rows numeric The number of rows to fit the plots onto
#' @export
render_all_lds_plots <- function(container,n_rows) {
  
  hm_list <- container$plots$all_lds_plots
  hm_legends <- container$plots$all_legends
  
  grid::grid.newpage()
  
  grid_n_col <- ceiling(length(hm_list)/n_rows)
  for (i in 1:length(hm_list)) {
    col_ndx <- i %% grid_n_col
    if (col_ndx == 0) {
      col_ndx <- grid_n_col
    }
    row_ndx <- n_rows - ceiling(i / grid_n_col) + 1
    
    if (length(hm_list)%%grid_n_col!=0 && row_ndx==1) {
      x_buffer <- ((1/grid_n_col)-.02) * (grid_n_col - (length(hm_list)%%grid_n_col))/2
      x_pos <- (col_ndx-1)*(1/grid_n_col) + x_buffer
    } else {
      x_pos <- (col_ndx-1)*(1/grid_n_col)
    }
    y_pos <- row_ndx/n_rows
    
    grid::pushViewport(grid::viewport(x = x_pos, y = y_pos, 
                                      width = (1/grid_n_col)-.02, height = 1/n_rows,
                                      just = c("left","top")))
    draw(hm_list[[i]],annotation_legend_list = hm_legends[[i]], 
         legend_grouping = "original",
         newpage=FALSE)
    grid::popViewport()
    
  }
  
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
#'
#' @return the project container with the plot in the slot
#' container$plots$donor_sig_genes$Factor#
#' @export
plot_donor_sig_genes <- function(container, factor_select, top_n_per_ctype, ctypes_use=NULL) {
  ## add catch in case they havent run jackstraw yet...
  
  # temporarily remove variance scaling
  orig_scale_decision <- container$experiment_params$scale_var
  container <- set_experiment_params(container, scale_var = FALSE)
  
  # form the tensor for specified cell types
  container <- collapse_by_donors(container, shuffle=FALSE)
  container <- form_tensor(container)
  
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
  ct_show <- as.factor(ct_show)
  
  ct_annot <- ComplexHeatmap::rowAnnotation(cell_types=anno_simple(ct_show),
                                            show_annotation_name=FALSE)

  # create the hmap
  col_fun = colorRamp2(c(min(donor_unfold_sub), 0, max(donor_unfold_sub)), c("blue", "white", "red"))

  myhmap <- Heatmap(donor_unfold_sub, name = "expr",
                    cluster_columns = FALSE,
                    cluster_rows = TRUE,
                    column_names_gp = gpar(fontsize = 8),
                    row_names_gp = gpar(fontsize = 10),
                    col = col_fun, top_annotation=ha, row_split = ct_show,
                    row_labels=rn_show,border=TRUE, show_column_names=TRUE,
                    left_annotation=ct_annot,show_row_dend = FALSE,
                    column_title = paste0('Factor ', factor_select,' Top Genes'),
                    column_title_gp = gpar(fontsize = 20, fontface = "bold"))

  # reset scale variance decision
  container <- set_experiment_params(container, scale_var = orig_scale_decision)

  container$plots$donor_sig_genes[[paste0('Factor',as.character(factor_select))]] <- myhmap
  return(container)
}









