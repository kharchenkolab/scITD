

#' Run fgsea on expression levels as summed after multiplying by donor scores
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The factor of interest
#' @param sig_thresh numeric The pvalue threshold for selecting gene sets
#' @param ctype character The cell type of interest
#' @param db_use character The database of gene sets to use. Database
#' options include "GO", "Reactome", "KEGG", and "BioCarta". More than
#' one database can be used. (default="GO")
#' @param collapse_paths logical Set to TRUE to collapse significant pathways with
#' considerable overlap (default=TRUE)
#'
#' @return data.frame of the fgsea results (including non-significant results)
#' @export
run_fgsea <- function(container, factor_select, sig_thresh, ctype, db_use="GO",
                      collapse_paths=TRUE) {
  donor_scores <- container$tucker_results[[1]]

  # select mean exp data for one cell type
  tnsr_slice <- container$scMinimal_ctype[[ctype]]$pseudobulk
  tnsr_slice <- scale(tnsr_slice, center=TRUE) # scaling to unit variance!
  # tnsr_slice <- tnsr_slice[rownames(donor_scores),]

  # ## testing use of gene significance...
  # sig_vectors <- get_significance_vectors(container,
  #                                         factor_select, container$experiment_params$ctypes_use)
  # sig_ctype <- sig_vectors[[ctype]]
  # min_nonzero <- min(sig_ctype[sig_ctype!=0])
  # sig_ctype[sig_ctype==0] <- min_nonzero
  # sig_ctype <- -log10(sig_ctype)
  # ##

  # get transformed expression for each gene by summing d_score * scaled exp
  exp_vals <- sapply(1:ncol(tnsr_slice), function(j) {
    exp_transform <- tnsr_slice[,j] * donor_scores[rownames(tnsr_slice),factor_select]
    de_val <- sum(exp_transform)

    # exp_transform <- tnsr_slice[,j] * donor_scores[rownames(tnsr_slice),factor_select]
    # de_val <- abs(sum(exp_transform))

    # exp_transform <- tnsr_slice[,j] * (donor_scores[rownames(tnsr_slice),factor_select]**2)
    # mysigns <- donor_scores[rownames(tnsr_slice),factor_select] < 0
    # exp_transform[mysigns] <- exp_transform[mysigns] * -1
    # de_val <- sum(exp_transform)

    # de_val <- cor(tnsr_slice[,j],donor_scores[rownames(tnsr_slice),factor_select])

    ## testing get a score based on significance from jackstraw that is signed appropriately
    # to avoid excessive # of ties, I'll use significance as a weighting factor
    # exp_transform <- tnsr_slice[,j] * donor_scores[rownames(tnsr_slice),factor_select]
    # de_val <- sum(exp_transform)
    # de_val <- de_val * sig_ctype[j]

    # # Trying to use lm slope as statistic
    # tmp <- as.data.frame(cbind(tnsr_slice[,j],donor_scores[rownames(tnsr_slice),factor_select]))
    # colnames(tmp) <- c('expr','dsc')
    # lmres <- lm(expr~dsc,data=tmp)
    # tmp2 <- summary(lmres)
    # de_val <- tmp2$coefficients['dsc','Estimate']

    return(de_val)

  })

  names(exp_vals) <- convert_gn(container,colnames(tnsr_slice))

  # remove duplicate genes
  ndx_remove <- duplicated(names(exp_vals)) | duplicated(names(exp_vals), fromLast = TRUE)
  exp_vals <- exp_vals[!ndx_remove]

  m_df <- data.frame()
  for (db in db_use) {
    if (db == "GO") {
      # select the GO Biological Processes group of gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C5", subcategory = "BP"))
    } else if (db == "Reactome") {
      # select the Reactome gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C2", subcategory = "CP:REACTOME"))
    } else if (db == "KEGG") {
      # select the KEGG gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C2", subcategory = "CP:KEGG"))
    } else if (db == "BioCarta") {
      # select the BioCarts gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C2", subcategory = "CP:BIOCARTA"))
    } else if (db == "Hallmark") {
      # select the BioCarts gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "H"))
    }
  }

  my_pathways <- split(m_df$gene_symbol, f = m_df$gs_name)
  # my_pathways <- split(m_df$gene_symbol, f = m_df$gs_exact_source)

  fgsea_res <- fgsea::fgsea(pathways = my_pathways,
                            stats = exp_vals,
                            minSize=15,
                            maxSize=500,
                            eps=0,
                            gseaParam=1,
                            nproc=container$experiment_params$ncores)
  # fgsea_res <- fgsea::fgsea(pathways = my_pathways,
  #                           stats = exp_vals,
  #                           minSize=15,
  #                           maxSize=500,
  #                           eps=0,
  #                           gseaParam=1,
  #                           scoreType = "pos",
  #                           nproc=container$experiment_params$ncores)
  # fgsea_res <- fgsea::fgsea(pathways = my_pathways,
  #                           stats = exp_vals,
  #                           minSize=15,
  #                           maxSize=400,
  #                           gseaParam=1,
  #                           nperm=10000,
  #                           nproc=container$experiment_params$ncores)

  fgsea_res <- fgsea_res[order(fgsea_res$padj, decreasing=FALSE),]

  fgsea_lim <- fgsea_res[fgsea_res$padj < sig_thresh,]
  # print(ctype)
  # print(fgsea_lim[,1:5])

  if (collapse_paths) {
    # collapse significant pathways
    fgsea_lim <- fgsea_res[fgsea_res$padj < sig_thresh,]
    # print(fgsea_lim)
    cpaths <- fgsea::collapsePathways(fgseaRes=fgsea_lim,pathways=my_pathways,stats=exp_vals,
                                      pval.threshold=.05,gseaParam=1)
    main_paths <- cpaths$mainPathways
  } else {
    main_paths <- NULL
  }


  return(list(fgsea_res,main_paths))
}



#' Compute enriched gene sets in a cell type for a factor using fgsea
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The factor of interest
#' @param ctype character The cell type of interest
#' @param db_use character The database of gene sets to use. Database
#' options include "GO", "Reactome", "KEGG", and "BioCarta". More than
#' one database can be used. (default="GO")
#'
#' @return data.frame of the fgsea results limited to genes passing the
#' significance threshold.
#' @export
run_fgsea_old <- function(container, factor_select, ctype,
                      db_use="GO", sig_thresh, collapse_paths) {

  set.seed(container$experiment_params$rand_seed)

  ldngs <- container$tucker_results[[2]]

  # prep the loadings matrix
  genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
  ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})

  sr_col <- ldngs[factor_select,]

  tmp_casted_num <- reshape_loadings(sr_col, genes, ctypes)

  ctype_lds <- tmp_casted_num[,ctype]

  names(ctype_lds) <- convert_gn(container, rownames(tmp_casted_num))

  # check for duplicate genes
  duplicates <- names(ctype_lds)[duplicated(names(ctype_lds))]
  genes_remove <- names(ctype_lds) %in% duplicates
  ctype_lds <- ctype_lds[!genes_remove]

  ## testing use of gene significance...
  sig_vectors <- get_significance_vectors(container,
                                          factor_select, container$experiment_params$ctypes_use)
  sig_ctype <- sig_vectors[[ctype]]
  min_nonzero <- min(sig_ctype[sig_ctype!=0])
  sig_ctype[sig_ctype==0] <- min_nonzero
  sig_ctype <- sqrt(-log10(sig_ctype))

  ctype_lds <- ctype_lds * sig_ctype[names(ctype_lds)]
  ##

  # # trying throwing out 0 value genes...
  # ctype_lds <- ctype_lds[ctype_lds!=0]

  m_df <- data.frame()
  for (db in db_use) {
    if (db == "GO") {
      # select the GO Biological Processes group of gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C5", subcategory = "BP"))
    } else if (db == "Reactome") {
      # select the Reactome gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C2", subcategory = "CP:REACTOME"))
    } else if (db == "KEGG") {
      # select the KEGG gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C2", subcategory = "CP:KEGG"))
    } else if (db == "BioCarta") {
      # select the BioCarts gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C2", subcategory = "CP:BIOCARTA"))
    }
  }

  my_pathways <- split(m_df$gene_symbol, f = m_df$gs_name)
  fgsea_res <- fgsea::fgsea(pathways = my_pathways,
                            stats = ctype_lds,
                            minSize=15,
                            maxSize=500,
                            eps=0,
                            gseaParam=1,
                            nPermSimple = 10000)

  fgsea_res <- fgsea_res[order(fgsea_res$padj, decreasing=FALSE),]

  return(list(fgsea_res,NULL))
}



#' Compute enriched gene sets among significant genes in a cell type for
#' a factor using hypergeometric test
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The factor of interest
#' @param ctype character The cell type of interest
#' @param up_down character Either "up" to compute enrichment among the significant
#' positive loading genes or "down" to compute enrichment among the significant
#' negative loading genes.
#' @param thresh numeric Pvalue significance threshold. Used for significant
#' gene selection (default=0.05)
#' @param db_use character The database of gene sets to use. Database
#' options include "GO", "Reactome", "KEGG", and "BioCarta". More than
#' one database can be used. (default="GO")
#'
#' @return pvalues for all tested gene sets
#' @export
run_hypergeometric_gsea <- function(container, factor_select, ctype, up_down,
                                     thresh=0.05, db_use="GO") {

  # make sure Tucker has been run
  if (is.null(container$tucker_results)) {
    stop('Run run_tucker_ica() first.')
  }

  ldngs <- container$tucker_results[[2]]

  # prep the loadings matrix
  genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
  ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})

  sr_col <- ldngs[factor_select,]

  tmp_casted_num <- reshape_loadings(sr_col, genes, ctypes)

  sig_vectors <- get_significance_vectors(container,
                                          factor_select, colnames(tmp_casted_num))
  # convert list to df
  sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

  # limit to just the genes in tmp_casted_num
  sig_df <- sig_df[rownames(tmp_casted_num),colnames(tmp_casted_num)]

  sig_df <- sig_df[,ctype,drop=FALSE]

  sig_genes <- rownames(sig_df)[sig_df<thresh]
  sig_genes <- convert_gn(container, sig_genes)
  sig_genes <- unique(sig_genes)

  if (up_down == "up") {
    sig_genes_updown <- rownames(tmp_casted_num)[tmp_casted_num[,ctype]>0]
    sig_genes_updown <- convert_gn(container, sig_genes_updown)
    sig_genes <- sig_genes[sig_genes %in% sig_genes_updown]
  } else if (up_down == "down") {
    sig_genes_updown <- rownames(tmp_casted_num)[tmp_casted_num[,ctype]<0]
    sig_genes_updown <- convert_gn(container, sig_genes_updown)
    sig_genes <- sig_genes[sig_genes %in% sig_genes_updown]
  }

  all_genes <- convert_gn(container, rownames(tmp_casted_num))
  total_num_genes <- length(all_genes)

  m_df <- data.frame()
  for (db in db_use) {
    if (db == "GO") {
      # select the GO Biological Processes group of gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                 category = "C5", subcategory = "BP"))
    } else if (db == "Reactome") {
      # select the Reactome gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                 category = "C2", subcategory = "CP:REACTOME"))
    } else if (db == "KEGG") {
      # select the KEGG gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                 category = "C2", subcategory = "CP:KEGG"))
    } else if (db == "BioCarta") {
      # select the BioCarts gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                 category = "C2", subcategory = "CP:BIOCARTA"))
    }
  }

  my_pathways = split(m_df$gene_symbol, f = m_df$gs_name)
  # my_pathways <- split(m_df$gene_symbol, f = m_df$gs_exact_source)

  pvals <- c()
  for (i in 1:length(my_pathways)) {
    pth <- my_pathways[[i]]
    pth_name <- names(my_pathways)[i]

    # A: total num genes in pathway in tmp_casted_num
    pth_in_df <- pth[which(pth %in% all_genes)]
    num_pth_in_df <- length(pth_in_df)

    # if set is too small continue
    if (num_pth_in_df < 15) {
      next
    }

    # B: number of genes from A in sig_genes
    num_in_sig <- sum(pth_in_df %in% sig_genes)

    # compute pvalue
    pval <- stats::phyper(num_in_sig-1, num_pth_in_df, total_num_genes - num_pth_in_df,
           length(sig_genes), lower.tail = FALSE) # I double checked this is right
    pvals[pth_name] <- pval
  }
  padj <- p.adjust(pvals,method='fdr')
  return(padj)
}

#' Run gsea separately for all cell types of one specified factor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The factor of interest
#' @param method character The method of gsea to use. Can either be "fgsea",
#' "fgsea_special or "hypergeometric". (default="fgsea")
#' @param thresh numeric Pvalue significance threshold to use. Will include gene sets in
#' resulting heatmap if pvalue is below this threshold for at least one cell type. (default=0.05)
#' @param db_use character The database of gene sets to use. Database
#' options include "GO", "Reactome", "KEGG", and "BioCarta". More than
#' one database can be used. (default="GO")
#' @param collapse_paths logical Set to TRUE to collapse significant pathways with
#' considerable overlap. Only used if method="fgsea". (default=TRUE)
#' @param reset_other_factor_plots logical Set to TRUE to set all other gsea plots to NULL (default=FALSE)
#' @param draw_plot logical Set to TRUE to show the plot. Plot is stored regardless. (default=TRUE)
#'
#' @return a heatmap plot of the gsea results in the slot
#' container$plots$gsea$FactorX
#' @export
run_gsea_one_factor <- function(container, factor_select, method="fgsea", thresh=0.05,
                                db_use="GO", collapse_paths=TRUE,
                                reset_other_factor_plots=FALSE, draw_plot=TRUE) {

  if (reset_other_factor_plots) {
    container$plots$gsea <- NULL
    container$gsea_results <- NULL
  }

  up_sets_all <- list()
  down_sets_all <- list()
  set_union <- c()
  ctypes_use <- container$experiment_params$ctypes_use
  for (ct in ctypes_use) {
    print(ct)
    if (method == 'fgsea') {
      fgsea_res <- run_fgsea(container, factor_select=factor_select,
                             sig_thresh=thresh, ctype=ct, db_use=db_use,
                             collapse_paths=collapse_paths)

      main_paths <- fgsea_res[[2]]
      fgsea_res <- fgsea_res[[1]]

      # remove results where NES is na
      fgsea_res <- fgsea_res[!is.na(fgsea_res$NES),]

      # keep separate track of positive/negative enriched sets
      up_sets_names <- fgsea_res$pathway[fgsea_res$NES > 0]
      up_sets <- fgsea_res$padj[fgsea_res$NES > 0]
      names(up_sets) <- up_sets_names
      down_sets_names <- fgsea_res$pathway[fgsea_res$NES < 0]
      down_sets <- fgsea_res$padj[fgsea_res$NES < 0]
      names(down_sets) <- down_sets_names
      up_sets_all[[ct]] <- up_sets
      down_sets_all[[ct]] <- down_sets
      set_union <- unique(c(set_union,main_paths))
      # store fgsea_res for access later on
      container$gsea_res_full[[paste0('Factor',factor_select)]][[ct]] <- fgsea_res
    } else if (method == 'hypergeometric') {
      gsea_res_up <- run_hypergeometric_gsea(container, factor_select=factor_select, ctype=ct,
                                             up_down='up', thresh=.1, db_use=db_use)
      gsea_res_down <- run_hypergeometric_gsea(container, factor_select=factor_select, ctype=ct,
                                               up_down='down', thresh=.1, db_use=db_use)

      up_sets_all[[ct]] <- gsea_res_up
      down_sets_all[[ct]] <- gsea_res_down
      set_union <- NULL
    }
  }

  # add results to container
  container$gsea_results[[as.character(factor_select)]] <- list('up'=up_sets_all,
                                                'down'=down_sets_all)

  # plot results
  myplot <- plot_gsea_hmap(container,factor_select,thresh,set_union=set_union)
  container$plots$gsea[[as.character(factor_select)]] <- myplot

  if (draw_plot) {
    draw(myplot,newpage=FALSE)
  }

  return(container)
}


#' Extract gene sets that are enriched selectively in one or more cell types for a given factor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The factor to investigate
#' @param these_ctypes_only character A vector of cell types for which to get gene sets that are
#' enriched in all of these and not in any other cell types
#' @param up_down character Set to "up" to get the gene sets for the positive loading genes. Set
#' to "down" to get the gene sets for the negative loadings genes.
#' @param thresh numeric Pvalue significance threshold for selecting enriched sets (default=0.05)
#'
#' @return a vector of pathways selectively enriched in the listed cell types
#' @export
get_intersecting_pathways <- function(container, factor_select, these_ctypes_only, up_down, thresh=0.05) {
  factor_select <- paste0('Factor',factor_select)
  ct_1_paths <- container$gsea_results[[factor_select]][[up_down]][[these_ctypes_only[1]]]
  intersect_pathways <- names(ct_1_paths)[ct_1_paths<thresh]
  if (length(these_ctypes_only) > 1) {
    for (i in 2:length(these_ctypes_only)) {
      cur_paths <- container$gsea_results[[factor_select]][[up_down]][[these_ctypes_only[i]]]
      sig_cur_paths <- names(cur_paths)[cur_paths<thresh]
      intersect_pathways <- intersect(intersect_pathways,
                                      sig_cur_paths)
    }
  }

  # now remove gene sets in any other pathway
  all_ctypes <- container$experiment_params$ctypes_use
  other_cts <- all_ctypes[!(all_ctypes %in% these_ctypes_only)]
  exclude_pathways <- c()
  for (i in 1:length(other_cts)) {
    exclude_pathways <- union(exclude_pathways,
                                    container$gsea_results[[factor_select]][[up_down]][[other_cts[i]]])
  }
  intersect_pathways_final <- intersect_pathways[!(intersect_pathways %in% exclude_pathways)]
  return(intersect_pathways_final)
}


#' Plot enriched gene sets from all cell types in a heatmap
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The factor to plot
#' @param thresh numeric Pvalue threshold to use for including gene sets in the heatmap
#' @param set_union character A vector of the main pathways from collapsePathways, and it is the
#' union of such pathways from all cell types of the factor (default=NULL)
#'
#' @return the heatmap plot
#' @export
plot_gsea_hmap <- function(container,factor_select,thresh,set_union=NULL) {
  factor_name <- paste0('Factor',as.character(factor_select))
  gsea_res <- container$gsea_results[[as.character(factor_select)]]

  df_list <- list()
  for (k in 1:length(gsea_res)) {
    up_down_sets <- gsea_res[[k]]

    # get unique gene sets
    all_sets <- c()
    for (i in 1:length(up_down_sets)) {
      all_sets <- c(all_sets, names(up_down_sets[[i]]))
    }
    all_sets <- unique(all_sets)
    all_sets <- all_sets[!is.na(all_sets)]

    res <- data.frame(matrix(1,ncol=length(up_down_sets),nrow = length(all_sets)))
    colnames(res) <- names(up_down_sets)
    rownames(res) <- all_sets

    for (i in 1:length(up_down_sets)) {
      ctype_res <- up_down_sets[[i]]
      ctype <- names(up_down_sets)[i]
      for (j in 1:length(ctype_res)) {
        res[names(ctype_res)[j],ctype] <- ctype_res[j]
      }
    }

    res_plot <- res[rowSums(res<thresh)>0,]

    if (nrow(res_plot) == 0) {
      next
    }

    if (!is.null(set_union)) {
      res_plot <- res_plot[rownames(res_plot)%in%set_union,]
    }

    df_list[[names(gsea_res)[k]]] <- res_plot
  }

  if (length(df_list)==0) {
    return(NULL)
  }

  # trying to build two separate hmaps and then concatenate them vertically
  for (k in 1:length(df_list)) {
    res_total <- df_list[[k]]

    # # parse gene set names at first underscore
    # rownames(res_total) <- sapply(rownames(res_total),function(x) {
    #   regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
    # })
    #
    # # cutoff gene set names
    # tmp_names <- sapply(rownames(res_total),function(x) {
    #   if (nchar(x) > 42) {
    #     return(paste0(substr(x,1,40),"..."))
    #   } else {
    #     return(x)
    #   }
    # })
    tmp_names <- rownames(res_total)

    if (k==1) {
      col_fun <- colorRamp2(c(.05, 0), c("white", "red"))
    } else {
      col_fun <- colorRamp2(c(.05, 0), c("white", "blue"))
    }

    myhmap <- Heatmap(as.matrix(res_total), name = paste0(names(df_list)[k],' pval'),
                      row_labels = tmp_names, row_title = names(df_list)[k],
                      show_row_dend = FALSE, show_column_dend = FALSE,
                      column_names_gp = gpar(fontsize = 10),
                      col = col_fun,
                      show_row_names = TRUE,
                      row_title_gp = gpar(fontsize = 14),
                      row_names_gp = gpar(fontsize = 6),
                      column_title = factor_name,
                      column_title_side = "top",
                      column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                      border=TRUE)
    if (k == 1) {
      hmap_list <- myhmap
    } else {
      hmap_list <- hmap_list %v% myhmap
    }
  }

  return(hmap_list)
}







##### making versions of these functions that work with GO names instead of IDs
#####

#' Plot enriched gene sets from all cell types in a heatmap
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The factor to plot
#' @param thresh numeric Pvalue threshold to use for including gene sets in the heatmap
#' @param set_union character A vector of the main pathways from collapsePathways, and it is the
#' union of such pathways from all cell types of the factor (default=NULL)
#'
#' @return the heatmap plot
#' @export
plot_gsea_hmap_w_similarity <- function(container,factor_select,direc,thresh) {
  factor_name <- paste0('Factor',as.character(factor_select))
  up_down_sets <- container$gsea_results[[as.character(factor_select)]][[direc]]

  # get unique gene sets
  all_sets <- c()
  for (i in 1:length(up_down_sets)) {
    all_sets <- c(all_sets, names(up_down_sets[[i]]))
  }
  all_sets <- unique(all_sets)
  all_sets <- all_sets[!is.na(all_sets)]

  res <- data.frame(matrix(1,ncol=length(up_down_sets),nrow = length(all_sets)))
  colnames(res) <- names(up_down_sets)
  rownames(res) <- all_sets

  # populate res with pvalues
  for (i in 1:length(up_down_sets)) {
    ctype_res <- up_down_sets[[i]]
    ctype <- names(up_down_sets)[i]
    for (j in 1:length(ctype_res)) {
      res[names(ctype_res)[j],ctype] <- ctype_res[j]
    }
  }

  # select only rows with at least one significant enrichment
  res_plot <- res[rowSums(res<thresh)>0,]

  if (nrow(res_plot) == 0) {
    return('no significant gene sets')
  }

  tmp_names <- rownames(res_plot)

  # convert from GO name to GO id for simplifyEnrichment
  gs <- msigdbr::msigdbr(species = "Homo sapiens",category = "C5", subcategory = "BP")
  gs <- gs[,c('gs_exact_source','gs_name')]
  gs <- as.data.frame(unique(gs))
  rownames(gs) <- gs$gs_name
  tmp_names <- gs[tmp_names,'gs_exact_source']

  mat <- simplifyEnrichment::GO_similarity(tmp_names,ont='BP')
  cl <- simplifyEnrichment::binary_cut(mat)
  sim_hmap_res <- ht_clusters(mat, cl, word_cloud_grob_param = list(max_width = 80))
  sim_hmap <- sim_hmap_res[[1]]
  ordering <- sim_hmap_res[[2]]

  col_fun <- colorRamp2(c(thresh, 0), c("white", "blue"))

  # show only some row names, taking reordering into account
  ndx_lab <- seq(from=1,to=nrow(res_plot),by=8)
  ndx_no_lab <- c(1:nrow(res_plot))[!(1:nrow(res_plot) %in% ndx_lab)]
  to_null <- ordering[ndx_no_lab]
  rlabs <- rownames(res_plot)
  rlabs[to_null] <- ''

  myhmap <- Heatmap(as.matrix(res_plot), name = paste0(direc,' pval'),
                    show_row_dend = FALSE, show_column_dend = FALSE,
                    column_names_gp = gpar(fontsize = 10),
                    col = col_fun,
                    row_order = ordering,
                    row_title_gp = gpar(fontsize = 14),
                    column_title = paste0('Factor ',factor_select," ",direc,' gene sets'),
                    column_title_side = "top",
                    column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                    border=TRUE,
                    width = unit(8, "cm"),
                    row_labels=rlabs,
                    row_names_side = "left",
                    row_names_gp = gpar(fontsize = 6))

  hm_list <- myhmap + sim_hmap

  draw(hm_list)

  decorate_heatmap_body("Similarity", {
    grid.rect(gp = gpar(fill = NA, col = "#404040"))
    cl = factor(cl, levels = unique(cl[ordering]))
    tbcl = table(cl)
    ncl = length(cl)
    x = cumsum(c(0, tbcl))/ncl
    grid.segments(x, 0, x, 1, default.units = "npc", gp = gpar(col = "#404040"))
    grid.segments(0, 1 - x, 1, 1 - x, default.units = "npc", gp = gpar(col = "#404040"))
  })

  decorate_heatmap_body(paste0(direc,' pval'), {
    grid.rect(gp = gpar(fill = NA, col = "#404040"))
    cl = factor(cl, levels = unique(cl[ordering]))
    tbcl = table(cl)
    ncl = length(cl)
    x = cumsum(c(0, tbcl))/ncl
    grid.segments(0, 1 - x, 1, 1 - x, default.units = "npc", gp = gpar(col = "#404040"))
  })

  # store cluster ordering and assignment if want to select later on
  container$gsea_last_info <- list(res_plot,cl,ordering)

  return(hm_list)
}

plot_gsea_sub <- function(container,factor_select,direc,thresh,clust_select) {
  res_plot <- container$gsea_last_info[[1]]
  cl <- container$gsea_last_info[[2]]
  ordering <- container$gsea_last_info[[3]]

  true_clust_order <- unique(cl[ordering])
  for (i in 1:length(cl)) {
    c_val <- cl[i]
    new_c_val <- which(true_clust_order == c_val)
    cl[i] <- new_c_val
  }

  # get GO gene set names present in the cluster of interest
  go_keep <- rownames(res_plot)[which(cl==clust_select)]

  # order res_plot by semantic similarity
  res_plot <- res_plot[ordering,]

  # keep only GO sets in cluster of interest
  res_plot_sub <- res_plot[rownames(res_plot) %in% go_keep,]

  col_fun <- colorRamp2(c(thresh, 0), c("white", "blue"))

  myhmap <- Heatmap(as.matrix(res_plot_sub), name = paste0(direc,' pval'),
                    show_row_dend = FALSE, show_column_dend = FALSE,
                    column_names_gp = gpar(fontsize = 10),
                    col = col_fun,
                    row_title_gp = gpar(fontsize = 12),
                    column_title = paste0('Factor ',factor_select," ",direc,' gene sets, cluster ', clust_select),
                    column_title_side = "top",
                    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                    border=TRUE,
                    width = unit(8, "cm"),
                    row_names_side = "left",
                    row_names_gp = gpar(fontsize = 6))

  return(myhmap)

}


# == title
# Visualize the similarity matrix and the clustering
#'@import grid
# == param
# -mat A similarity matrix.
# -cl Cluster labels inferred from the similarity matrix, e.g. from `cluster_terms` or `binary_cut`.
# -dend Used internally.
# -col A vector of colors that map from 0 to the 95^th percentile of the similarity values.
# -draw_word_cloud Whether to draw the word clouds.
# -term The full name or the description of the corresponding GO IDs.
# -min_term Minimal number of functional terms in a cluster. All the clusters
#     with size less than ``min_term`` are all merged into one separated cluster in the heatmap.
# -order_by_size Whether to reorder clusters by their sizes. The cluster
#      that is merged from small clusters (size < ``min_term``) is always put to the bottom of the heatmap.
# -cluster_slices Whether to cluster slices.
# -exclude_words Words that are excluded in the word cloud.
# -max_words Maximal number of words visualized in the word cloud.
# -word_cloud_grob_param A list of graphic parameters passed to `word_cloud_grob`.
# -fontsize_range The range of the font size. The value should be a numeric vector with length two.
#       The minimal font size is mapped to word frequency value of 1 and the maximal font size is mapped
#       to the maximal word frequency. The font size interlopation is linear.
# -column_title Column title for the heatmap.
# -ht_list A list of additional heatmaps added to the left of the similarity heatmap.
# -use_raster Whether to write the heatmap as a raster image.
# -... Other arguments passed to `ComplexHeatmap::draw,HeatmapList-method`.
#
# == value
# A `ComplexHeatmap::HeatmapList-class` object.
#
# == example
# mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds",
#     package = "simplifyEnrichment"))
# cl = binary_cut(mat)
# ht_clusters(mat, cl, word_cloud_grob_param = list(max_width = 80))
# ht_clusters(mat, cl, word_cloud_grob_param = list(max_width = 80),
#     order_by_size = TRUE)
ht_clusters = function(mat, cl, dend = NULL, col = c("white", "red"),
                       draw_word_cloud = simplifyEnrichment:::is_GO_id(rownames(mat)[1]) || !is.null(term),
                       term = NULL, min_term = 5,
                       order_by_size = FALSE, cluster_slices = FALSE,
                       exclude_words = character(0), max_words = 10,
                       word_cloud_grob_param = list(), fontsize_range = c(4, 16),
                       column_title = NULL, ht_list = NULL, use_raster = TRUE, ...) {

  if(length(col) == 1) col = c("white", rgb(t(col2rgb(col)), maxColorValue = 255))
  col_fun = colorRamp2(seq(0, quantile(mat, 0.95), length = length(col)), col)
  if(!is.null(dend)) {
    ht = Heatmap(mat, col = col_fun,
                 name = "Similarity", column_title = column_title,
                 show_row_names = FALSE, show_column_names = FALSE,
                 cluster_rows = dend, cluster_columns = dend,
                 show_row_dend = TRUE, show_column_dend = FALSE,
                 row_dend_width = unit(4, "cm"),
                 border = "#404040", row_title = NULL,
                 use_raster = use_raster,
                 width = unit(8, "cm"))
    draw(ht)
    return(invisible(NULL))
  } else {
    if(inherits(cl, "try-error")) {
      grid.newpage()
      pushViewport(viewport())
      grid.text("Clustering has an error.")
      popViewport()
      return(invisible(NULL))
    }

    # if(!is.factor(cl)) cl = factor(cl, levels = unique(cl))
    cl = as.vector(cl)
    cl_tb = table(cl)
    cl[as.character(cl) %in% names(cl_tb[cl_tb < min_term])] = 0
    cl = factor(cl, levels = c(setdiff(sort(cl), 0), 0))

    if(order_by_size) {
      cl = factor(cl, levels = c(setdiff(names(sort(table(cl), decreasing = TRUE)), 0), 0))
    }
    # od2 = order.dendrogram(dend_env$dend)
    od2 = unlist(lapply(levels(cl), function(le) {
      l = cl == le
      if(sum(l) <= 1) {
        return(which(l))
      } else {
        mm = mat[l, l, drop = FALSE]
        which(l)[hclust(stats::dist(mm))$order]
      }
    }))
    ht = Heatmap(mat, col = col_fun,
                 name = "Similarity", column_title = column_title,
                 show_row_names = FALSE, show_column_names = FALSE,
                 show_row_dend = FALSE, show_column_dend = FALSE,
                 row_order = od2, column_order = od2,
                 border = "#404040", row_title = NULL,
                 use_raster = use_raster,
                 width = unit(8, "cm"))

    if(is.null(term)) {
      if(is.null(rownames(mat))) {
        draw_word_cloud = FALSE
      } else if(!grepl("^GO:[0-9]+$", rownames(mat)[1])) {
        draw_word_cloud = FALSE
      }
    }

    if(draw_word_cloud) {
      go_id = rownames(mat)
      if(!is.null(term)) {
        if(length(term) != length(go_id)) {
          stop_wrap("Length of `term` should be the same as the nrow of `mat`.")
        }
        id2term = structure(term, names = go_id)
      }
      keywords = tapply(go_id, cl, function(term_id) {
        if(is.null(term)) {
          suppressMessages(suppressWarnings(df <- simplifyEnrichment:::count_word(term_id, exclude_words = exclude_words)))
        } else {
          suppressMessages(suppressWarnings(df <- simplifyEnrichment:::count_word(term = id2term[term_id], exclude_words = exclude_words)))
        }
        df = df[df$freq > 1, , drop = FALSE]
        if(nrow(df) > max_words) {
          df = df[order(df$freq, decreasing = TRUE)[seq_len(max_words)], ]
        }
        df
      })
      keywords = keywords[names(keywords) != "0"]
      keywords = keywords[vapply(keywords, nrow, 0) > 0]

      align_to = split(seq_len(nrow(mat)), cl)
      align_to = align_to[names(align_to) != "0"]
      align_to = align_to[names(align_to) %in% names(keywords)]

      word_cloud_grob_param = word_cloud_grob_param[setdiff(names(word_cloud_grob_param), c("text", "fontsize"))]
      pdf(NULL)
      oe = try({
        gbl <- lapply(names(align_to), function(nm) {
          kw = rev(keywords[[nm]][, 1])
          freq = rev(keywords[[nm]][, 2])
          fontsize = scale_fontsize(freq, rg = c(1, max(10, freq)), fs = fontsize_range)

          lt = c(list(text = kw, fontsize = fontsize), word_cloud_grob_param)
          do.call(simplifyEnrichment:::word_cloud_grob, lt)
        })
        names(gbl) = names(align_to)

        margin = unit(8, "pt")
        gbl_h = lapply(gbl, function(x) convertHeight(grobHeight(x), "cm") + margin)
        gbl_h = do.call(unit.c, gbl_h)

        gbl_w = lapply(gbl, function(x) convertWidth(grobWidth(x), "cm"))
        gbl_w = do.call(unit.c, gbl_w)
        gbl_w = max(gbl_w) + margin

      }, silent = TRUE)
      dev.off()
      if(inherits(oe, "try-error")) {
        stop(oe)
      }

      panel_fun = function(index, nm) {
        pushViewport(viewport())
        grid.rect(gp = gpar(fill = "#DDDDDD", col = NA))
        grid.lines(c(0, 1, 1, 0), c(0, 0, 1, 1), gp = gpar(col = "#AAAAAA"), default.units = "npc")
        pushViewport(viewport(width = unit(1, "npc") - margin, height = unit(1, "npc") - margin))
        gb = gbl[[nm]]
        gb$vp$x = gb$vp$width*0.5
        gb$vp$y = gb$vp$height*0.5
        grid.draw(gb)
        popViewport()
        popViewport()
      }

      ht = ht + rowAnnotation(keywords = anno_link(align_to = align_to, which = "row", panel_fun = panel_fun,
                                                   size = gbl_h, gap = unit(2, "mm"), width = gbl_w + unit(5, "mm"),
                                                   link_gp = gpar(fill = "#DDDDDD", col = "#AAAAAA"), internal_line = FALSE))
    } else {
      if(any(cl == "0")) {
        ht = ht + Heatmap(ifelse(cl == "0", "< 5", ">= 5"), col = c("< 5" = "darkgreen", ">= 5" = "white"), width = unit(1, "mm"),
                          heatmap_legend_param = list(title = "", at = "< 5", labels = "Small clusters"),
                          show_column_names = FALSE)
      }
    }
  }

  return(list(ht,od2))
}

# == title
# Scale font size
#
# == param
# -x A numeric vector.
# -rg The range.
# -fs Range of the font size.
#
# == detaisl
# It is a linear interpolation.
#
# == value
# A numeric vector.
#
# == example
# x = runif(10, min = 1, max = 20)
# # scale x to fontsize 4 to 16.
# scale_fontsize(x)
scale_fontsize = function(x, rg = c(1, 30), fs = c(4, 16)) {
  k = (fs[2] - fs[1])/(rg[2] - rg[1])
  b = fs[2] - k*rg[2]
  y = k*x + b
  y[y < fs[1]] = fs[1]
  y[y > fs[2]] = fs[2]
  round(y)
}



# should be able to accommodate both up and down enriched gene sets
plot_select_sets <- function(container, factors_all, sets_plot, thresh=.05, color_sets=NULL, cl_rows=FALSE) {
  hm_list <- NULL
  for (factor_select in factors_all) {
    factor_name <- paste0('Factor',as.character(factor_select))


    gsea_res <- container[["gsea_res_full"]][[factor_name]]

    res <- data.frame(matrix(1,ncol=length(gsea_res),nrow = length(sets_plot)))
    colnames(res) <- names(gsea_res)
    rownames(res) <- sets_plot

    # populate res with pvalues
    for (i in 1:length(sets_plot)) {
      myset <- sets_plot[i]
      for (j in 1:length(gsea_res)) {
        ct <- names(gsea_res)[j]
        myres <- gsea_res[[j]]
        if (myset %in% myres$pathway) {
          # see if NES is negative
          is_neg <- myres$NES[myres$pathway==myset] < 0
          if (is_neg) {
            res[myset,ct] <- log10(myres$padj[myres$pathway==myset])
          } else {
            res[myset,ct] <- -log10(myres$padj[myres$pathway==myset])
          }
        }
      }
    }

    nrn <- rownames(res)
    # make set names multi line if too long!
    for (j in 1:length(nrn)) {
      nm <- nrn[j]
      max_char <- nchar(nm)
      if (nchar(nm) > 40) {
        # cut at underscore
        u_loc <- stringr::str_locate_all(pattern ='_', nm)[[1]]
        ndx_chop <- max(u_loc[,'start'][u_loc[,'start'] < 40])
        nrn[j] <- paste0(substr(nm,1,ndx_chop),'\n',substr(nm,ndx_chop+1,max_char))
      }
    }
    rownames(res) <- nrn

    # order columns the same as corresponding loadings plot
    col_ordering <- colnames(container[["plots"]][["lds_plots_data"]][[as.character(factor_select)]])
    res <- res[,col_ordering]

    res_disc <- res
    res_disc[res>0] <- 'enriched up'
    res_disc[res<0] <- 'enriched down'
    res_disc[abs(res)<(-log10(.05))] <- 'NS'
    colors = structure(c('#FF3333','#3333FF','#E0E0E0'), names = c("enriched up", "enriched down", "NS"))

    hm_legend1 <- Legend(labels = c("enriched up", "enriched down", "NS"), title = "enr direction", legend_gp = gpar(fill = c('#FF3333','#3333FF','#E0E0E0')))
    hm_legend2 <- Legend(labels = c('padj < 0.05','padj < 0.01','padj < 0.001'), title = "significance", type = "points", pch = c("*","**","***"))
    pd <- packLegend(hm_legend1, hm_legend2, direction = "vertical")

    # cluster rows if specified
    if (cl_rows) {
      clust_ord <- hclust(dist(res), method = "single")$order
      res <- res[clust_ord,]
      res_disc <- res_disc[clust_ord,]
      color_sets <- color_sets[clust_ord]
    }

    # height should depend on number of sets
    myheight <- (3/5) * nrow(res)
    myhmap <- Heatmap(as.matrix(res_disc), name = 'signed -log10(padj)',
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      show_row_dend = FALSE, show_column_dend = FALSE,
                      column_names_gp = gpar(fontsize = 12),
                      col = colors,
                      row_title_gp = gpar(fontsize = 12),
                      column_title = 'Cell Types',
                      column_title_side = "bottom",
                      column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                      row_title = 'Gene Sets',
                      row_title_side = "left",
                      border=TRUE,
                      row_names_side = "right",
                      row_names_gp = gpar(fontsize = 8, col = color_sets),
                      show_heatmap_legend = FALSE,
                      width = unit(10, "cm"),
                      height = unit(myheight, "cm"),
                      cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                        if (abs(res[i,j]) > -log10(.001)) {
                          grid.text('***', x, y, gp = gpar(fontface='bold'))
                        } else if (abs(res[i,j]) > -log10(.01)) {
                          grid.text('**', x, y, gp = gpar(fontface='bold'))
                        } else if (abs(res[i,j]) > -log10(.05)) {
                          grid.text('*', x, y, gp = gpar(fontface='bold'))
                        } else {
                          grid.text('', x, y)
                        }
                      })

    hm_list <- hm_list + myhmap
  }
  draw(hm_list,heatmap_legend_list = pd,
       heatmap_legend_side = 'left',
       legend_grouping = "original")
  return(list(hm_list,pd))
}





























