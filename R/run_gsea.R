
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
run_fgsea <- function(container, factor_select, ctype,
                      db_use="GO") {

  # make sure Tucker has been run
  if (is.null(container$tucker_results)) {
    stop('Run run_tucker_ica() first.')
  }

  # convert ensg symbols to gene symbols
  if (is.null(container$gn_convert)) {
    stop('Gene symbols are not present in your data and no gene name conversion was provided')
  }

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
  
  # trying throwing out 0 value genes...
  ctype_lds <- ctype_lds[ctype_lds!=0]
  

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
                          gseaParam=3)

  fgsea_res <- fgsea_res[order(fgsea_res$padj, decreasing=FALSE),]

  return(fgsea_res)
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
#' @return pvalues for the significantly enriched gene sets
#' @export
run_hypergeometric_gsea <- function(container, factor_select, ctype, up_down,
                                     thresh=0.05, db_use="GO") {

  # make sure Tucker has been run
  if (is.null(container$tucker_results)) {
    stop('Run run_tucker_ica() first.')
  }

  if (is.null(container$gn_convert)) {
    stop('Gene symbols are not present in your data and no gene name conversion was provided')
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
           length(sig_genes), lower.tail = FALSE)
    pvals[pth_name] <- pval
  }
  padj <- p.adjust(pvals,method='fdr')
  # padj <- padj[padj<thresh]
  return(padj)
}

#' Run gsea separately for all cell types of one specified factor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The factor of interest
#' @param method character The method of gsea to use. Can either be "fgsea" or
#' "hypergeometric". (default="fgsea")
#' @param thresh numeric Pvalue significance threshold to use. Will include gene sets in
#' resulting heatmap if pvalue is below this threshold for at least one cell type. (default=0.05)
#' @param db_use character The database of gene sets to use. Database
#' options include "GO", "Reactome", "KEGG", and "BioCarta". More than
#' one database can be used. (default="GO")
#'
#' @return a heatmap plot of the gsea results
#' @export
run_gsea_one_factor <- function(container, factor_select, method="fgsea", thresh=0.05,
                                 db_use="GO") {
  up_sets_all <- list()
  down_sets_all <- list()
  ctypes_use <- container$experiment_params$ctypes_use
  for (ct in ctypes_use) {
    if (method == 'fgsea') {
      fgsea_res <- run_fgsea(container, factor_select=factor_select, ctype=ct,
                             db_use=db_use)

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
    } else if (method == 'hypergeometric') {
      gsea_res_up <- run_hypergeometric_gsea(container, factor_select=factor_select, ctype=ct,
                                             up_down='up', thresh=thresh, db_use=db_use)
      gsea_res_down <- run_hypergeometric_gsea(container, factor_select=factor_select, ctype=ct,
                                               up_down='down', thresh=thresh, db_use=db_use)

      up_sets_all[[ct]] <- gsea_res_up
      down_sets_all[[ct]] <- gsea_res_down
    }
  }


  # plot results
  plot_up <- plot_gsea_hmap(up_sets_all,thresh)
  plot_down <- plot_gsea_hmap(down_sets_all,thresh)

  # add results to container
  factor_name <- paste0('Factor', as.character(factor_select))
  container$gsea_results[[factor_name]] <- list('up'=up_sets_all,
                                                'down'=down_sets_all)
  container$plots$gsea[[factor_name]] <- list('up'=plot_up,
                                              'down'=plot_down)

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
#' @param up_down_sets list The gsea results with either only the up genes sets or down gene sets
#' @param thresh numeric Pvalue threshold to use for including gene sets in the heatmap
#'
#' @return the heatmap plot
#' @export
plot_gsea_hmap <- function(up_down_sets,thresh) {

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
      sum(is.na(ctype_res))
      res[names(ctype_res)[j],ctype] <- ctype_res[j]
    }
  }

  res_plot <- res[rowSums(res<thresh)>0,]

  if (nrow(res_plot) == 0) {
    return(NULL)
  }

  # parse gene set names at first underscore
  rownames(res_plot) <- sapply(rownames(res_plot),function(x) {
    regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
  })

  # cutoff gene set names
  tmp_names <- sapply(rownames(res_plot),function(x) {
    if (nchar(x) > 42) {
      return(paste0(substr(x,1,40),"..."))
    } else {
      return(x)
    }
  })

  # ensure repeats given unique name
  rownames(res_plot) <- make.names(tmp_names, unique = TRUE)

  col_fun <- colorRamp2(c(.05, 0), c("white", "green"))

  myhmap <- Heatmap(as.matrix(res_plot), name = "padj",
                    show_row_dend = FALSE, show_column_dend = FALSE,
                    column_names_gp = gpar(fontsize = 10),
                    col = col_fun, row_title = "Gene Sets",
                    show_row_names = TRUE,
                    row_title_gp = gpar(fontsize = 14),
                    row_names_gp = gpar(fontsize = 6))
  return(myhmap)
}


























