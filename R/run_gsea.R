
#' Compute enriched gene sets in a cell type for a factor using fgsea
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The factor of interest
#' @param ctype character The cell type of interest
#' @param thresh numeric Pvalue significance threshold (default=0.05)
#' @param db_use character The database of gene sets to use. Database
#' options include "GO", "Reactome", "KEGG", and "BioCarta". More than
#' one database can be used. (default="GO")
#' @param num_iter numeric The number of random shufflings to perform
#' (default=10000)
#' @param print_res logical If TRUE, prints out the top up/down gene sets
#' (default=TRUE)
#'
#' @return data.frame of the fgsea results limited to genes passing the
#' significance threshold.
#' @export
run_fgsea <- function(container, factor_select, ctype, thresh=0.05,
                      db_use="GO", num_iter=10000, print_res=T) {

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

  fgsea_res <- fgsea::fgsea(pathways = my_pathways,
                          stats = ctype_lds,
                          minSize=15,
                          maxSize=500,
                          nperm=num_iter)

  fgsea_res <- fgsea_res[fgsea_res$padj < thresh,]
  fgsea_res <- fgsea_res[order(fgsea_res$NES, decreasing=T),]

  if (print_res) {
    if (nrow(fgsea_res) >= 10) {
      tmp <- fgsea_res
      tmp$pathway <- sapply(tmp$pathway,function(x) {
        if (nchar(x) > 40) {
          return(paste0(substr(x,1,40),"..."))
        } else {
          return(x)
        }
      })

      print('Enriched in positive loading genes')
      printout1 <- tmp[order(tmp$NES, decreasing=T)[1:10],c('pathway', 'padj', 'NES')]
      print(printout1[printout1$NES>0,])
      print('')
      print('Enriched in negative loading genes')
      printout2 <- tmp[order(tmp$NES, decreasing=F)[1:10],c('pathway', 'padj', 'NES')]
      print(printout2[printout2$NES<0,])
    }
  }

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
#' @param thresh numeric Pvalue significance threshold. Used for both significant
#' gene selection as well as selection of significant gene sets. (default=0.05)
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

  sig_df <- sig_df[,ctype,drop=F]

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
           length(sig_genes), lower.tail = F)
    pvals[pth_name] <- pval
  }
  padj <- p.adjust(pvals,method='fdr')
  padj <- padj[padj<thresh]
  return(padj)
}


#' Run gsea for all cell types and all factors and generate venn diagrams to
#' visualize overlap between enriched genes sets between different cell types
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param method character The method of gsea to use. Can either be "fgsea" or
#' "hypergeometric". (default="fgsea")
#' @param thresh numeric Pvalue significance threshold (default=0.05)
#' @param db_use character The database of gene sets to use. Database
#' options include "GO", "Reactome", "KEGG", and "BioCarta". More than
#' one database can be used. (default="GO")
#' @param num_iter numeric The number of random shufflings to perform. Only
#' applies if using fgsea method (default=10000)
#'
#' @return project container with venn diagrams for each factor added in
#' container$plots slots and gsea results added in container$gsea_results
#' @export
run_gsea_all_factors <- function(container, method="fgsea", thresh=0.05,
                                 db_use="GO", num_iter=10000) {
  for (i in 1:ncol(container$tucker_results[[1]])) {
    up_sets_all <- list()
    down_sets_all <- list()
    ctypes_use <- container$experiment_params$ctypes_use
    for (ct in ctypes_use) {
      if (method == 'fgsea') {
        fgsea_res <- run_fgsea(container, factor_select=i, ctype=ct,
                               thresh=thresh, db_use=db_use, num_iter=num_iter,
                               print_res=F)

        # keep separate track of positive/negative enriched sets
        up_sets <- fgsea_res$pathway[fgsea_res$NES > 0]
        down_sets <- fgsea_res$pathway[fgsea_res$NES < 0]
        up_sets_all[[ct]] <- up_sets
        down_sets_all[[ct]] <- down_sets
      } else if (method == 'hypergeometric') {
        gsea_res_up <- run_hypergeometric_gsea(container, factor_select=i, ctype=ct,
                                            up_down='up', thresh=thresh, db_use=db_use)
        gsea_res_down <- run_hypergeometric_gsea(container, factor_select=i, ctype=ct,
                                               up_down='down', thresh=thresh, db_use=db_use)

        up_sets_all[[ct]] <- names(gsea_res_up)
        down_sets_all[[ct]] <- names(gsea_res_down)
      }
    }

    # plot venn diagram
    if (check_for_all_null(up_sets_all)) {
      up_plot <- venn::venn(x = up_sets_all, box=F, ggplot=T)
    } else {
      up_plot <- NULL
    }

    if (check_for_all_null(down_sets_all)) {
      down_plot <- venn::venn(x = down_sets_all, box=F, ggplot=T)
    } else {
      down_plot <- NULL
    }

    # add results and venn diagram to container
    factor_name <- paste0('Factor', as.character(i))
    container$gsea_results[[factor_name]] <- list('up'=up_sets_all,
                                                  'down'=down_sets_all)
    container$plots$gsea[[factor_name]] <- list('up'=up_plot,
                                                'down'=down_plot)

  }

  return(container)
}

#' Checks if there are any non-NULL elements in a list
#'
#' @param mylist list The list to check
#'
#' @return TRUE if there exist non-NULL elements or FALSE if all are NULL
#' @export
check_for_all_null <- function(mylist) {
  result <- F
  for (i in 1:length(mylist)) {
    if (length(mylist[[i]]) > 0) {
      result <- T
    }
  }
  return(result)
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
#'
#' @return a vector of pathways selectively enriched in the listed cell types
#' @export
get_intersecting_pathways <- function(container, factor_select, these_ctypes_only, up_down) {
  factor_select <- paste0('Factor',factor_select)
  intersect_pathways <- container$gsea_results[[factor_select]][[up_down]][[these_ctypes_only[1]]]
  if (length(these_ctypes_only) > 1) {
    for (i in 2:length(these_ctypes_only)) {
      intersect_pathways <- intersect(intersect_pathways,
                                      container$gsea_results[[factor_select]][[up_down]][[these_ctypes_only[i]]])
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






























