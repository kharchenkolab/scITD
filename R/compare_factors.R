

#' Title
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param f_compare numeric The number of the two factors to compare
#' @param direction character Vector specifying the directions to compare for each
#' factor. Can either be 'up' or 'down' for each factor.
#' @param compare_type character Set to either 'same' or 'different' to get
#' genes that are significant in both or only significant in one or the other,
#' respectively.
#' @param sig_thresh numeric The significance p-value threshold to use (default=0.05)
#'
#' @return The comparison heatmap in container$plots$comparisons$<f1>_<f2>
#' @export
compare_factors <- function(container, f_compare, direction, compare_type,
                            sig_thresh=0.05) {
  f1 <- f_compare[1]
  f2 <- f_compare[2]

  ldngs <- container$tucker_results[[2]]

  # break down a factor from the loadings matrix
  genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
  ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})

  # get f1 sig gene loadings
  sr_col <- ldngs[f1,]
  tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)
  sig_vectors <- get_significance_vectors(container,
                                          f1, colnames(tmp_casted_num))
  sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
  sig_df <- sig_df[rownames(tmp_casted_num),colnames(tmp_casted_num)]
  tmp_casted_num1 <- tmp_casted_num[rowSums(sig_df < sig_thresh) > 0,]
  tmp_casted_num1[sig_df[rownames(tmp_casted_num1),colnames(tmp_casted_num1)] > sig_thresh] <- 0
  if (direction[1]=='up') {
    tmp_casted_num1[tmp_casted_num1<0] <- 0
  } else if (direction[1]=='down') {
    tmp_casted_num1[tmp_casted_num1>0] <- 0
  }
  tmp_casted_num1 <- tmp_casted_num1[rowSums(tmp_casted_num1!=0)>0,]

  # get f2 sig gene loadings
  sr_col <- ldngs[f2,]
  tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)
  sig_vectors <- get_significance_vectors(container,
                                          f2, colnames(tmp_casted_num))
  sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
  sig_df <- sig_df[rownames(tmp_casted_num),colnames(tmp_casted_num)]
  tmp_casted_num2 <- tmp_casted_num[rowSums(sig_df < sig_thresh) > 0,]
  tmp_casted_num2[sig_df[rownames(tmp_casted_num2),colnames(tmp_casted_num2)] > sig_thresh] <- 0
  if (direction[2]=='up') {
    tmp_casted_num2[tmp_casted_num2<0] <- 0
  } else if (direction[2]=='down') {
    tmp_casted_num2[tmp_casted_num2>0] <- 0
  }
  tmp_casted_num2 <- tmp_casted_num2[rowSums(tmp_casted_num2!=0)>0,]


  if (compare_type=='same') {
    g_in_both <- intersect(rownames(tmp_casted_num1),rownames(tmp_casted_num2))

    tmp_casted_num1 <- tmp_casted_num1[g_in_both,]
    tmp_casted_num2 <- tmp_casted_num2[g_in_both,]

    # loop through each cell type
    res <- c()
    for (i in 1:ncol(tmp_casted_num1)) {
      c1 <- tmp_casted_num1[,i]
      c2 <- tmp_casted_num2[,i]

      tmp <- (c1!=0) & (c2!=0)

      tmp[tmp] <- 1
      tmp[!tmp] <- 0

      res <- cbind(res,tmp)
    }

    colnames(res) <- colnames(tmp_casted_num1)

    r_order <- hclust(dist(res),method='median')$order
    c_order <- hclust(dist(t(res)),method='median')$order

    res <- res[r_order,c_order]

    # remove rows that are all not DE both
    res <- res[rowSums(res>0)>0,]

    container$factor_comparison <- as.matrix(res)

    res[res==0] <- "not DE both"
    res[res==1] <- "DE both"

    # colors = structure(c('white','red'), names = c("1", "2"))
    colors = structure(c('white','red'), names = c("not DE both", "DE both"))

    if (nrow(res) > 75) {
      callout_ndx <- sample(1:nrow(res),75)
      gene_callouts <- rownames(res)[callout_ndx]
      myannot <- rowAnnotation(callouts = anno_mark(at = callout_ndx, which='row',
                                                    labels = gene_callouts))
      show_genes <- FALSE
    } else {
      show_genes <- TRUE
      myannot <- NULL
    }

    hmap <- Heatmap(res, name='sig both',
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            col = colors, border=TRUE,
            right_annotation=myannot,
            show_row_names=show_genes)

  } else if (compare_type=='different') {
    # get union of gene names
    g_in_both <- union(rownames(tmp_casted_num1),rownames(tmp_casted_num2))

    # make a df to store results
    res <- data.frame(matrix(ncol=ncol(tmp_casted_num1),nrow=length(g_in_both)))
    rownames(res) <- g_in_both
    colnames(res) <- colnames(tmp_casted_num1)

    # put the union genes into the tmp_casted matrices if not already there
    genes_not_df1 <- g_in_both[!(g_in_both %in% rownames(tmp_casted_num1))]
    genes_not_df2 <- g_in_both[!(g_in_both %in% rownames(tmp_casted_num2))]

    df_add1 <- matrix(0,ncol=ncol(tmp_casted_num1),nrow=length(genes_not_df1))
    rownames(df_add1) <- genes_not_df1
    colnames(df_add1) <- colnames(tmp_casted_num1)

    df_add2 <- matrix(0,ncol=ncol(tmp_casted_num2),nrow=length(genes_not_df2))
    rownames(df_add2) <- genes_not_df2
    colnames(df_add2) <- colnames(tmp_casted_num2)

    tmp_casted_num1 <- rbind(tmp_casted_num1,df_add1)
    tmp_casted_num2 <- rbind(tmp_casted_num2,df_add2)

    tmp_casted_num1 <- tmp_casted_num1[rownames(res),]
    tmp_casted_num2 <- tmp_casted_num2[rownames(res),]

    equal_tester <- function(v1,v2) {
      if (v1!=0 && v2==0) {
        return(1)
      } else if (v1==0 && v2!=0) {
        return(-1)
      }
      return(0)
    }

    for (j in 1:ncol(res)) {
      for (i in 1:nrow(res)) {
        res[i,j] <- equal_tester(tmp_casted_num1[i,j],tmp_casted_num2[i,j])
      }
    }

    r_order <- hclust(dist(res),method='median')$order
    c_order <- hclust(dist(t(res)),method='median')$order

    res <- res[r_order,c_order]

    # remove rows that are all not DE both
    res <- res[rowSums(res!=0)>0,]

    container$factor_comparison <- as.matrix(res)

    res[res==0] <- "DE both/neither"
    res[res==1] <- paste0("DE Factor ",f1)
    res[res==-1] <- paste0("DE Factor ",f2)

    colors = structure(c('white','red','blue'), names = c("DE both/neither",
                                                          paste0("DE Factor ",f1),
                                                          paste0("DE Factor ",f2)))

    if (nrow(res) > 75) {
      callout_ndx <- sample(1:nrow(res),75)
      gene_callouts <- rownames(res)[callout_ndx]
      myannot <- rowAnnotation(callouts = anno_mark(at = callout_ndx, which='row',
                                                    labels = gene_callouts))
      show_genes <- FALSE
    } else {
      show_genes <- TRUE
      myannot <- NULL
    }

    hmap <- Heatmap(as.matrix(res), name='diff sig',
                    show_row_dend = FALSE,
                    show_column_dend = FALSE,
                    col = colors, border=TRUE,
                    right_annotation=myannot,
                    show_row_names=show_genes)

  }

  container$plots$comparisons[[paste0(f1,'_',f2)]] <- hmap

  return(container)

}



get_compare_go_enrich <- function(container,ctype,direc,db_use='GO') {
  comp <- container$factor_comparison
  comp <- as.data.frame(comp)

  sig_genes <- rownames(comp)[comp[[ctype]]==direc]

  all_genes <- container$tensor_data[[2]]

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
    pval <- stats::phyper(num_in_sig-1, num_pth_in_df, length(all_genes) - num_pth_in_df,
                          length(sig_genes), lower.tail = FALSE)
    pvals[pth_name] <- pval
  }
  padj <- p.adjust(pvals,method='fdr')

  return(padj)
}


































