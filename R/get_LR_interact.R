utils::globalVariables(c("myx", "myy"))

#' Prepare data for LR analysis and get soft thresholds to use for gene modules
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param lr_pairs data.frame Data of ligand-receptor pairs. First column should
#' be ligands and second column should be one or more receptors separated by an
#' underscore such as receptor1_receptor2 in the case that multiple receptors are
#' required for signaling.
#' @param norm_method character The normalization method to use on the pseudobulked
#' count data. Set to 'regular' to do standard normalization of dividing by
#' library size. Set to 'trim' to use edgeR trim-mean normalization, whereby counts
#' are divided by library size times a normalization factor. (default='trim')
#' @param scale_factor numeric The number that gets multiplied by fractional counts
#' during normalization of the pseudobulked data (default=10000)
#' @param var_scale_power numeric Exponent of normalized variance that is
#' used for variance scaling. Variance for each gene
#' is initially set to unit variance across donors (for a given cell type).
#' Variance for each gene is then scaled by multiplying the unit scaled values
#' by each gene's normalized variance (where the effect of the mean-variance
#' dependence is taken into account) to the exponent specified here.
#' If NULL, uses var_scale_power from container$experiment_params. (default=.5)
#' @param batch_var character A batch variable from metadata to remove (default=NULL)
#'
#' @return The project container with added container$scale_pb_extra slot that contains
#' the tensor with additional ligands and receptors
#' @export
prep_LR_interact <- function(container, lr_pairs, norm_method='trim', scale_factor=10000,
                             var_scale_power=0.5, batch_var=NULL) {
  # store original pseudobulk matrices because they will be altered
  orig_pb <- list()
  for (ct in container$experiment_params$ctypes_use) {
    orig_pb[[ct]] <- container$scMinimal_ctype[[ct]]$pseudobulk
  }

  container <- get_pseudobulk(container)
  container <- normalize_pseudobulk(container, method=norm_method, scale_factor=scale_factor)

  # reduce to vargenes + ligands and receptors
  all_lig <- unique(lr_pairs[,1])
  l_mask <- all_lig %in% colnames(container$scMinimal_ctype[[1]]$pseudobulk)
  lig_add <- all_lig[l_mask]

  all_rec <- lapply(lr_pairs[,2], function(x) {
    return(strsplit(x,split='_')[[1]])
  })
  all_rec <- unlist(all_rec)
  all_rec <- unique(all_rec)
  r_mask <- all_rec %in% colnames(container$scMinimal_ctype[[1]]$pseudobulk)
  rec_add <- all_rec[r_mask]

  # adding rest of ligands and receptors to the pseudobulked data
  vargenes <- container$all_vargenes
  vargenes <- unique(c(vargenes,lig_add,rec_add))
  for (ct in container$experiment_params$ctypes_use) {
    scMinimal <- container$scMinimal_ctype[[ct]]
    scMinimal$pseudobulk <- scMinimal$pseudobulk[,vargenes]
  }

  # need to save pseudobulked normalized data before scaling!!
  no_scale_pb_extra <- list()
  for (ct in container$experiment_params$ctypes_use) {
    no_scale_pb_extra[[ct]] <- container$scMinimal_ctype[[ct]]$pseudobulk
  }
  container$no_scale_pb_extra <- no_scale_pb_extra

  container <- scale_variance(container,var_scale_power=var_scale_power)

  if (!is.null(batch_var)) {
    container <- apply_combat(container,batch_var=batch_var)
  }

  # put new scaled pb data with added genes separate slot
  scale_pb_extra <- list()
  for (ct in container$experiment_params$ctypes_use) {
    scale_pb_extra[[ct]] <- container$scMinimal_ctype[[ct]]$pseudobulk
  }
  container$scale_pb_extra <- scale_pb_extra

  # restore original pseudobulk data in its correct slot
  for (ct in container$experiment_params$ctypes_use) {
    container$scMinimal_ctype[[ct]]$pseudobulk <- orig_pb[[ct]]
  }

  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  for (ct in container$experiment_params$ctypes_use) {

    datExpr <- container$scale_pb_extra[[ct]]

    # Call the network topology analysis function
    sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed",)
  }
  return(container)
}

#' Compute WGCNA gene modules for each cell type
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param sft_thresh numeric A vector indicating the soft threshold to use for
#' each cell type. Length should be the same as container$experiment_params$ctypes_use
#'
#' @return The project container with gene modules added
#' @export
get_gene_modules <- function(container,sft_thresh) {
  ctypes_use <- container$experiment_params$ctypes_use
  cor <- WGCNA::cor # use cor() from wgcna

  for (i in 1:length(ctypes_use)) {
    ct <- ctypes_use[i]
    # get scaled expression data for the cell type
    datExpr <- container$scale_pb_extra[[ct]]

    # get gene modules, used to use cut height of .25
    net <- WGCNA::blockwiseModules(datExpr, power = sft_thresh[i], maxBlockSize = 10000,
                                   TOMType = "signed", networkType = "signed", minModuleSize = 15,
                                   reassignThreshold = 0, mergeCutHeight = 0.25,
                                   numericLabels = TRUE, pamRespectsDendro = FALSE,
                                   saveTOMs = FALSE,
                                   verbose = 3)

    MEs <- net$MEs
    col_ndx <- sapply(colnames(MEs),function(x) {
      as.numeric(strsplit(x,split='ME')[[1]][[2]])
    })
    MEs <- MEs[,order(col_ndx)]
    MEs[,1] <- NULL

    # store modules
    container$module_eigengenes[[ct]] <- MEs
    container$module_genes[[ct]] <- net$colors

  }
  cor <- stats::cor # reset cor function

  return(container)
}



#' Compute and plot the LR interactions for one factor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param lr_pairs data.frame Data of ligand-receptor pairs. First column should
#' be ligands and second column should be one or more receptors separated by an
#' underscore such as receptor1_receptor2 in the case that multiple receptors are
#' required for signaling.
#' @param sig_thresh numeric The p-value significance threshold to use for module-
#' factor associations and ligand-factor associations (default=0.05)
#' @param percentile_exp_rec numeric The percentile above which the top donors expressing the
#' ligand all must be expressing the receptor (default=0.75)
#' @param add_ld_fact_sig logical Set to TRUE to append a heatmap showing significance
#' of associations between each ligand hit and each factor (default=TRUE)
#' @param ncores numeric The number of cores to use (default=container$experiment_params$ncores)
#'
#' @return The results heatmap(s)
#' @export
compute_LR_interact <- function(container, lr_pairs, sig_thresh=0.05,
                                percentile_exp_rec=0.75, add_ld_fact_sig=TRUE, 
                                ncores=container$experiment_params$ncores) {
  all_eg <- container[["module_eigengenes"]]
  all_lig <- unique(lr_pairs[,1])
  ctypes_use <- container$experiment_params$ctypes_use

  ## make matrix to store results
  lig_ct_rec_names <- sapply(ctypes_use,function(x) {
    tmp <- sapply(1:nrow(lr_pairs),function(y) {
      return(paste0(lr_pairs[y,1],'_',x,'_',lr_pairs[y,2]))
    })
  })

  mod_ct_names <- sapply(ctypes_use,function(x) {
    tmp <- sapply(1:length(all_eg[[x]]),function(y) {
      return(paste0(x,'_m',y))
    })
  })
  mod_ct_names <- unlist(mod_ct_names)
  names(mod_ct_names) <- NULL
  myres_mat <- matrix(NA,nrow=length(lig_ct_rec_names),ncol=length(mod_ct_names))
  colnames(myres_mat) <- mod_ct_names
  rownames(myres_mat) <- lig_ct_rec_names

  ## loop through lig_ct_rec combos
  myres <- sccore::plapply(1:length(lig_ct_rec_names), function(lcr_ndx) {
    ct_mod_sig <- c()
    lig_ct_rec <- lig_ct_rec_names[lcr_ndx]
    split_name <- strsplit(lig_ct_rec,split='_')[[1]]
    lig <- split_name[[1]]
    source_ct <- split_name[[2]]
    n_rec_comps <- length(split_name) - 2
    rec_elements <- split_name[3:(2+n_rec_comps)]

    # getting ligand expression in source ctype
    if (!(lig %in% colnames(container$scale_pb_extra[[source_ct]]))) {
      return(NA)
    }

    lig_ct_exp <- container$scale_pb_extra[[source_ct]][,lig]

    # checking ligand expressed to some minimal extent in n% of donors
    if (sum(lig_ct_exp>.2)<(.01*length(lig_ct_exp))) {
      return(NA)
    }

    # # to simply check that ligand is expressed in any donor
    # if (sum(lig_ct_exp!=0)==0) {
    #   return(NA)
    # }

    # loop through target ctypes
    for (target_ct in ctypes_use) {
      if (target_ct == source_ct) {
        next
      }

      # check if rec elements in data
      counts <- container$scMinimal_ctype[[target_ct]]$count_data
      if (sum(rec_elements %in% rownames(counts))!=length(rec_elements)) {
        return(NA)
      }

      # check if rec elements are all present in target ct
      rec_pres <- check_rec_pres(container,lig_ct_exp,rec_elements,target_ct,percentile_exp_rec)

      if (!rec_pres) {
        next
      }

      MEs <- all_eg[[target_ct]]

      # loop through modules for target ct and calculate associations
      for (mod_ndx in 1:ncol(MEs)) {
        tmp <- cbind.data.frame(MEs[names(lig_ct_exp),mod_ndx],lig_ct_exp)
        colnames(tmp) <- c('ME','l_exp')
        lmres <- lm(ME~l_exp,data=tmp)
        lmres <- summary(lmres)
        pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
        ct_mod_sig[paste0(lig,"_",source_ct,'_',target_ct,"_m",mod_ndx)] <- pval
      }
    }
    if (length(ct_mod_sig)==0) {
      return(NA)
    } else {
      return(ct_mod_sig)
    }
  }, mc.preschedule=TRUE,n.cores=ncores,progress=TRUE)

  # copy results and name
  myres_saved <- myres

  # unlist and remove duplicates
  myres2 <- unlist(myres)
  names_keep <- unique(names(myres2))
  myres2 <- myres2[names_keep]
  myres2 <- myres2[!is.na(myres2)]

  # adjust p-values
  myres2 <- p.adjust(myres2,method='fdr')

  # fill results matrix with adjusted p-values
  names(myres) <- lig_ct_rec_names
  for (i in 1:length(myres)) {
    lig_ct_rec <- names(myres)[i]
    if (!is.na(myres[[i]])) {
      for (j in 1:length(myres[[i]])) {
        lig_source_target_mod <- names(myres[[i]][j])
        split_nm <- strsplit(lig_source_target_mod,split='_')[[1]]
        mod_nm <- paste0(split_nm[[3]],"_",split_nm[[4]])
        myres_mat[lig_ct_rec,mod_nm] <- myres2[lig_source_target_mod]
      }
    }
  }

  # make all na to be 1
  myres_mat[is.na(myres_mat)] <- 1

  # store raw results
  container$lr_res <- myres_mat

  # reduce to rows/columns with at least one significant hit
  myres_mat <- myres_mat[rowSums(myres_mat<sig_thresh)>0,]
  myres_mat <- myres_mat[,colSums(myres_mat<.001)>0]

  # log transform values
  myres_mat <- -log10(myres_mat)

  # get split indices
  rs <- sapply(rownames(myres_mat),function(x){
    strsplit(x,split='_')[[1]][[2]]
  })
  cs <- sapply(colnames(myres_mat),function(x){
    strsplit(x,split='_')[[1]][[1]]
  })

  # if have same source and target ctypes, order them same
  rs <- factor(rs,levels=unique(rs))
  cs <- factor(cs,levels=levels(rs))

  # put na values back where source ct == target ct
  myres_mat[outer(rs, cs, "==")] <- NA

  # make new rownames without source ctype in middle
  new_rnames <- sapply(rownames(myres_mat),function(x){
    lig_ct_rec_name <- strsplit(x,split='_')[[1]]
    lig <- lig_ct_rec_name[[1]]
    source <- lig_ct_rec_name[[2]]
    n_rec_comps <- length(lig_ct_rec_name) - 2
    rec_elements <- lig_ct_rec_name[3:(2+n_rec_comps)]
    rec_nm <- rec_elements[1]
    if (n_rec_comps>1) {
      for (j in 2:length(rec_elements)) {
        rec_nm <- paste0(rec_nm,"_",rec_elements[j])
      }
    }
    return(paste0(lig,"_",rec_nm))
  })
  names(new_rnames) <- NULL

  # make heatmap
  col_fun = colorRamp2(c(0, -log10(.1), 10), c("white", "white", "red"))
  myhmap1 <- Heatmap(as.matrix(myres_mat), name='lig_mod -log10(padj)',
                     row_names_side='left', column_names_side='top',
                     show_row_dend=FALSE,
                     show_column_dend=FALSE,
                     column_names_gp = gpar(fontsize = 9),
                     row_names_gp = gpar(fontsize = 9),
                     col=col_fun,
                     border=TRUE,
                     row_split = rs,
                     column_split = cs,
                     cluster_row_slices = FALSE,
                     cluster_column_slices = FALSE,
                     na_col = "gray", column_names_rot = 55,
                     row_labels=new_rnames)

  if (add_ld_fact_sig) {
    fact_res <- matrix(nrow=nrow(myres_mat),ncol=ncol(container$tucker_results[[1]]))
    colnames(fact_res) <- sapply(1:ncol(container$tucker_results[[1]]),function(x){
      paste0('Factor_',x)
    })
    for (i in 1:ncol(container$tucker_results[[1]])) {
      for (j in 1:nrow(myres_mat)) {
        lig <- strsplit(rownames(myres_mat)[j],split='_')[[1]][[1]]
        ct <- strsplit(rownames(myres_mat)[j],split='_')[[1]][[2]]

        lig_ct_exp <- container$scale_pb_extra[[ct]][,lig]

        if (sum(lig_ct_exp!=0)==0) {
          pval <- NA
        } else {
          tmp <- cbind.data.frame(container$tucker_results[[1]][names(lig_ct_exp),i],lig_ct_exp)
          colnames(tmp) <- c('dsc','l_exp')
          lmres <- lm(dsc~l_exp,data=tmp)
          lmres <- summary(lmres)
          pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
          fact_res[j,i] <- pval
        }
      }
    }
    # adjust p-values and log-transform
    fact_res2 <- matrix(p.adjust(fact_res,method='fdr'),ncol=ncol(fact_res),nrow=nrow(fact_res))
    colnames(fact_res2) <- colnames(fact_res)
    fact_res2[is.na(fact_res2)] <- 1
    fact_res2 <- fact_res2[,colSums(fact_res2<sig_thresh)>0]
    fact_res2 <- -log10(fact_res2)
    col_fun = colorRamp2(c(0, -log10(.1), 10), c("white", "white", "purple"))
    myhmap2 <- Heatmap(fact_res2, name='lig_factor -log10(padj)',
                       show_row_dend=FALSE,
                       show_column_dend=FALSE,
                       cluster_columns = FALSE,
                       column_names_gp = gpar(fontsize = 9),
                       row_names_gp = gpar(fontsize = 9),
                       col=col_fun, column_names_rot = 55,
                       border=TRUE)
    myhmap1 <- myhmap1 + myhmap2
  }
  return(myhmap1)
}


#' Helper function to check whether receptor is present in target cell type
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param lig_ct_exp numeric Scaled expression for a ligand in the source cell type
#' @param rec_elements character One or more components of a receptor complex
#' @param target_ct character The name of the target cell type
#' @param percentile_exp_rec numeric The percentile of ligand expression above which
#' all donors need to have at least 5 cells expressing the receptor.
#'
#' @return A logical indicating whether receptor is present or not
check_rec_pres <- function(container,lig_ct_exp,rec_elements,target_ct,percentile_exp_rec) {
  nth_quantile <- quantile(lig_ct_exp, probs = c(percentile_exp_rec))
  d_above <- names(lig_ct_exp)[lig_ct_exp > nth_quantile]

  meta <- container$scMinimal_ctype[[target_ct]]$metadata
  counts <- container$scMinimal_ctype[[target_ct]]$count_data
  for (d in d_above) {
    cells_keep <- rownames(meta)[meta$donors==d]
    counts_sub <- counts[rec_elements,cells_keep,drop=FALSE]
    express_cell_counts <- colSums(counts_sub>0)
    num_cells_expressing <- sum(express_cell_counts==length(rec_elements))
    if (num_cells_expressing < 5) {
      return(FALSE)
    }
  }
  return(TRUE)
}




#' Compute gene sets that are enriched within specified gene co-regulatory modules.
#' Uses a hypergeometric test for over-representation. Used in plot_multi_module_enr.
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ctype character The name of cell type for the cell type module to test
#' @param mod_select numeric The module number for the cell type module to test
#' @param db_use character The database of gene sets to use. Database
#' options include "GO", "Reactome", "KEGG", "BioCarta", "Hallmark", "TF", and
#' "immuno". More than one database can be used. (default="GO")
#' @param adjust_pval logical Set to TRUE to apply FDR correction (default=TRUE)
#'
#' @return p-values for the tested gene sets
#' @export
get_module_enr <- function(container,ctype,mod_select,db_use='GO',adjust_pval=TRUE) {
  mod_genes_all <- container$module_genes[[ctype]] #vector of cluster assignments for each gene
  all_genes <- names(mod_genes_all)
  total_num_genes <- length(all_genes)
  mod_genes <- names(mod_genes_all)[mod_genes_all==mod_select]

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
    } else if (db == "TF") {
      # select the BioCarts gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C3", subcategory = "TFT:GTRD"))
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C3", subcategory = "TFT:TFT_Legacy"))

      # limit it to just sets ending in 'TARGET_GENES'
      target_gene_label <- sapply(m_df$gs_name, function(x) {
        return(grepl('TARGET_GENES', x, fixed = TRUE))
      })
      m_df <- m_df[target_gene_label,]

    } else if (db == "immuno") {
      # select the BioCarts gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C7"))
    }
  }

  my_pathways = split(m_df$gene_symbol, f = m_df$gs_name)
  # my_pathways <- split(m_df$gene_symbol, f = m_df$gs_exact_source)


  pvals <- c()
  for (i in 1:length(my_pathways)) {
    pth <- my_pathways[[i]]
    pth_name <- names(my_pathways)[i]

    # A: total num genes in pathway in tmp_casted_num
    pth_in_df <- unique(pth[which(pth %in% all_genes)])
    num_pth_in_df <- length(pth_in_df)

    # if set is too small continue
    if (num_pth_in_df < 15) {
      next
    }

    # B: number of genes from A in mod_genes
    num_in_sig <- sum(pth_in_df %in% mod_genes)

    # compute pvalue
    pval <- stats::phyper(num_in_sig-1, num_pth_in_df, total_num_genes - num_pth_in_df,
                          length(mod_genes), lower.tail = FALSE) # I double checked this is right
    pvals[pth_name] <- pval
  }

  if (adjust_pval) {
    pvals <- p.adjust(pvals,method='fdr')
  }

  return(pvals)
}

#' Generate gene set x ct_module heatmap showing significantly enriched sets
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ctypes character A vector of cell type names corresponding to the module
#' numbers in mod_select, specifying the modules to compute enrichment for
#' @param modules numeric A vector of module numbers corresponding to the cell
#' types in ctype, specifying the modules to compute enrichment for
#' @param sig_thresh numeric P-value threshold for results to include. Only shows
#' a given gene set if at least one module has a result lower than the threshold.
#' (default=0.05)
#' @param db_use character The database of gene sets to use. Database
#' options include "GO", "Reactome", "KEGG", "BioCarta", "Hallmark", "TF", and
#' "immuno". More than one database can be used. (default="GO")
#' @param max_plt_pval max pvalue shown on plot, but not used to remove rows like
#' sig_thresh (default=.1)
#' @param h_w numeric Vector specifying height and width (defualt=NULL)
#'
#' @return the heatmap plot of enrichment results
#' @export
plot_multi_module_enr <- function(container, ctypes, modules, sig_thresh=0.05, db_use='TF', max_plt_pval=.1, h_w=NULL) {
  ct_mod <- sapply(1:length(ctypes), function(x) {paste0(ctypes[x],"_",modules[x])})

  mod_res <- list()
  for (i in 1:length(ctypes)) {
    ct <- ctypes[i]
    mymod <- modules[i]
    mod_pvals <- get_module_enr(container, ct, mymod, db_use=db_use, adjust_pval=FALSE)
    mod_res[[ct_mod[i]]] <- mod_pvals
  }

  # get all gene sets tested
  all_gsets <- c()
  for (i in 1:length(mod_res)) {
    gsets <- names(mod_res[[i]])
    all_gsets <- c(all_gsets,gsets)
  }
  all_gsets <- unique(all_gsets)

  # make and populate matrix of all pvals
  myres <- data.frame(matrix(ncol=length(modules),nrow=length(all_gsets)))
  colnames(myres) <- ct_mod
  rownames(myres) <- all_gsets
  for (i in 1:length(mod_res)) {
    mod_pv <- mod_res[[i]]
    myres[names(mod_pv),names(mod_res)[i]] <- mod_pv
  }

  # # make and populate matrix of all pvals
  # myres <- data.frame(matrix(1,ncol=length(modules),nrow=length(mod_res[[1]])))
  # colnames(myres) <- ct_mod
  # rownames(myres) <- names(mod_res[[1]])
  # for (i in 1:length(mod_res)) {
  #   mod_pv <- mod_res[[i]]
  #   myres[names(mod_pv),names(mod_res)[i]] <- mod_pv
  # }

  # adjust pvals
  myres2 <- c(as.matrix(myres))
  myres2 <- p.adjust(myres2,method='fdr')
  # myres2 <- matrix(myres2,ncol=length(modules),nrow=length(mod_res[[1]]))
  myres2 <- matrix(myres2,ncol=ncol(myres),nrow=nrow(myres))
  rownames(myres2) <- rownames(myres)
  colnames(myres2) <- colnames(myres)

  # make NA elements go to 1
  myres2[is.na(myres2)] <- 1

  # limit to just TF sets with a significant result
  ndx_keep <- which(rowSums(myres2<sig_thresh)>0)
  myres2 <- myres2[ndx_keep,,drop=FALSE]

  # myres2 <- myres2[rowSums(myres2<sig_thresh)>0,,drop=FALSE]

  if (nrow(myres2)<1) {
    message('no significant gene sets')
    return(NULL)
  }

  nrn <- rownames(myres2)
  # make set names multi line if too long!
  for (j in 1:length(nrn)) {
    nm <- nrn[j]
    max_char <- nchar(nm)
    if (nchar(nm) > 55) {
      # cut at underscore
      u_loc <- stringr::str_locate_all(pattern ='_', nm)[[1]]
      ndx_chop <- max(u_loc[,'start'][u_loc[,'start'] < 55])
      nrn[j] <- paste0(substr(nm,1,ndx_chop),'\n',substr(nm,ndx_chop+1,max_char))
    }
  }
  rownames(myres2) <- nrn

  # plot
  col_fun <- colorRamp2(c(max_plt_pval, 0), c("white", "green"))
  # myhmap <- Heatmap(myres2,name='adj pval',
  #         show_row_dend=FALSE,
  #         show_column_dend=FALSE,
  #         col=col_fun,
  #         row_names_gp = gpar(fontsize = 6),
  #         border=TRUE,
  #         row_names_side='right',
  #         width = unit(6, "cm"),
  #         column_title = 'Co-expression Modules',
  #         column_title_gp = gpar(fontsize = 12),
  #         column_title_side = "bottom")

  if (!is.null(h_w)) {
    myhmap <- Heatmap(myres2,name='adj pval',
                      show_row_dend=FALSE,
                      show_column_dend=FALSE,
                      col=col_fun,
                      row_names_gp = gpar(fontsize = 6.5),
                      border=TRUE,
                      row_names_side='right',
                      column_title = 'Co-expression Modules',
                      column_title_gp = gpar(fontsize = 12),
                      column_title_side = "bottom",
                      clustering_method_rows = "single",
                      width = unit(h_w[2], "cm"), height = unit(h_w[1], "cm"))
  } else {
    myhmap <- Heatmap(myres2,name='adj pval',
                      show_row_dend=FALSE,
                      show_column_dend=FALSE,
                      col=col_fun,
                      row_names_gp = gpar(fontsize = 6.5),
                      border=TRUE,
                      row_names_side='right',
                      column_title = 'Co-expression Modules',
                      column_title_gp = gpar(fontsize = 12),
                      column_title_side = "bottom",
                      clustering_method_rows = "single")
  }


  return(myhmap)
}



#' Plot trio of associations between ligand expression, module levels, and
#' factor scores
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The factor to use
#' @param mod_ct character The name of the cell type for the corresponding module
#' @param mod numeric The number of the corresponding module
#' @param lig_ct character The name of the cell type where the ligand is expressed
#' @param lig character The name of the ligand to use
#'
#' @return plots of the three associations
#' @export
plot_mod_and_lig <- function(container,factor_select,mod_ct,mod,lig_ct,lig) {

  dsc <- container$tucker_results[[1]][,factor_select]
  lig_exp <- container$scale_pb_extra[[lig_ct]][,lig]
  MEs <- container[["module_eigengenes"]][[mod_ct]]
  ME <- MEs[,mod]
  names(ME) <- rownames(MEs)

  tmp <- as.data.frame(cbind(dsc[names(ME)],ME,lig_exp[names(ME)]))
  colnames(tmp) <- c('dsc','ME','lig_exp')

  lmres <- lm(lig_exp~dsc,data=tmp)
  line_range <- seq(min(tmp$dsc),max(tmp$dsc),.001)
  line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
  line_df <- cbind.data.frame(line_range,line_dat)
  colnames(line_df) <- c('myx','myy')

  mycor1 <- cor(tmp$dsc,tmp$lig_exp)
  p1 <- ggplot(tmp,aes(x=dsc,y=lig_exp)) +
    geom_point(alpha = 0.3,pch=19,size=2) +
    geom_line(data=line_df,aes(x=myx,y=myy)) +
    xlab(paste0('Factor ',factor_select,' donor score')) +
    ylab(paste0(lig,' expression in ',lig_ct)) +
    annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
             label=paste0('pearson r = ',round(mycor1,digits=3))) +
    theme_bw()

  lmres <- lm(ME~dsc,data=tmp)
  line_range <- seq(min(tmp$dsc),max(tmp$dsc),.001)
  line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
  line_df <- cbind.data.frame(line_range,line_dat)
  colnames(line_df) <- c('myx','myy')

  mycor2 <- cor(tmp$dsc,tmp$ME)
  p2 <- ggplot(tmp,aes(x=dsc,y=ME)) +
    geom_point(alpha = 0.3,pch=19,size=2) +
    geom_line(data=line_df,aes(x=myx,y=myy)) +
    xlab(paste0('Factor ',factor_select,' donor score')) +
    ylab(paste0(mod_ct,'_',mod,' module expression')) +
    annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
             label=paste0('pearson r = ',round(mycor2,digits=3))) +
    theme_bw()

  lmres <- lm(ME~lig_exp,data=tmp)
  line_range <- seq(min(tmp$lig_exp),max(tmp$lig_exp),.001)
  line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
  line_df <- cbind.data.frame(line_range,line_dat)
  colnames(line_df) <- c('myx','myy')

  mycor3 <- cor(tmp$ME,tmp$lig_exp)
  p3 <- ggplot(tmp,aes(x=lig_exp,y=ME)) +
    geom_point(alpha = 0.3,pch=19,size=2) +
    geom_line(data=line_df,aes(x=myx,y=myy)) +
    xlab(paste0(lig,' expression in ',lig_ct)) +
    ylab(paste0(mod_ct,'_',mod,' module expression')) +
    annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
             label=paste0('pearson r = ',round(mycor3,digits=3))) +
    theme_bw()

  fig_top <- cowplot::plot_grid(plotlist=list(p1,p2), ncol=2)
  fig_bot <- cowplot::plot_grid(plotlist=list(NULL,p3,NULL), ncol=3, rel_widths = c(.25, .5, .25))
  fig <- cowplot::plot_grid(plotlist=list(fig_top,fig_bot), ncol=1)
  return(fig)
}















































