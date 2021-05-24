

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
                             var_scale_power=.5, batch_var=NULL) {
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
    print(ct)

    datExpr <- container$scale_pb_extra[[ct]]

    # Call the network topology analysis function
    sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "unsigned",)
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
                                   TOMType = "unsigned", networkType = "unsigned", minModuleSize = 15,
                                   reassignThreshold = 0, mergeCutHeight = 0.35,
                                   numericLabels = TRUE, pamRespectsDendro = FALSE,
                                   saveTOMs = FALSE,
                                   verbose = 3)
    # net <- WGCNA::blockwiseModules(datExpr, power = sft_thresh[i], maxBlockSize = 10000,
    #                                TOMType = "unsigned", networkType = "unsigned", minModuleSize = 15,
    #                                reassignThreshold = 0, mergeCutHeight = 0.25,
    #                                numericLabels = TRUE, pamRespectsDendro = FALSE,
    #                                saveTOMs = FALSE,
    #                                verbose = 3)

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
#' @param factor_select numeric The factor to get associated LR-modules for
#' @param sig_thresh numeric The p-value significance threshold to use for module-
#' factor associations and ligand-factor associations (default=0.05)
#' @param percentile_exp_rec numeric The percentile above which the top donors expressing the
#' ligand all must be expressing the receptor (default=0.75)
#' @param show_rec_sig logical Set to TRUE to append a heatmap showing the significance
#' of using the receptor expression in predicting module expression (default=TRUE)
#'
#' @return The project container with the plot added in container$plots$lr_analysis$factor_num
#' @export
compute_LR_interact <- function(container, lr_pairs, factor_select, sig_thresh=0.05,
                                percentile_exp_rec=0.75, show_rec_sig=TRUE) {
  ctypes_use <- container$experiment_params$ctypes_use
  dsc <- container$tucker_results[[1]][,factor_select]

  # get list of ligands with associated expression with the factor
  all_lig <- unique(lr_pairs[,1])
  sig_vecs <- get_significance_vectors(container,factor_select,ctypes_use)
  sig_df <- t(as.data.frame(do.call(rbind, sig_vecs)))
  all_lig_mask <- all_lig %in% rownames(sig_df)
  ligs_use <- all_lig[all_lig_mask]
  ct_sig_ligs <- list()
  for (ct in ctypes_use) {
    ct_sig_ligs[[ct]] <- c()
    for (l in ligs_use) {
      if (sig_df[l,ct] < sig_thresh) {
        ct_sig_ligs[[ct]] <- c(ct_sig_ligs[[ct]],l)
      }
    }
  }

  # for all newly added ligands (ones not in sig_df but in pseudbulk) need to compute
  # lm pvalues to get significance
  lig_rest <- colnames(container$scale_pb_extra[[1]])[!(colnames(container$scale_pb_extra[[1]]) %in% rownames(sig_df))]
  lig_rest <- lig_rest[lig_rest %in% all_lig]
  tmp_ct_pvals <- list()
  for (ct in ctypes_use) {
    pb <- container$scale_pb_extra[[ct]]
    for (l in lig_rest) {
      tmp <- as.data.frame(cbind(dsc,pb[names(dsc),l]))
      colnames(tmp) <- c('dsc','expr')
      lmres <- lm(dsc~expr,data=tmp)
      lmres <- summary(lmres)
      pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
      tmp_ct_pvals[[paste0(l,'_',ct)]] <- pval
    }
  }
  tmp_ct_pvals <- p.adjust(tmp_ct_pvals,method='fdr')
  for (i in 1:length(tmp_ct_pvals)) {
    pval <- tmp_ct_pvals[i]
    l <- strsplit(names(pval),split='_')[[1]][[1]]
    ct <- strsplit(names(pval),split='_')[[1]][[2]]
    if (pval < sig_thresh) {
      ct_sig_ligs[[ct]] <- c(ct_sig_ligs[[ct]],l)
    }
  }

  ligs_test <- unlist(ct_sig_ligs)
  ligs_test <- unique(ligs_test)

  # make vector thats inverse of ct_sig_ligs, so lig_ct (ct where it's expressed)
  ligs_ct_test <- lapply(ligs_test,function(x){
    combos <- c()
    for (ct in ctypes_use) {
      if (x %in% ct_sig_ligs[[ct]]) {
        combos <- c(combos,paste0(x,'_',ct))
      }
    }
    return(combos)
  })
  ligs_ct_test <- unlist(ligs_ct_test)
  myres <- NULL
  aovres <- NULL
  mod_fact_r_all <- list()
  ME_pvals <- list()
  ct_rec_pres <- list()
  for (i in 1:length(ctypes_use)) {
    ct <- ctypes_use[i]

    # check if receptor is expressed in the ctype (for each associated ligand)
    no_scale_pb_extra <- container$no_scale_pb_extra[[ct]]

    rec_pres <- c()
    for (lig in ligs_ct_test) {
      # get top nth percentile of donors expressing the ligand
      ligand <- strsplit(lig,split='_')[[1]][[1]]
      ct_exp <- strsplit(lig,split='_')[[1]][[2]]
      lig_ct_exp <- container$scale_pb_extra[[ct_exp]][,ligand]
      nth_quantile <- quantile(lig_ct_exp, probs = c(percentile_exp_rec))
      d_above <- names(lig_ct_exp)[lig_ct_exp > nth_quantile]

      rec_groups <- lr_pairs[lr_pairs[,1]==ligand,2]
      for (r in rec_groups) {
        r_comps <- strsplit(r,split='_')[[1]]
        # if all receptor components are in expression matrix...
        if (sum(r_comps %in% colnames(no_scale_pb_extra))==length(r_comps)) {
          # how many of top donors have 0 expression for any receptor componenet
          mysum <- sum(rowSums(no_scale_pb_extra[d_above,r_comps,drop=FALSE]==0) > 0)
          if (mysum==0) {
            rec_pres <- c(rec_pres,paste0(lig,'_',r))
          }
        }
      }
    }

    # make or add to df for storing results
    if (is.null(myres)) {
      myres <- data.frame(matrix(ncol=0,nrow=length(rec_pres)))
      rownames(myres) <- rec_pres
      aovres <- data.frame(matrix(1,ncol=0,nrow=length(rec_pres)))
      rownames(aovres) <- rec_pres
    } else {
      rec_already_in <- rec_pres %in% rownames(myres)
      if (sum(rec_already_in)!=length(rec_pres)) {
        l_ct_r_add <- rec_pres[!rec_already_in]
        nr_cur <- nrow(myres)
        from <- nr_cur + 1
        to <- nr_cur + length(l_ct_r_add)
        myres[from:to,] <- 0
        rownames(myres)[from:to] <- l_ct_r_add
        aovres[from:to,] <- 1
        rownames(aovres)[from:to] <- l_ct_r_add
      }
    }

    # store receptor presence in each ctype
    ct_rec_pres[[ct]] <- rec_pres

    # extract stored module eigengenes for the ctype
    MEs <- container$module_eigengenes[[ct]]
    for (j in 1:ncol(MEs)) {
      ME <- MEs[,j]
      names(ME) <- rownames(MEs)

      # calculate significance of association with factor as well as Rsq
      tmp <- as.data.frame(cbind(dsc[names(ME)],ME))
      colnames(tmp) <- c('dsc','eg')
      lmres <- lm(dsc~eg,data=tmp)
      lmres <- summary(lmres)
      pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
      mod_fact_r <- cor(dsc[names(ME)],ME)
      mod_fact_r_all[[paste0(ct,'_',j)]] <- mod_fact_r
      ME_pvals[[paste0(ct,"_",as.character(j))]] <- pval
    }
  }
  # fdr correct ME-factor pvals
  ME_pvals <- p.adjust(ME_pvals,method='fdr')

  # for significant modules, test correlations with significant ligands
  for (i in 1:length(ME_pvals)) {
    pval <- ME_pvals[i]
    ct_mod <- names(ME_pvals[i])
    ct <- strsplit(ct_mod,split='_')[[1]][[1]]
    mod <- as.numeric(strsplit(ct_mod,split='_')[[1]][[2]])
    ME <- container$module_eigengenes[[ct]][,mod]
    names(ME) <- rownames(container$module_eigengenes[[ct]])
    ## if significant...
    if (pval < sig_thresh) { # should add bonferroni correction here too
      # add module column to results df if not already there
      if (!(ct_mod %in% colnames(myres))) {
        myres[,ncol(myres)+1] <- 0
        colnames(myres)[ncol(myres)] <- ct_mod
        aovres[,ncol(aovres)+1] <- 1
        colnames(aovres)[ncol(aovres)] <- ct_mod
      }

      # calculate r with significant ligands (that have receptor(s) all present in the ctype)
      rec_pres <- ct_rec_pres[[ct]]
      for (l_ct_r in rec_pres) {
        l_ct_r_splt <- strsplit(l_ct_r,split='_')[[1]]
        ligand <- l_ct_r_splt[[1]]
        ligand_ct <- l_ct_r_splt[[2]]
        rec <- l_ct_r_splt[3:length(l_ct_r_splt)]
        lig_exp <- container$scale_pb_extra[[ligand_ct]][,ligand]
        lig_mod_cor <- cor(lig_exp[names(ME)],ME)
        myres[l_ct_r,ct_mod] <- lig_mod_cor
        rec_exp <- container$scale_pb_extra[[ct]][,rec,drop=FALSE]


        # test whether receptor levels help with prediction
        tmp <- as.data.frame(cbind(lig_exp[names(ME)],ME,rec_exp[names(ME),]))
        colnames(tmp)[1:2] <- c('lig','eg')
        colnames(tmp)[3:ncol(tmp)] <- sapply(1:length(rec),function(x) {
          paste0('rec_',x)
        })
        lm1 <- lm(eg~lig,data=tmp)
        base_formula <- 'eg ~ lig'
        for (k in 1:length(rec)) {
          base_formula <- paste0(base_formula,' + rec_',k)
        }
        base_formula <- as.formula(base_formula)
        lm2 <- lm(base_formula,data=tmp)
        anova_res <- anova(lm1,lm2)
        anova_pval <- anova_res$`Pr(>F)`[2]
        if (ligand %in% rec && ligand_ct == ct) {
          aovres[l_ct_r,ct_mod] <- 1
        } else {
          aovres[l_ct_r,ct_mod] <- anova_pval
        }

        # ensure receptors not in target module
        tmp <- container$module_genes[[ct]]
        tmp_mod <- tmp[tmp==mod]
        if (sum(rec %in% names(tmp_mod))>0) {
          aovres[l_ct_r,ct_mod] <- 1
        }
      }
    }
  }

  lig_dsc_cor_all <- list()
  for (i in 1:nrow(myres)) {
    l_ct_r <- rownames(myres)[i]
    l_ct_r_splt <- strsplit(l_ct_r,split='_')[[1]]
    ligand <- l_ct_r_splt[[1]]
    ligand_ct <- l_ct_r_splt[[2]]
    lig_exp <- container$scale_pb_extra[[ligand_ct]][,ligand]
    lig_dsc_cor <- cor(dsc[names(lig_exp)],lig_exp)
    lig_dsc_cor_all[[l_ct_r]] <- lig_dsc_cor
  }
  lig_dsc_cor_all <- unlist(lig_dsc_cor_all)
  mod_fact_r_all <- unlist(mod_fact_r_all)
  mod_fact_r_all <- mod_fact_r_all[colnames(myres)]

  # plot results for the factor!
  # dont forget to adjust aov pvals and take -log10
  aovres <- as.matrix(aovres)
  aovres2 <- c(aovres)
  aovres2 <- p.adjust(aovres2,method='fdr')
  aovres2 <- matrix(aovres2,nrow=nrow(aovres),ncol=ncol(aovres))
  colnames(aovres2) <- colnames(aovres)
  rownames(aovres2) <- rownames(aovres)
  aovres2 <- -log10(aovres2)

  hc <- hclust(dist(myres))
  hc <- hc[["order"]]
  vc <- hclust(dist(t(myres)))
  vc <- vc[["order"]]

  col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  la <- ComplexHeatmap::rowAnnotation(lig_dsc_cor = lig_dsc_cor_all,col=list(lig_dsc_cor=col_fun),
                                      show_annotation_name=FALSE)
  ta <- ComplexHeatmap::HeatmapAnnotation(mod_dsc_cor = mod_fact_r_all,col=list(mod_dsc_cor=col_fun),
                                          show_annotation_name=FALSE)
  myhmap1 <- Heatmap(as.matrix(myres), name='ligand-mod cor',
                     row_names_side='left', column_names_side='top',
                     show_row_dend=FALSE,
                     show_column_dend=FALSE,
                     column_names_gp = gpar(fontsize = 8),
                     row_names_gp = gpar(fontsize = 8),
                     col=col_fun,
                     left_annotation=la, top_annotation=ta,
                     row_order=hc, column_order=vc,
                     border=TRUE)

  col_fun2 = colorRamp2(c(0, -log10(.05), max(aovres2)), c("white", "white", "green"))
  col_fun2 = colorRamp2(c(0, -log10(.05), 5), c("white", "white", "green"))
  myhmap2 <- Heatmap(aovres2, name='rec_sig',
                     column_names_side='top',
                     show_row_dend=FALSE,
                     show_column_dend=FALSE,
                     column_names_gp = gpar(fontsize = 8),
                     row_names_gp = gpar(fontsize = 6),
                     col=col_fun2,
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     show_row_names=FALSE,
                     row_order=hc, column_order=vc,
                     top_annotation=ta,
                     border=TRUE)

  if (show_rec_sig) {
    hmlist <- myhmap1 + myhmap2
  } else {
    hmlist <- myhmap1
  }

  container$plots$lr_analysis[[paste0('Factor',factor_select)]] <- hmlist
  container$lr_res_raw[[paste0('Factor',factor_select)]] <- myres
  return(container)
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
#'
#' @return the heatmap plot of enrichment results
#' @export
plot_multi_module_enr <- function(container, ctypes, modules, sig_thresh=0.05, db_use='TF') {
  ct_mod <- sapply(1:length(ctypes), function(x) {paste0(ctypes[x],"_",modules[x])})

  mod_res <- list()
  for (i in 1:length(ctypes)) {
    ct <- ctypes[i]
    mymod <- modules[i]
    mod_pvals <- get_module_enr(container, ct, mymod, db_use=db_use, adjust_pval=FALSE)
    mod_res[[ct_mod[i]]] <- mod_pvals
  }

  # make and populate matrix of all pvals
  myres <- data.frame(matrix(1,ncol=length(modules),nrow=length(mod_res[[1]])))
  colnames(myres) <- ct_mod
  rownames(myres) <- names(mod_res[[1]])
  for (i in 1:length(mod_res)) {
    mod_pv <- mod_res[[i]]
    myres[names(mod_pv),names(mod_res)[i]] <- mod_pv
  }

  # adjust pvals
  myres2 <- c(as.matrix(myres))
  myres2 <- p.adjust(myres2,method='fdr')
  myres2 <- matrix(myres2,ncol=length(modules),nrow=length(mod_res[[1]]))
  rownames(myres2) <- rownames(myres)
  colnames(myres2) <- colnames(myres)

  # limit to just TF sets with a significant result
  myres2 <- myres2[rowSums(myres2<sig_thresh)>0,]

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
  col_fun <- colorRamp2(c(.1, 0), c("white", "green"))
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
  lig_exp <- container$scMinimal_ctype[[lig_ct]]$pseudobulk[,lig]
  MEs <- container[["module_eigengenes"]][[mod_ct]]
  ME <- MEs[,mod]
  names(ME) <- rownames(MEs)

  tmp <- as.data.frame(cbind(dsc[names(ME)],ME,lig_exp[names(ME)]))
  colnames(tmp) <- c('dsc','ME','lig_exp')

  mycor1 <- cor(tmp$dsc,tmp$lig_exp)
  p1 <- ggplot(tmp,aes(x=dsc,y=lig_exp)) +
    geom_point() +
    xlab(paste0('Factor',factor_select,' donor score')) +
    ylab(paste0(lig,' expression in ',lig_ct)) +
    annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
             label=paste0('pearson r = ',round(mycor1,digits=3)))

  mycor2 <- cor(tmp$dsc,tmp$ME)
  p2 <- ggplot(tmp,aes(x=dsc,y=ME)) +
    geom_point() +
    xlab(paste0('Factor ',factor_select,' donor score')) +
    ylab(paste0(mod_ct,'_',mod,' module expression')) +
    annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
             label=paste0('pearson r = ',round(mycor2,digits=3)))

  mycor3 <- cor(tmp$ME,tmp$lig_exp)
  p3 <- ggplot(tmp,aes(x=lig_exp,y=ME)) +
    geom_point() +
    xlab(paste0(lig,' expression in ',lig_ct)) +
    ylab(paste0(mod_ct,'_',mod,' module expression')) +
    annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
             label=paste0('pearson r = ',round(mycor3,digits=3)))

  fig_top <- cowplot::plot_grid(plotlist=list(p1,p2), ncol=2)
  fig_bot <- cowplot::plot_grid(plotlist=list(NULL,p3,NULL), ncol=3, rel_widths = c(.25, .5, .25))
  fig <- cowplot::plot_grid(plotlist=list(fig_top,fig_bot), ncol=1)
  return(fig)
}










