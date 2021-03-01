

get_LR_interact<- function(container,lr_pairs,factor_select) {

  # make sure lr_pairs are unique
  lr_pairs <- unique(lr_pairs)

  ctypes_use <- container$experiment_params$ctypes_use

  # extract significance of all genes in all ctypes and put in list
  sig_vectors <- get_significance_vectors(container,
                                          factor_select, ctypes_use)
  # convert list to df
  sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

  # set 0 pvals to the min nonzero pval and take -log10
  min_nonz <- min(sig_df[sig_df!=0])
  sig_df[sig_df==0] <- min_nonz
  sig_df <- -log10(sig_df)

  # sign sig_df by loading
  ldngs <- container$tucker_results[[2]]
  genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
  ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})
  sr_col <- ldngs[factor_select,]
  tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)
  tmp_casted_num <- tmp_casted_num[rownames(sig_df),colnames(sig_df)]
  neg_mask <- tmp_casted_num < 0
  sig_df[neg_mask] <- sig_df[neg_mask] * -1

  # loop through lr pairs
  perms <- c()
  for (ct1 in ctypes_use) {
    for (ct2 in ctypes_use) {
      perms <- c(perms,paste0(ct1,"_",ct2))
    }
  }
  myres <- data.frame(matrix(nrow=0,ncol=length(perms)))
  colnames(myres) <- perms
  rndx <- 0
  for (i in 1:nrow(lr_pairs)) {
    lr <- lr_pairs[i,]
    lig <- lr[[1]]
    rec <- lr[[2]]
    if (lig %in% rownames(sig_df) && rec %in% rownames(sig_df)) {
      rndx <- rndx + 1
      for (j in 1:ncol(sig_df)) {
        v1 <- sig_df[lig,j]
        for (k in 1:ncol(sig_df)) {
          v2 <- sig_df[rec,k]
          cn <- paste0(colnames(sig_df)[j],'_',colnames(sig_df)[k])
          if (v1 * v2 > 0) {
            vs <- c(v1,v2)
            if (vs[1]==vs[2]) {
              min_ndx <- 1
            } else {
              min_ndx <- which(abs(vs)==min(abs(vs)))
            }
            is_neg <- vs[min_ndx] < 0
            if (length(is_neg)>1) {
              print(is_neg)
              print(vs)
            }
            if (is_neg) {
              myres[rndx,cn] <- -1 * min(abs(vs))
            } else {
              myres[rndx,cn] <- min(abs(vs))
            }
          } else {
            myres[rndx,cn] <- 0
          }
          rownames(myres)[rndx] <- paste0(lig,'_',rec)
        }
      }
    }
  }

  # remove rows lr pairs where all nonsignificant
  sig_mask <- abs(myres) > -log10(.01)
  ndx_keep <- rowSums(sig_mask) > 0
  myres <- myres[ndx_keep,]

  col_fun = colorRamp2(c(-max(abs(myres)), log10(.05), 0, -log10(.05), max(myres)), c("blue", "white", "white", "white", "red"))
  hmap <- Heatmap(as.matrix(myres), name='log10(fdr)\nsigned',
                  col=col_fun,
                  show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  border=TRUE,
                  row_names_side='left',
                  row_title='LR Pairs (ligand_receptor)',
                  column_title='Expressing Cell Type',
                  column_title_side = "bottom",
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 8))

  draw(hmap, padding = unit(c(2, 2, 10, 2), "mm")) # add space for titles
  decorate_heatmap_body("log10(fdr)\nsigned", {
    grid::grid.text(paste0("Factor ",factor_select,' LR Pairs'), y = unit(1, "npc") + unit(2, "mm"), just = "bottom")
  })

  return(hmap)
}



# starting with ligands that were significantly variable and present in subsetted dataset
# if not getting very many hits I can incorporate more in later on
# need to access pseudobulked data that hasn't been scaled to set a threshold for receptor presence
# if we don't have a requirement for receptors to be overexpressed then I should incorporate new
# ones in even if not significantly variable
# will need to have scaled/batch corrected values for the receptors to test whether their values
# help with prediction of target module response
# lr_pairs df needs to have ligands in first column and receptors in second. If multiple receptors required, then
# each receptor should be separated by an underscore
prep_LR_interact <- function(container, lr_pairs, norm_method, scale_factor, var_scale_power, batch_var) {
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
  saved_pb <- list()
  for (ct in container$experiment_params$ctypes_use) {
    saved_pb[[ct]] <- container$scMinimal_ctype[[ct]]$pseudobulk
  }
  container$no_scale_pb <- saved_pb

  container <- scale_variance(container,var_scale_power=var_scale_power)

  container <- apply_combat(container,batch_var=batch_var)

  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  for (ct in container$experiment_params$ctypes_use) {
    print(ct)

    datExpr <- container$scMinimal_ctype[[ct]]$pseudobulk

    # Call the network topology analysis function
    sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  }
  return(container)
}

# see line 2051 lupus analysis v2 for prep of new lr_pairs
sft_thresh <- c(3,3,2,2,2,2,2)
compute_LR_interact <- function(container, lr_pairs, sft_thresh, factor_select, sig_thresh=0.05) {
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
  lig_rest <- colnames(container$scMinimal_ctype[[1]]$pseudobulk)[!(colnames(container$scMinimal_ctype[[1]]$pseudobulk) %in% rownames(sig_df))]
  lig_rest <- lig_rest[lig_rest %in% all_lig]
  for (ct in ctypes_use) {
    pb <- container$scMinimal_ctype[[ct]]$pseudobulk
    for (l in lig_rest) {
      tmp <- as.data.frame(cbind(dsc,pb[names(dsc),l]))
      colnames(tmp) <- c('dsc','expr')
      lmres <- lm(dsc~expr,data=tmp)
      lmres <- summary(lmres)
      pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
      if (pval < (sig_thresh/length(lig_rest))) { # using bonferroni correction here
      # if (pval < sig_thresh) {
        print(l)
        print(ct)
        print('')
        ct_sig_ligs[[ct]] <- c(ct_sig_ligs[[ct]],l)
      }
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
  for (i in 1:length(ctypes_use)) {
    ct <- ctypes_use[i]

    # check if receptor is expressed in the ctype (for each associated ligand)
    no_scale_pb <- container$no_scale_pb[[ct]]

    rec_pres <- c()
    for (lig in ligs_ct_test) {
      # get top nth percentile of donors expressing the ligand
      ligand <- strsplit(lig,split='_')[[1]][[1]]
      ct_exp <- strsplit(lig,split='_')[[1]][[2]]
      lig_ct_exp <- container$scMinimal_ctype[[ct_exp]]$pseudobulk[,ligand]
      nth_quantile <- quantile(lig_ct_exp, probs = c(.75))
      d_above <- names(lig_ct_exp)[lig_ct_exp > nth_quantile]

      rec_groups <- lr_pairs[lr_pairs[,1]==ligand,2]
      for (r in rec_groups) {
        r_comps <- strsplit(r,split='_')[[1]]
        # if all receptor components are in expression matrix...
        if (sum(r_comps %in% colnames(no_scale_pb))==length(r_comps)) {
          # how many of top donors have 0 expression for any receptor componenet
          mysum <- sum(rowSums(no_scale_pb[d_above,r_comps,drop=FALSE]==0) > 0)
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

    datExpr <- container$scMinimal_ctype[[ct]]$pseudobulk
    cor <- WGCNA::cor # use cor() from wgcna
    net <- WGCNA::blockwiseModules(datExpr, power = sft_thresh[i], maxBlockSize = 10000,
                           TOMType = "unsigned", minModuleSize = 15,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           saveTOMs = FALSE,
                           verbose = 3)
    cor <- stats::cor # reset cor function

    MEs <- net$MEs
    col_ndx <- sapply(colnames(MEs),function(x) {
      as.numeric(strsplit(x,split='ME')[[1]][[2]])
    })
    MEs <- MEs[,order(col_ndx)]
    MEs[,1] <- NULL
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

      ## if significant...
      if (pval < sig_thresh) { # should add bonferroni correction here too
        # add module column to results df if not already there
        if (!(paste0(ct,'_',j) %in% colnames(myres))) {
          myres[,ncol(myres)+1] <- 0
          colnames(myres)[ncol(myres)] <- paste0(ct,'_',j)
          aovres[,ncol(aovres)+1] <- 1
          colnames(aovres)[ncol(aovres)] <- paste0(ct,'_',j)
        }

        # calculate r with significant ligands (that have receptor(s) all present in the ctype)
        for (l_ct_r in rec_pres) {
          l_ct_r_splt <- strsplit(l_ct_r,split='_')[[1]]
          ligand <- l_ct_r_splt[[1]]
          ligand_ct <- l_ct_r_splt[[2]]
          rec <- l_ct_r_splt[3:length(l_ct_r_splt)]
          lig_exp <- container$scMinimal_ctype[[ligand_ct]]$pseudobulk[,ligand]
          lig_mod_cor <- cor(lig_exp[names(ME)],ME)
          myres[l_ct_r,paste0(ct,'_',j)] <- lig_mod_cor

          # test whether receptor levels help with prediction
          tmp <- as.data.frame(cbind(lig_exp[names(ME)],ME,datExpr[names(ME),rec]))
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
            aovres[l_ct_r,paste0(ct,'_',j)] <- 1
          } else {
            aovres[l_ct_r,paste0(ct,'_',j)] <- anova_pval
          }
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
    lig_exp <- container$scMinimal_ctype[[ligand_ct]]$pseudobulk[,ligand]
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
                     row_names_gp = gpar(fontsize = 8),
                     col=col_fun2,
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     show_row_names=FALSE,
                     row_order=hc, column_order=vc,
                     border=TRUE)
  hmlist <- myhmap1 + myhmap2

}

# need to reset data by rerunning tensor formation with the given parameters






sft_thresh <- c(3,3,2,2,2,2,2)
get_gene_modules <- function(container,sft_thresh) {
  ctypes_use <- container$experiment_params$ctypes_use
  cor <- WGCNA::cor # use cor() from wgcna

  for (i in 1:length(ctypes_use)) {
    ct <- ctypes_use[i]
    # get scaled expression data for the cell type
    datExpr <- container$scMinimal_ctype[[ct]]$pseudobulk

    # get gene modules
    net <- WGCNA::blockwiseModules(datExpr, power = sft_thresh[i], maxBlockSize = 10000,
                                   TOMType = "unsigned", minModuleSize = 15,
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

# see line 2051 lupus analysis v2 for prep of new lr_pairs
compute_LR_interact_v2 <- function(container, lr_pairs, factor_select, sig_thresh=0.05) {
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
  lig_rest <- colnames(container$scMinimal_ctype[[1]]$pseudobulk)[!(colnames(container$scMinimal_ctype[[1]]$pseudobulk) %in% rownames(sig_df))]
  lig_rest <- lig_rest[lig_rest %in% all_lig]
  tmp_ct_pvals <- list()
  for (ct in ctypes_use) {
    pb <- container$scMinimal_ctype[[ct]]$pseudobulk
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
    no_scale_pb <- container$no_scale_pb[[ct]]

    rec_pres <- c()
    for (lig in ligs_ct_test) {
      # get top nth percentile of donors expressing the ligand
      ligand <- strsplit(lig,split='_')[[1]][[1]]
      ct_exp <- strsplit(lig,split='_')[[1]][[2]]
      lig_ct_exp <- container$scMinimal_ctype[[ct_exp]]$pseudobulk[,ligand]
      nth_quantile <- quantile(lig_ct_exp, probs = c(.75))
      d_above <- names(lig_ct_exp)[lig_ct_exp > nth_quantile]

      rec_groups <- lr_pairs[lr_pairs[,1]==ligand,2]
      for (r in rec_groups) {
        r_comps <- strsplit(r,split='_')[[1]]
        # if all receptor components are in expression matrix...
        if (sum(r_comps %in% colnames(no_scale_pb))==length(r_comps)) {
          # how many of top donors have 0 expression for any receptor componenet
          mysum <- sum(rowSums(no_scale_pb[d_above,r_comps,drop=FALSE]==0) > 0)
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
        lig_exp <- container$scMinimal_ctype[[ligand_ct]]$pseudobulk[,ligand]
        lig_mod_cor <- cor(lig_exp[names(ME)],ME)
        myres[l_ct_r,ct_mod] <- lig_mod_cor
        rec_exp <- container$scMinimal_ctype[[ct]]$pseudobulk[,rec,drop=FALSE]


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
      }
    }
  }

  lig_dsc_cor_all <- list()
  for (i in 1:nrow(myres)) {
    l_ct_r <- rownames(myres)[i]
    l_ct_r_splt <- strsplit(l_ct_r,split='_')[[1]]
    ligand <- l_ct_r_splt[[1]]
    ligand_ct <- l_ct_r_splt[[2]]
    lig_exp <- container$scMinimal_ctype[[ligand_ct]]$pseudobulk[,ligand]
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
                     row_names_gp = gpar(fontsize = 8),
                     col=col_fun2,
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     show_row_names=FALSE,
                     row_order=hc, column_order=vc,
                     border=TRUE)
  hmlist <- myhmap1 + myhmap2

}





















