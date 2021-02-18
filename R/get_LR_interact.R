

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


