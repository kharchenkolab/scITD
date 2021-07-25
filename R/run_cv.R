






run_cv <- function(container, ranks, n_splits=5) {
  # store original tucker results
  dnr_full <- container$tucker_results[[1]]
  lds_full <- container$tucker_results[[2]]

  donors <- rownames(container[["scMinimal_ctype"]][[1]][["pseudobulk"]])
  n_donors <- length(donors)
  donor_ndx_all <- c(1:n_donors)

  # sample donors into splits
  split_ndx <- createFolds(donors, k = n_splits, list = TRUE, returnTrain = TRUE)

  # save full tensor data
  full_tensor_data <- container$tensor_data

  # store dscore results
  all_dsc_res <- list()

  # loop through splits
  for (i in 1:length(split_ndx)) {
    d_train_ndx <- split_ndx[[i]]
    d_test_ndx <- donor_ndx_all[!(donor_ndx_all %in% d_train_ndx)]

    # reduce tensor data to just the train donors
    container[["tensor_data"]][[4]] <- container[["tensor_data"]][[4]][d_train_ndx,,]
    container[["tensor_data"]][[1]] <- container[["tensor_data"]][[1]][d_train_ndx]

    # run tucker and rotation
    container <- run_tucker_ica(container, ranks=ranks,
                                tucker_type = 'regular', rotation_type = 'ica')
    # container <- get_lm_pvals(container)

    donor_mat <- container$tucker_results[[1]]
    ldngs <- container$tucker_results[[2]]

    all_dsc_res[[i]] <- donor_mat

    ldngs <- t(ldngs)

    # swapping variance from loadings to donor scores
    all_rss <- c()
    for (j in 1:ncol(ldngs)) {
      rss <- sqrt(sum(ldngs[,j]**2))
      all_rss <- c(all_rss,rss)
    }
    ldngs <- sweep(ldngs,2,all_rss,FUN='/')

    for (j in 1:ncol(donor_mat)) {
      donor_mat[,j] <- donor_mat[,j] * all_rss[j]
    }

    # ## reducing loadings to only significant genes
    # # need to get gene significance and reshape association results
    # # container <- get_lm_pvals(container)
    # sig_res <- matrix(ncol=ncol(ldngs),nrow=nrow(ldngs))
    # rownames(sig_res) <- rownames(ldngs)
    # for (j in 1:ncol(ldngs)) {
    #   sig_vectors <- get_significance_vectors(container,
    #                                           factor_select=j, container$experiment_params$ctypes_use)
    #   # convert list to df
    #   sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
    #
    #   tmp <- sapply(1:ncol(sig_df),function(x){
    #     ct <- colnames(sig_df)[x]
    #     tmp <- sapply(rownames(sig_df),function(y){
    #       paste0(ct,":",y)
    #     })
    #   })
    #
    #   sig_uf <- c(sig_df)
    #   names(sig_uf) <- tmp
    #   sig_res[names(sig_uf),j] <- sig_uf
    # }
    #
    # # now set non-significant loadings to 0
    # ldngs[sig_res > .05] <- 0
    # ##

    # project scores for test donors
    test_tnsr <- rTensor::as.tensor(full_tensor_data[[4]][d_test_ndx,,])
    test_tnsr_unfold <- rTensor::k_unfold(test_tnsr,1)@data
    dscores_approx <- test_tnsr_unfold %*% ldngs

    ldngs <- t(ldngs)

    ## loop through donor factors to calculate cumulative reconstruction errors
    for (j in 1:ncol(dscores_approx)) {
      factor_use <- 1:j
      # calculate the reconstruction for test donors
      recon <- dscores_approx[,factor_use,drop=FALSE] %*% ldngs[factor_use,,drop=FALSE]
      recon_tnsr <- rTensor::k_fold(recon,m=1,modes=test_tnsr@modes)

      # calculate error from using just a single factor
      unexp_var <- (rTensor::fnorm(recon_tnsr - test_tnsr)**2) / (rTensor::fnorm(test_tnsr)**2)
      exp_var <- (1 - unexp_var) * 100
      print(exp_var)
    }
    # reset tensor data
    container$tensor_data <- full_tensor_data
    print('done')
  }

  # reset original tucker results
  container$tucker_results[[1]] <- dnr_full
  container$tucker_results[[2]] <- lds_full

  return(all_dsc_res)
}


run_cv_cors <- function(container, ranks, n_splits=5) {
  dnr_full <- container$tucker_results[[1]]
  lds_full <- container$tucker_results[[2]]

  donors <- rownames(container[["scMinimal_ctype"]][[1]][["pseudobulk"]])
  n_donors <- length(donors)
  donor_ndx_all <- c(1:n_donors)

  # sample donors into splits
  split_ndx <- createFolds(donors, k = n_splits, list = TRUE, returnTrain = TRUE)

  # save full tensor data
  full_tensor_data <- container$tensor_data

  # store dscore results
  res_list <- list()


  # loop through splits
  for (i in 1:length(split_ndx)) {
    d_train_ndx <- split_ndx[[i]]
    d_test_ndx <- donor_ndx_all[!(donor_ndx_all %in% d_train_ndx)]

    # reduce tensor data to just the train donors
    container[["tensor_data"]][[4]] <- container[["tensor_data"]][[4]][d_train_ndx,,]
    container[["tensor_data"]][[1]] <- container[["tensor_data"]][[1]][d_train_ndx]

    # run tucker and rotation
    container <- run_tucker_ica(container, ranks=ranks,
                                tucker_type = 'regular', rotation_type = 'ica')
    # container <- run_tucker_ica(container, ranks=ranks,
    #                             tucker_type = 'sparse', rotation_type = 'ica',
    #                             sparsity=sqrt(2))

    donor_mat <- container$tucker_results[[1]]
    ldngs <- container$tucker_results[[2]]

    d_max <- get_max_correlations(dnr_full,donor_mat,res_use='dscores')
    l_max <- get_max_correlations(lds_full,ldngs,res_use='loadings')

    # store results
    res_list[[i]] <- list(d_max,l_max)

    # reset tensor data
    container$tensor_data <- full_tensor_data





    # ## trying with projected data
    # ldngs <- t(ldngs)
    #
    # # swapping variance from loadings to donor scores
    # all_rss <- c()
    # for (j in 1:ncol(ldngs)) {
    #   rss <- sqrt(sum(ldngs[,j]**2))
    #   all_rss <- c(all_rss,rss)
    # }
    # ldngs <- sweep(ldngs,2,all_rss,FUN='/')
    #
    # for (j in 1:ncol(donor_mat)) {
    #   donor_mat[,j] <- donor_mat[,j] * all_rss[j]
    # }
    # test_tnsr <- rTensor::as.tensor(full_tensor_data[[4]][d_test_ndx,,])
    # test_tnsr_unfold <- rTensor::k_unfold(test_tnsr,1)@data
    # dscores_approx <- test_tnsr_unfold %*% ldngs
    # rownames(dscores_approx) <- rownames(dnr_full)[d_test_ndx]
    #
    # ldngs <- t(ldngs)
    #
    # d_max <- get_max_correlations(dnr_full,dscores_approx,res_use='dscores')
    # l_max <- get_max_correlations(lds_full,ldngs,res_use='loadings')
    #
    # # store results
    # res_list[[i]] <- list(d_max,l_max)
    #
    # # reset tensor data
    # container$tensor_data <- full_tensor_data
  }

  stability_results <- do.call(rbind.data.frame, res_list)

  factor_indicator <- sapply(1:n_splits,function(i) {
    return(1:ncol(dnr_full))
  })

  stability_results <- cbind.data.frame(stability_results,c(factor_indicator))
  colnames(stability_results) <- c('dscores','ldngs','factor')

  # save results in container
  container$stability_results <- stability_results

  # plot results
  stability_plot_dsc <- plot_stability_results(container,plt_data='dsc')
  stability_plot_lds <- plot_stability_results(container,plt_data='lds')
  container$plots$stability_plot_dsc <- stability_plot_dsc
  container$plots$stability_plot_lds <- stability_plot_lds

  # reset original tucker results
  container$tucker_results[[1]] <- dnr_full
  container$tucker_results[[2]] <- lds_full

  return(container)
}

get_max_correlations <- function(res_full, res_sub, res_use) {
  if (res_use == 'loadings') {

    # get gene_ctype combos present in both decompositions
    gc_use <- intersect(colnames(res_full),colnames(res_sub))
    res_full <- res_full[,gc_use]
    res_sub <- res_sub[,gc_use]

    cormat <- cor(t(res_full), t(res_sub))
  } else if (res_use == 'dscores') {
    # make sure donors ordered the same
    res_full <- res_full[rownames(res_sub),]

    # calculate correlations
    cormat <- cor(res_full, res_sub)
  }

  # calculate average max correlation between factors
  max_cors <- apply(abs(cormat), MARGIN=1, FUN=max)

  return(max_cors)
}

plot_stability_results <- function(container,plt_data) {
  sr <- container$stability_results

  if (plt_data=='lds') {
    p <- ggplot(sr, aes(x=as.factor(factor), y=as.numeric(ldngs))) +
      geom_point(alpha=.5) +
      xlab('Factor') +
      ylab('Loadings correlation') +
      ggtitle('Loadings stability') +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(c(0,1))
  } else if (plt_data=='dsc') {
    p <- ggplot(sr, aes(x=as.factor(factor), y=as.numeric(dscores))) +
      geom_point(alpha=.5) +
      xlab('Factor') +
      ylab('Donor scores correlation') +
      ggtitle('Donor scores stability') +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(c(0,1))
  }


  return(p)
}

run_cv_v2 <- function(container, ranks, n_splits=5) {
  donors <- rownames(container[["scMinimal_ctype"]][[1]][["pseudobulk"]])
  n_donors <- length(donors)
  donor_ndx_all <- c(1:n_donors)

  # sample donors into splits
  split_ndx <- createFolds(donors, k = n_splits, list = TRUE, returnTrain = TRUE)

  # save full tensor data
  full_tensor_data <- container$tensor_data

  # initialize reconstruction data
  all_recon <- list()
  for (i in 1:ranks[1]) {
    all_recon[[i]] <- full_tensor_data[[4]]
  }

  # loop through splits
  for (i in 1:length(split_ndx)) {
    d_train_ndx <- split_ndx[[i]]
    d_test_ndx <- donor_ndx_all[!(donor_ndx_all %in% d_train_ndx)]

    # reduce tensor data to just the train donors
    container[["tensor_data"]][[4]] <- container[["tensor_data"]][[4]][d_train_ndx,,]
    container[["tensor_data"]][[1]] <- container[["tensor_data"]][[1]][d_train_ndx]

    # run tucker and rotation
    container <- run_tucker_ica(container, ranks=ranks,
                                tucker_type = 'regular', rotation_type = 'ica')

    donor_mat <- container$tucker_results[[1]]
    ldngs <- container$tucker_results[[2]]

    ldngs <- t(ldngs)

    # swapping variance from loadings to donor scores
    all_rss <- c()
    for (j in 1:ncol(ldngs)) {
      rss <- sqrt(sum(ldngs[,j]**2))
      all_rss <- c(all_rss,rss)
    }
    ldngs <- sweep(ldngs,2,all_rss,FUN='/')

    for (j in 1:ncol(donor_mat)) {
      donor_mat[,j] <- donor_mat[,j] * all_rss[j]
    }

    # project scores for test donors
    test_tnsr <- rTensor::as.tensor(full_tensor_data[[4]][d_test_ndx,,])
    test_tnsr_unfold <- rTensor::k_unfold(test_tnsr,1)@data
    dscores_approx <- test_tnsr_unfold %*% ldngs

    ldngs <- t(ldngs)

    ## loop through donor factors to calculate cumulative reconstruction errors
    for (j in 1:ncol(dscores_approx)) {
      factor_use <- 1:j
      # calculate the reconstruction for test donors
      recon <- dscores_approx[,factor_use,drop=FALSE] %*% ldngs[factor_use,,drop=FALSE]
      recon_tnsr <- rTensor::k_fold(recon,m=1,modes=test_tnsr@modes)

      # replace data in all_recon with the recon data
      all_recon[[j]][d_test_ndx,,] <- recon_tnsr@data

    }
    # reset tensor data
    container$tensor_data <- full_tensor_data
  }

  # calculate rec errors at each num ranks
  for (j in 1:ranks[1]) {

    # replace data in all_recon with the recon data
    recon_tnsr <- all_recon[[j]]

    # calculate error from using just a single factor
    unexp_var <- (rTensor::fnorm(recon_tnsr - rTensor::as.tensor(full_tensor_data[[4]]))**2) / (rTensor::fnorm(rTensor::as.tensor(full_tensor_data[[4]]))**2)
    exp_var <- (1 - unexp_var) * 100
    print(exp_var)
  }
}




run_cor_bootstrap <- function(container,n_iterations=100) {
  dnr_full <- container$tucker_results[[1]]
  lds_full <- container$tucker_results[[2]]

  donors <- rownames(container[["scMinimal_ctype"]][[1]][["pseudobulk"]])
  n_donors <- length(donors)
  donor_ndx_all <- c(1:n_donors)

  # save full tensor data
  full_tensor_data <- container$tensor_data

  # store dscore results
  res_list <- list()
  test_mxs1 <- c()
  test_mxs2 <- c()
  for (i in 1:n_iterations) {
    # sample donors into splits
    bsamp <- sample(donor_ndx_all,n_donors,TRUE)

    # reduce tensor data to just the train donors
    container[["tensor_data"]][[4]] <- container[["tensor_data"]][[4]][bsamp,,]
    container[["tensor_data"]][[1]] <- container[["tensor_data"]][[1]][bsamp]

    # run tucker and rotation
    container <- run_tucker_ica(container, ranks=ranks,
                                tucker_type = 'regular', rotation_type = 'ica')

    donor_mat <- container$tucker_results[[1]]
    ldngs <- container$tucker_results[[2]]

    d_max <- get_max_correlations(dnr_full,donor_mat,res_use='dscores')
    l_max <- get_max_correlations(lds_full,ldngs,res_use='loadings')

    # store results
    test_mxs1 <- c(test_mxs1,d_max[1])
    test_mxs2 <- c(test_mxs2,d_max[5])

    # reset tensor data
    container$tensor_data <- full_tensor_data

  }


#
#   stability_results <- do.call(rbind.data.frame, res_list)
#
#   factor_indicator <- sapply(1:n_splits,function(i) {
#     return(1:ncol(dnr_full))
#   })
#
#   stability_results <- cbind.data.frame(stability_results,c(factor_indicator))
#   colnames(stability_results) <- c('dscores','ldngs','factor')
#
#   # save results in container
#   container$stability_results <- stability_results
#
#   # plot results
#   stability_plot_dsc <- plot_stability_results(container,plt_data='dsc')
#   stability_plot_lds <- plot_stability_results(container,plt_data='lds')
#   container$plots$stability_plot_dsc <- stability_plot_dsc
#   container$plots$stability_plot_lds <- stability_plot_lds
#
#   # reset original tucker results
#   container$tucker_results[[1]] <- dnr_full
#   container$tucker_results[[2]] <- lds_full
#
#   return(container)
}








