



run_stability_analysis <- function(container, ranks, subset_type, sub_prop=.75, n_iterations=100, tucker_type='regular', sparsity=sqrt(2)) {
  dnr_full <- container$tucker_results[[1]]
  lds_full <- container$tucker_results[[2]]

  donors <- rownames(container[["scMinimal_ctype"]][[1]][["pseudobulk"]])
  n_donors <- length(donors)
  donor_ndx_all <- c(1:n_donors)

  # save full tensor data
  full_tensor_data <- container$tensor_data

  # store dscore results
  # res_list <- list()
  # all_bsamps <- list()
  # saf_cors <- c()
  # all_decomps <- list()

  # for (i in 1:n_iterations) {
  res_list <- mclapply(1:n_iterations, function(x) {
    # sample donors
    if (subset_type=='subset') {
      bsamp <- sample(donor_ndx_all,round(n_donors*sub_prop),FALSE)
    } else if (subset_type=='bootstrap') {
      bsamp <- sample(donor_ndx_all,n_donors,TRUE)
    }
    # all_bsamps[[i]] <- bsamp

    # reduce tensor data to just the train donors
    container[["tensor_data"]][[4]] <- container[["tensor_data"]][[4]][bsamp,,]
    container[["tensor_data"]][[1]] <- container[["tensor_data"]][[1]][bsamp]

    # run tucker and rotation
    container <- run_tucker_ica(container, ranks=ranks,
                                tucker_type = tucker_type, rotation_type = 'ica',
                                sparsity=sparsity)

    donor_mat <- container$tucker_results[[1]]
    ldngs <- container$tucker_results[[2]]
    # all_decomps[[i]] <- donor_mat

    d_max <- get_max_correlations(dnr_full,donor_mat,res_use='dscores')
    l_max <- get_max_correlations(lds_full,ldngs,res_use='loadings')
    # saf_cors <- c(saf_cors,d_max[6])

    # reset tensor data
    container$tensor_data <- full_tensor_data

    return(list(d_max,l_max))

  }, mc.cores=container$experiment_params$ncores)

  stability_results <- do.call(rbind.data.frame, res_list)

  factor_indicator <- sapply(1:n_iterations,function(i) {
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
      # geom_point(alpha=.5) +
      # geom_violin() +
      geom_boxplot() +
      # geom_dotplot(binwidth = .0075, dotsize = .75, method='histodot', binaxis = 'y', stackdir='center') +
      xlab('Factor') +
      ylab('Loadings correlation') +
      ggtitle('Loadings stability') +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(c(0,1))
  } else if (plt_data=='dsc') {
    p <- ggplot(sr, aes(x=as.factor(factor), y=as.numeric(dscores))) +
      # geom_point(alpha=.5) +
      # geom_violin() +
      geom_boxplot() +
      # geom_dotplot(binwidth = .0075, dotsize = .75, method='histodot', binaxis = 'y', stackdir='center') +
      xlab('Factor') +
      ylab('Donor scores correlation') +
      ggtitle('Donor scores stability') +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(c(0,1))
  }
  return(p)
}





