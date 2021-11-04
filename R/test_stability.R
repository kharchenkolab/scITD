
utils::globalVariables(c("ldngs", "dscores"))



#' Test stability of a decomposition by subsampling or bootstrapping donors. Note that
#' running this function will replace the decomposition in the project container
#' with one resulting from the tucker parameters entered here.
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ranks numeric The number of donor, gene, and cell type ranks, respectively,
#' to decompose to using Tucker decomposition.
#' @param tucker_type character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
#' with sparsity constraints (default='regular')
#' @param rotation_type character Set to 'hybrid' to optimize loadings via our hybrid
#' method (see paper for details). Set to 'ica_dsc' to perform ICA rotation
#' on resulting donor factor matrix. Set to 'ica_lds' to optimize loadings by the
#' ICA rotation. (default='hybrid')
#' @param sparsity numeric To use with sparse tucker. Higher indicates more sparse (default=sqrt(2))
#' @param subset_type character Set to either 'subset' or 'bootstrap' (default='subset')
#' @param sub_prop numeric The proportion of donors to keep when using subset_type='subset' (default=.75)
#' @param n_iterations numeric The number of iterations to perform (default=100)
#' @param ncores numeric The number of cores to use (default=container$experiment_params$ncores)
#'
#' @return The project container with the donor scores stability plot in
#' container$plots$stability_plot_dsc and the loadings stability plot in
#' container$plots$stability_plot_lds
#' @export
#' 
#' @examples
#' test_container <- run_stability_analysis(test_container, ranks=c(2,4),
#' tucker_type='regular', rotation_type='hybrid', subset_type='subset', 
#' sub_prop=0.75, n_iterations=5, ncores=1)
run_stability_analysis <- function(container, ranks, tucker_type='regular',
                                   rotation_type='hybrid',  sparsity=sqrt(2),
                                   subset_type='subset', sub_prop=0.75,
                                   n_iterations=100, ncores=container$experiment_params$ncores) {

  ## run tucker with the above parameters in case they changed them
  container <- run_tucker_ica(container, ranks = ranks,
                              tucker_type = tucker_type,
                              rotation_type = rotation_type)
  # pca_unfolded(pbmc_container,2)


  dnr_full <- container$tucker_results[[1]]
  lds_full <- container$tucker_results[[2]]

  donors <- rownames(container[["scMinimal_ctype"]][[1]][["pseudobulk"]])
  n_donors <- length(donors)
  donor_ndx_all <- c(1:n_donors)

  # save full tensor data
  full_tensor_data <- container$tensor_data

  res_list <- sccore::plapply(1:n_iterations, function(x) {
    # sample donors
    if (subset_type=='subset') {
      bsamp <- sample(donor_ndx_all,round(n_donors*sub_prop),FALSE)
    } else if (subset_type=='bootstrap') {
      bsamp <- sample(donor_ndx_all,n_donors,TRUE)
    }

    # reduce tensor data to just the train donors
    container[["tensor_data"]][[4]] <- container[["tensor_data"]][[4]][bsamp,,]
    container[["tensor_data"]][[1]] <- container[["tensor_data"]][[1]][bsamp]

    # run tucker and rotation
    container <- run_tucker_ica(container, ranks=ranks,
                                tucker_type = tucker_type,
                                rotation_type = rotation_type,
                                sparsity=sparsity)

    donor_mat <- container$tucker_results[[1]]
    ldngs <- container$tucker_results[[2]]

    d_max <- get_max_correlations(dnr_full,donor_mat,res_use='dscores')
    l_max <- get_max_correlations(lds_full,ldngs,res_use='loadings')

    # reset tensor data
    container$tensor_data <- full_tensor_data

    return(list(d_max,l_max))

  }, mc.preschedule=TRUE, n.cores=ncores, progress=TRUE)

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


#' Computes the max correlation between each factor of the decomposition done using
#' the whole dataset to each factor computed using the subsampled/bootstrapped dataset
#'
#' @param res_full matrix Either the donor scores or loadings matrix from the original
#' decomposition
#' @param res_sub matrix Either the donor scores or loadings matrix from the new
#' decomposition
#' @param res_use character Can either be 'loadings' or 'dscores' and should correspond
#' with the data matrix used
#'
#' @return a vector of the max correlations for each original factor
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

  # calculate max correlations from original factors to new ones
  max_cors <- apply(abs(cormat), MARGIN=1, FUN=max)

  return(max_cors)
}

#' Generate a plot for either the donor scores or loadings stability test
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param plt_data character Either 'lds' or 'dsc' and indicates which plot to make
#'
#' @return the plot
plot_stability_results <- function(container,plt_data) {
  sr <- container$stability_results

  if (plt_data=='lds') {
    p <- ggplot(sr, aes(x=as.factor(factor), y=as.numeric(ldngs))) +
      geom_boxplot() +
      xlab('Factor') +
      ylab('Loadings correlation') +
      ggtitle('Loadings stability') +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(c(0,1))
  } else if (plt_data=='dsc') {
    p <- ggplot(sr, aes(x=as.factor(factor), y=as.numeric(dscores))) +
      geom_boxplot() +
      xlab('Factor') +
      ylab('Donor scores correlation') +
      ggtitle('Donor scores stability') +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(c(0,1))
  }
  return(p)
}





