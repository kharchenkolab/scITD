
#' Compute stability of Tucker decomposition by comparing results to those computed
#' on the original dataset downsampled randomly
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param downsample_ratio numeric A decimal between 0-1 that indicates the fraction
#' of cells to retain in the subsampled datasets (default=0.9)
#' @param n_iter numeric The number of downsampling iterations to run (default=500)
#'
#' @return the container with a plot of the results in container$plots$stability_plot
#' and the raw results in container$stability_results
#' @export
run_stability_analysis <- function(container, downsample_ratio=0.9, n_iter=500) {

  # make sure Tucker has been run
  if (is.null(container$tucker_results)) {
    stop('Run run_tucker_ica() first.')
  }

  # set random seed
  RNGkind("L'Ecuyer-CMRG")
  set.seed(container$experiment_params$rand_seed)

  ncores <- container$experiment_params$ncores

  # get tucker loadings to compare subsampled results to
  lds_full <- container$tucker_results[[2]]
  dnr_full <- container$tucker_results[[1]]

  vargenes <- container$all_vargenes
  donors_keep <- container$tensor_data[[1]]
  scMinimal <- subset_scMinimal(container$scMinimal_full, make_copy=TRUE,
                                donors_use = donors_keep, genes_use = vargenes)
  all_avmax_lds_cors <- c()
  all_avmax_dnr_cors <- c()
  stability_results <- mclapply(1:n_iter, function(i) {
    # subsample gene by cell matrix
    cells <- colnames(scMinimal$data_sparse)
    cells_keep <- sample(cells,length(cells) * downsample_ratio)

    # create new container
    scMinimal_sub <- subset_scMinimal(scMinimal, make_copy=TRUE,
                                      cells_use=cells_keep)
    container_sub <- make_new_container(scMinimal_sub,
                                        ctypes_use=container$experiment_params$ctypes_use,
                                        scale_var=container$experiment_params$scale_var,
                                        var_scale_power=container$experiment_params$var_scale_power,
                                        rotate_modes=container$experiment_params$rotate_modes,
                                        ranks=container$experiment_params$ranks,
                                        ncores=container$experiment_params$ncores,
                                        rand_seed=container$experiment_params$rand_seed)

    # create donor av ctype matrices (we already limited it to vargenes above)
    container_sub <- get_ctype_data(container_sub, make_clean=FALSE)

    # run tucker
    container_sub <- run_tucker_ica(container_sub)

    # get correlation of loadings to full dataset
    lds_sub <- container_sub$tucker_results[[2]]
    dnr_sub <- container_sub$tucker_results[[1]]
    avmax_lds_cor <- get_lds_avmax_correlations(lds_full, lds_sub)
    avmax_dnr_cor <- get_dnr_avmax_correlations(dnr_full, dnr_sub)
    return(list(avmax_dnr_cor, avmax_lds_cor))
  }, mc.cores = ncores)

  stability_results <- do.call(rbind.data.frame, stability_results)

  # save results in container
  container$stability_results <- list(stability_results[,1], stability_results[,2])

  # plot results
  stability_plot <- plot_stability_results(container)
  container$plots$stability_plot <- stability_plot

  return(container)
}

#' Compute correlations between loadings from full dataset decomposition to those of
#' the loadings from the subsampled dataset decomposition
#'
#' @param lds_full matrix Loadings matrix from Tucker decomposition on full dataset
#' @param lds_sub matrix Loadings matrix from Tucker decomposition on subsampled dataset
#'
#' @return the average of all maximum correlations between factors of the two decompositions
#' @export
get_lds_avmax_correlations <- function(lds_full, lds_sub) {
  cormat <- cor(t(lds_full), t(lds_sub))
  max_cors <- apply(abs(cormat), MARGIN=2, FUN=max)
  avmax_cor <- mean(max_cors)
  return(avmax_cor)
}

#' Compute correlations between donor scores from full dataset decomposition to those of
#' the donor scores from the subsampled dataset decomposition
#'
#' @param dnr_full matrix Donor scores matrix from Tucker decomposition on full dataset
#' @param dnr_sub matrix Donor scores matrix from Tucker decomposition on subsampled dataset
#'
#' @return the average of all maximum correlations between factors of the two decompositions
#' @export
get_dnr_avmax_correlations <- function(dnr_full, dnr_sub) {
  # make sure donors ordered the same
  dnr_sub <- dnr_sub[rownames(dnr_full),]

  # calculate correlations
  cormat <- cor(dnr_full, dnr_sub)
  max_cors <- apply(abs(cormat), MARGIN=2, FUN=max)
  avmax_cor <- mean(max_cors)
  return(avmax_cor)
}

#' Plot stability results
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#'
#' @return the plot of stability results
#' @export
plot_stability_results <- function(container) {
  sr <- container$stability_results
  all_res <- c(cbind(sr[[1]],sr[[2]]))
  ana_type <- c(cbind(rep('dnr_scores',length(sr[[1]])),
                    rep('loadings',length(sr[[2]]))))
  stability_res <- as.data.frame(cbind(all_res, ana_type))

  p <- ggplot(stability_res, aes(x=as.factor(ana_type), y=as.numeric(as.character(all_res)))) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', method = 'histodot',
                 dotsize = .70, binwidth = .003) +
    scale_y_continuous(breaks=seq(0,1,.1), limits=c(0,1)) +
    scale_x_discrete(labels= c('Donor Scores', 'Loadings')) +
    ylab('Mean Max Correlation (All Factors)') +
    xlab('') +
    ggtitle('Tucker Decomposition \nFull vs Subsampled Dataset') +
    theme(plot.title = element_text(hjust = 0.5))

  return(p)
}
























