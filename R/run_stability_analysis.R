
#' Compute stability of Tucker decomposition by comparing results to those computed
#' on the original dataset downsampled randomly
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ranks numeric The number of donor, gene, and cell type ranks, respectively,
#' to decompose to using Tucker decomposition.
#' @param downsample_ratio numeric A decimal between 0-1 that indicates the fraction
#' of cells to retain in the subsampled datasets (default=0.9)
#' @param n_iter numeric The number of downsampling iterations to run (default=500)
#' @param norm_method character The normalization method to use on the pseudobulked
#' count data. Set to 'regular' to do standard normalization of dividing by
#' library size. Set to 'trim' to use edgeR trim-mean normalization, whereby counts
#' are divided by library size times a normalization factor. (default='trim')
#' @param scale_factor numeric The number that gets multiplied by fractional counts
#' during normalization of the pseudobulked data (default=10000)
#' @param batch_var character A batch variable from metadata to remove (default=NULL)
#' @param scale_var logical TRUE to scale the gene expression variance across donors
#' for each cell type. If FALSE then all genes are scaled to unit variance across
#' donors for each cell type. (default=TRUE)
#' @param var_scale_power numeric Exponent of normalized variance that is
#' used for variance scaling. Variance for each gene
#' is initially set to unit variance across donors (for a given cell type).
#' Variance for each gene is then scaled by multiplying the unit scaled values
#' by each gene's normalized variance (where the effect of the mean-variance
#' dependence is taken into account) to the exponent specified here.
#' If NULL, uses var_scale_power from container$experiment_params. (default=.5)
#' @param tucker_type character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
#' with sparsity constraints (default='regular')
#' @param rotation_type character Set to 'ica' to perform ICA rotation on resulting donor factor
#' matrix and loadings. Otherwise set to 'varimax' to perform varimax rotation. (default='ica')
#'
#' @return the container with a plot of the results in container$plots$stability_plot
#' and the raw results in container$stability_results
#' @export
run_stability_analysis <- function(container, ranks, downsample_ratio=0.9, n_iter=500,
                                   norm_method='trim', scale_factor=10000,
                                   batch_var=NULL, scale_var=TRUE,
                                   var_scale_power=.5, tucker_type='regular',
                                   rotation_type='ica') {

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

  # save vargenes so they don't have to be recalculated
  all_vargenes <- container$all_vargenes

  # subset full counts matrix by donors used in original decomposition
  donors_keep <- container$tensor_data[[1]]
  scMinimal <- subset_scMinimal(container$scMinimal_full, in_place=FALSE,
                                donors_use = donors_keep)

  stability_results <- mclapply(1:n_iter, function(i) {
  # stability_results <- lapply(1:n_iter, function(i) {
    # subsample gene by cell matrix
    cells <- colnames(scMinimal$count_data)
    cells_keep <- sample(cells,length(cells) * downsample_ratio)

    # create new container
    scMinimal_sub <- subset_scMinimal(scMinimal, in_place=FALSE,
                                      cells_use=cells_keep)

    container_sub <- make_new_container(params=container$experiment_params, scMinimal=scMinimal_sub)

    # run through tensor formation (no donor cleaning or new vargenes selection...)
    container_sub <- parse_data_by_ctypes(container_sub)
    container_sub <- clean_data(container_sub, donor_min_cells=0, gene_min_cells=5)
    container_sub$all_vargenes <- all_vargenes[all_vargenes %in% rownames(container_sub$scMinimal_ctype[[1]]$count_data)]
    container_sub <- get_pseudobulk(container_sub)
    container_sub <- normalize_pseudobulk(container_sub, method=norm_method, scale_factor=scale_factor)
    container_sub <- get_normalized_variance(container_sub)
    container_sub <- reduce_to_vargenes(container_sub)
    if (!is.null(batch_var)) {
      invisible(utils::capture.output(
        container_sub <- apply_combat(container_sub,batch_var=batch_var)
      ))
    }
    if (scale_var) {
      container_sub <- scale_variance(container_sub,var_scale_power=var_scale_power)
    }
    container_sub <- stack_tensor(container_sub)

    # run tucker
    container_sub <- run_tucker_ica(container_sub, ranks, tucker_type=tucker_type, rotation_type=rotation_type)

    # get correlation of loadings to full dataset
    lds_sub <- container_sub$tucker_results[[2]]
    dnr_sub <- container_sub$tucker_results[[1]]
    avmax_lds_cor <- get_avmax_correlations(lds_full,lds_sub,res_use='loadings')
    avmax_dnr_cor <- get_avmax_correlations(dnr_full,dnr_sub,res_use='dscores')

    return(list(avmax_dnr_cor, avmax_lds_cor))
  # })
  }, mc.cores = ncores)

  stability_results <- do.call(rbind.data.frame, stability_results)

  # save results in container
  container$stability_results <- list(stability_results[,1], stability_results[,2])

  # plot results
  stability_plot <- plot_stability_results(container)
  container$plots$stability_plot <- stability_plot

  return(container)
}


#' Compute correlations between factors from full dataset decomposition to those of
#' from the subsampled dataset decomposition
#'
#' @param res_full matrix Either loadings or donor scores matrix from
#' Tucker decomposition on full dataset
#' @param res_sub matrix Either loadings or donor scores matrix from
#' Tucker decomposition on subsampled dataset
#' @param res_use character Set to either 'loadings' or 'dscores' to indicate
#' the data type that has been entered
#'
#' @return the average of all maximum correlations between factors of the two decompositions
#' @export
get_avmax_correlations <- function(res_full, res_sub, res_use) {
  if (res_use == 'loadings') {

    # get gene_ctype combos present in both decompositions
    gc_use <- intersect(colnames(res_full),colnames(res_sub))
    res_full <- res_full[,gc_use]
    res_sub <- res_sub[,gc_use]

    cormat <- cor(t(res_full), t(res_sub))
  } else if (res_use == 'dscores') {
    # make sure donors ordered the same
    res_sub <- res_sub[rownames(res_full),]

    # calculate correlations
    cormat <- cor(res_full, res_sub)
  }

  # calculate average max correlation between factors
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
























