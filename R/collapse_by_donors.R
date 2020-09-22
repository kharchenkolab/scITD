
#' Computes mean expression for each gene in each cell type per donor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param shuffle logical If TRUE shuffles cell to donor linkages (default=FALSE)
#'
#' @return an updated container with the collapsed gene by donor average matrices
#' in the slots of container$scMinimal_ctype
#' @export
collapse_by_donors <- function(container, shuffle=F) {
  for (ct in container$experiment_params$ctypes_use) {
    scMinimal <- container$scMinimal_ctype[[ct]]
    dge_sparse <- Matrix::t(scMinimal$data_sparse)

    if (shuffle) {
      donor_meta <- as.factor(sample(scMinimal$metadata$donors))
    } else {
      donor_meta <- as.factor(scMinimal$metadata$donors)
    }
    donor_means <- get_means(dge_sparse, donor_meta, table(donor_meta))

    # remove first row because it's all NaN
    donor_means <- donor_means[2:nrow(donor_means),]

    scMinimal$data_means <- donor_means
  }
  return(container)
}





