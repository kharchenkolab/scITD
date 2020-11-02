
#' Run Tucker decomposition and apply ICA transformation
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ranks numeric The number of donor, gene, and cell type ranks, respectively,
#' to decompose to using Tucker decomposition. If NULL, uses ranks in container$experiment_params
#' field. (default=NULL)
#' @param shuffle logical If TRUE, randomly shuffles cell to donor linkages, resulting in a random
#' tensor (default=FALSE)
#'
#' @return container with results of the decomposition in container$tucker_results
#' @export
run_tucker_ica <- function(container, ranks=NULL, shuffle=FALSE) {

  # check that var_scale_power has been set if scale_var is TRUE
  if (container$experiment_params$scale_var && is.null(container$experiment_params$var_scale_power)) {
    stop("Need to set variance scaling power parameter, var_scale_power. Use set_experiment_params()")
  }

  if (!is.null(ranks)) {
    # set ranks param in experiment params if specified here
    container <- set_experiment_params(container, ranks = ranks)
  } else {
    # make sure ranks parameter has been set
    if (is.null(container$experiment_params$ranks)) {
      stop("Need to set ranks parameter.")
    }
    ranks <- container$experiment_params$ranks
  }

  # form the tensor for specified cell types
  container <- collapse_by_donors(container, shuffle=shuffle)
  container <- form_tensor(container)

  # run tucker with ica on the tensor
  rotate_modes <- container$experiment_params$rotate_modes
  tensor_data <- container$tensor_data
  tucker_res <- tucker_ica_helper(tensor_data, ranks, rotate_modes)
  container$tucker_results <- tucker_res

  return(container)
}

#' Tucker helper function that actually does the decomposition
#' @importFrom GPArotation GPForth GPFoblq
#'
#' @param tensor_data list The tensor data including donor, gene, and cell type labels
#' as well as the tensor array itself
#' @param ranks numeric The number of donor, gene, and cell type ranks, respectively,
#' to decompose to using Tucker decomposition.
#' @param rotate_modes character The names of the tensor modes to rotate with
#' ICA during Tucker decomposition. Can include 'donors', 'genes', and/or 'ctypes'
#'
#' @return list of results for tucker decomposition with donor scores matrix in first
#' element and loadings matrix in second element
#' @export
tucker_ica_helper <- function(tensor_data, ranks, rotate_modes) {
  # extract tensor and labels
  donor_nm <- tensor_data[[1]]
  gene_nm  <- tensor_data[[2]]
  ctype_nm  <- tensor_data[[3]]
  tnsr <- tensor_data[[4]]

  # # run tucker
  # invisible(utils::capture.output(
  #   tucker_decomp <- rTensor::tucker(rTensor::as.tensor(tnsr), ranks=ranks)
  # ))

  # # trying nonnegative tucker
  # tnsr <- tnsr + min(tnsr)
  # invisible(utils::capture.output(
  #   tucker_decomp <- rTensor::tucker.nonneg(rTensor::as.tensor(tnsr), ranks=ranks)
  # ))

  # run sparse tucker
  invisible(utils::capture.output(
    tucker_decomp <- tucker_sparse(rTensor::as.tensor(tnsr), ranks=ranks)
  ))
  print('yup')
  
  gene_by_factors <- tucker_decomp$U[[2]]
  rownames(gene_by_factors) <- gene_nm
  ctype_by_factors <- tucker_decomp$U[[3]]
  rownames(ctype_by_factors) <- ctype_nm
  donor_mat <- tucker_decomp$U[[1]]
  rownames(donor_mat) <- donor_nm

  if (ranks[1]>1) {
    if ('donors' %in% rotate_modes) {
      # rotate donors matrix by ICA
      # donor_mat <- ica::icafast(donor_mat,ranks[1],center=TRUE,alg='def')$S
      donor_mat <- GPForth(donor_mat, method = 'varimax')[[1]]
      # donor_mat <- GPFoblq(donor_mat)[[1]]
    }
    if ('genes' %in% rotate_modes) {
      # rotate donors matrix by ICA
      # gene_by_factors <- ica::icafast(gene_by_factors,ranks[2],center=TRUE,alg='def')$S
      gene_by_factors <- GPForth(gene_by_factors, method = 'quartimax')[[1]]
    }
    if ('ctypes' %in% rotate_modes) {
      # rotate donors matrix by ICA
      # ctype_by_factors <- ica::icafast(ctype_by_factors,ranks[3],center=TRUE,alg='def')$S
      ctype_by_factors <- GPForth(ctype_by_factors, method = 'quartimax')[[1]]
    }
  }

  # compute kronecker product
  kron_prod <- kronecker(ctype_by_factors,gene_by_factors,make.dimnames = TRUE)

  # generate rotated core tensor unfolded along donor dimension
  core_new <- t(as.matrix(donor_mat)) %*% rTensor::k_unfold(rTensor::as.tensor(tnsr),1)@data %*% kron_prod

  # compute loadings matrix with rotated core tensor
  ldngs <- core_new %*% t(kron_prod)

  return(list(donor_mat,ldngs))
}





























