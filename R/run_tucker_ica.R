
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
  
  # check that tucker_type and rotation_type parameters have been set
  if (is.null(container$experiment_params$tucker_type)  ) {
    stop("Need to set tucker_type parameter in experiment params first. Use set_experiment_params()")
  } else {
    tucker_type <- container$experiment_params$tucker_type
  }
  
  if (is.null(container$experiment_params$rotation_type)  ) {
    stop("Need to set rotation_type parameter in experiment params first. Use set_experiment_params()")
  } else {
    rotation_type <- container$experiment_params$rotation_type
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
  tensor_data <- container$tensor_data
  tucker_res <- tucker_ica_helper(tensor_data, ranks, tucker_type,
                                  rotation_type)
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
#' @param tucker_type character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
#' with sparsity constraints
#' @param rotation_type character Set to 'ica' to perform ICA rotation on resulting donor factor
#' matrix and loadings. Otherwise set to 'varimax' to perform varimax rotation. Set to 'none' to
#' not do any rotation after tucker.
#'
#' @return list of results for tucker decomposition with donor scores matrix in first
#' element and loadings matrix in second element
#' @export
tucker_ica_helper <- function(tensor_data, ranks, tucker_type, rotation_type) {
  # extract tensor and labels
  donor_nm <- tensor_data[[1]]
  gene_nm  <- tensor_data[[2]]
  ctype_nm  <- tensor_data[[3]]
  tnsr <- tensor_data[[4]]

  if (tucker_type=='regular') {
    # run regular tucker
    invisible(utils::capture.output(
      tucker_decomp <- rTensor::tucker(rTensor::as.tensor(tnsr), ranks=ranks)
    ))
  } else if (tucker_type=='sparse') {
    # run sparse tucker
    invisible(utils::capture.output(
      tucker_decomp <- tucker_sparse(rTensor::as.tensor(tnsr), ranks=ranks)
    ))
  }

  gene_by_factors <- tucker_decomp$U[[2]]
  rownames(gene_by_factors) <- gene_nm
  ctype_by_factors <- tucker_decomp$U[[3]]
  rownames(ctype_by_factors) <- ctype_nm
  donor_mat <- tucker_decomp$U[[1]]
  rownames(donor_mat) <- donor_nm

  if (ranks[1]>1) {
    # rotate donors matrix
    if (rotation_type == 'ica') {
      donor_mat <- ica::icafast(donor_mat,ranks[1],center=FALSE,alg='def')$S
      
      # make all vectors length 1 as ICA didn't preserve this
      all_rss <- c()
      for (j in 1:ncol(donor_mat)) {
        rss <- sqrt(sum(donor_mat[,j]**2))
        all_rss <- c(all_rss,rss)
      }
      donor_mat <- sweep(donor_mat,2,all_rss,FUN='/')
    } else if (rotation_type=='varimax') {
      donor_mat <- GPForth(donor_mat, method = 'varimax')[[1]]
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





























