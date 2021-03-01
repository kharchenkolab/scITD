
#' Run Tucker decomposition and apply ICA transformation
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ranks numeric The number of donor, gene, and cell type ranks, respectively,
#' to decompose to using Tucker decomposition.
#' @param tucker_type character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
#' with sparsity constraints (default='regular')
#' @param rotation_type character Set to 'ica' to perform ICA rotation on resulting donor factor
#' matrix and loadings. Otherwise set to 'varimax' to perform varimax rotation. (default='ica')
#'
#' @return container with results of the decomposition in container$tucker_results
#' @export
run_tucker_ica <- function(container, ranks, tucker_type='regular', rotation_type='ica') {

  if (is.null(container$tensor_data)) {
    stop("need to run form_tensor() first")
  }

  # run tucker with ica on the tensor
  tensor_data <- container$tensor_data
  tucker_res <- tucker_ica_helper(tensor_data, ranks, tucker_type, rotation_type)
  container$tucker_results <- tucker_res

  # reorder factors by explained variance
  explained_variances <- c()
  for (i in 1:ranks[1]) {
    exp_var <- get_factor_exp_var(container,i)
    explained_variances[i] <- exp_var
  }
  container$tucker_results[[1]] <- container$tucker_results[[1]][,order(explained_variances,decreasing=TRUE)]
  container$tucker_results[[2]] <- container$tucker_results[[2]][order(explained_variances,decreasing=TRUE),]
  container$exp_var <- explained_variances[order(explained_variances,decreasing=TRUE)]

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

run_gene_centric_tucker <- function(tensor_data, ranks, tucker_type, rotation_type) {

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
      # need to rotate both donor and gene factor matrices
      donor_mat <- ica::icafast(donor_mat,ranks[1],center=FALSE,alg='def')$S

      # make all vectors length 1 as ICA didn't preserve this
      all_rss <- c()
      for (j in 1:ncol(donor_mat)) {
        rss <- sqrt(sum(donor_mat[,j]**2))
        all_rss <- c(all_rss,rss)
      }
      donor_mat <- sweep(donor_mat,2,all_rss,FUN='/')

      gene_by_factors <- ica::icafast(gene_by_factors,ranks[2],center=FALSE,alg='def')$S

      # make all vectors length 1 as ICA didn't preserve this
      all_rss <- c()
      for (j in 1:ncol(gene_by_factors)) {
        rss <- sqrt(sum(gene_by_factors[,j]**2))
        all_rss <- c(all_rss,rss)
      }
      gene_by_factors <- sweep(gene_by_factors,2,all_rss,FUN='/')
    } else if (rotation_type=='varimax') {
      donor_mat <- GPForth(donor_mat, method = 'varimax')[[1]]
      gene_by_factors <- GPForth(gene_by_factors, method = 'varimax')[[1]]
    }
  }

  # compute kronecker product
  kron_prod <- kronecker(ctype_by_factors,donor_mat,make.dimnames = TRUE)

  # generate rotated core tensor unfolded along donor dimension
  core_new <- t(as.matrix(gene_by_factors)) %*% rTensor::k_unfold(rTensor::as.tensor(tnsr),2)@data %*% kron_prod

  # compute loadings matrix with rotated core tensor
  ldngs <- core_new %*% t(kron_prod)

  return(list(gene_by_factors,ldngs))
}


#' Get explained variance for each cell type for one factor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_use numeric The factor to investigate
#'
#' @return explained variance for each cell type in a list
get_factor_exp_var <- function(container, factor_use) {
  tnsr <- rTensor::as.tensor(container$tensor_data[[4]])
  donor_mat <- container$tucker_results[[1]]
  ldngs <- container$tucker_results[[2]]


  recon <- donor_mat[,factor_use,drop=FALSE] %*% ldngs[factor_use,,drop=FALSE]
  recon_tnsr <- rTensor::k_fold(recon,m=1,modes=tnsr@modes)

  # calculate error from using just a single factor
  unexp_var <- (rTensor::fnorm(recon_tnsr - tnsr)**2) / (rTensor::fnorm(tnsr)**2)
  exp_var <- (1 - unexp_var) * 100

  return(exp_var)
}

#' Get explained variance for each cell type for one factor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_use numeric The factor to get variance explained for
#' @param ctype character The cell type to get variance explained for
#'
#' @return explained variance for each cell type in a list
get_ctype_exp_var <- function(container, factor_use, ctype) {
  tnsr <- rTensor::as.tensor(container$tensor_data[[4]])
  donor_mat <- container$tucker_results[[1]]
  ldngs <- container$tucker_results[[2]]


  ctype_ndx <- which(container$tensor_data[[3]]==ctype)
  recon1 <- donor_mat[,factor_use,drop=FALSE] %*% ldngs[factor_use,,drop=FALSE]
  recon1 <- rTensor::k_fold(recon1,m=1,modes=tnsr@modes)

  # The reconstruction should be the original tensor with reconstruction for the one ctype
  recon2 <- tnsr
  recon2[,,ctype_ndx] <- recon2[,,ctype_ndx] - recon1[,,ctype_ndx]


  # calculate error from using just a single factor
  unexp_var <- (rTensor::fnorm(recon2)**2) / (rTensor::fnorm(tnsr)**2)
  exp_var <- (1 - unexp_var) * 100

  return(exp_var)
}

























