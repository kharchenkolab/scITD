
#' Run Tucker decomposition and apply ICA transformation
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ranks numeric The number of donor factors and gene factors, respectively,
#' to decompose the data into. Since we rearrange the standard output of
#' the Tucker decomposition to be 'donor centric', the number of donor factors will
#' also be the total number of main factors that can be used for downstream analysis.
#' The number of gene factors will only impact the quality of the decomposition.
#' @param tucker_type character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
#' with sparsity constraints. The 'sparse' method is still under development, so we recommend
#' using 'regular'. (default='regular')
#' @param rotation_type character Set to 'hybrid' to optimize loadings via our hybrid
#' method (see paper for details). Set to 'ica_dsc' to perform ICA rotation
#' on resulting donor factor matrix. Set to 'ica_lds' to optimize loadings by the
#' ICA rotation. (default='hybrid')
#' @param sparsity numeric To use with sparse tucker. Higher indicates more sparse (default=sqrt(2))
#'
#' @return container with results of the decomposition in container$tucker_results
#' @export
#'
#' @examples
#' test_container <- run_tucker_ica(test_container,ranks=c(2,4))
run_tucker_ica <- function(container, ranks, tucker_type='regular', rotation_type='hybrid', sparsity=sqrt(2)) {

  if (is.null(container$tensor_data)) {
    stop("need to run form_tensor() first")
  }

  # check that they picked a valid tucker_type and rotation type
  if (!(tucker_type %in% c('regular','sparse'))) {
    stop("tucker_type can only be 'regular' or 'sparse'")
  }
  if (!(rotation_type %in% c('hybrid','ica_dsc','ica_lds'))) {
    stop("rotation_type can only be 'hybrid', 'ica_dsc', or 'ica_lds'")
  }

  # get tensor data
  tensor_data <- container$tensor_data

  # set 3rd ranks position to be the number of cell types
  ranks[3] <- length(tensor_data[[3]])

  # run tucker with rotation on the tensor
  tucker_res <- tucker_ica_helper(tensor_data, ranks, tucker_type, rotation_type, sparsity)
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
#'
#' @param tensor_data list The tensor data including donor, gene, and cell type labels
#' as well as the tensor array itself
#' @param ranks numeric The number of donor and gene factors respectively,
#' to decompose to using Tucker decomposition.
#' @param tucker_type character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
#' with sparsity constraints
#' @param rotation_type character Set to 'hybrid' to optimize loadings via our hybrid
#' method (see paper for details). Set to 'ica_dsc' to perform ICA rotation
#' on resulting donor factor matrix. Set to 'ica_lds' to optimize loadings by the
#' ICA rotation.
#' @param sparsity numeric Higher indicates more sparse
#'
#' @return list of results for tucker decomposition with donor scores matrix in first
#' element and loadings matrix in second element
#' @export
tucker_ica_helper <- function(tensor_data, ranks, tucker_type, rotation_type, sparsity) {

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
      tucker_decomp <- tucker_sparse(rTensor::as.tensor(tnsr), ranks=ranks, sparsity=sparsity)
    ))
  }

  gene_by_factors <- tucker_decomp$U[[2]]
  rownames(gene_by_factors) <- gene_nm
  ctype_by_factors <- tucker_decomp$U[[3]]
  rownames(ctype_by_factors) <- ctype_nm
  donor_mat <- tucker_decomp$U[[1]]
  rownames(donor_mat) <- donor_nm

  if (rotation_type=='hybrid') {
    ## rotate gene matrix
    # gene_by_factors <- ica::icafast(gene_by_factors,ranks[2],center=FALSE,
    #                                 alg='def',maxit = 200,tol = 1e-15)$S
    gene_by_factors <- ica::icafast(gene_by_factors,ranks[2],center=FALSE,alg='def')$S

    # make all vectors length 1 as ICA didn't preserve this
    all_rss <- c()
    for (j in 1:ncol(gene_by_factors)) {
      rss <- sqrt(sum(gene_by_factors[,j]**2))
      all_rss <- c(all_rss,rss)
    }
    gene_by_factors <- sweep(gene_by_factors,2,all_rss,FUN='/')

    # setting ctype factor matrix to identity
    ctype_by_factors <- diag(ncol(tucker_decomp$U[[3]]))
    rownames(ctype_by_factors) <- ctype_nm

    # compute kronecker product
    kron_prod <- kronecker(ctype_by_factors,gene_by_factors,make.dimnames = TRUE)

    # generate counter-rotated core tensor unfolded along donor dimension
    core_new <- t(as.matrix(donor_mat)) %*% rTensor::k_unfold(rTensor::as.tensor(tnsr),1)@data %*% kron_prod

    ## optimize core by rotating it with varimax
    vari_res <- stats::varimax(t(core_new),eps = 1e-15)
    # vari_res <- stats::varimax(t(core_new))
    core_new <- t(t(core_new) %*% vari_res$rotmat)
    ldngs <- core_new %*% t(kron_prod)

    # counter-rotate donor scores matrix
    donor_mat <- donor_mat %*% solve(t(vari_res$rotmat))

  } else if (rotation_type=='ica_dsc') {
    if (ranks[1]>1) {
      ## rotate donor scores matrix by ICA
      # donor_mat <- ica::icafast(donor_mat,ranks[1],center=FALSE,alg='def',
      #                           maxit = 2000,tol = 1e-20)$S
      donor_mat <- ica::icafast(donor_mat,ranks[1],center=FALSE,alg='def')$S

      # make all vectors length 1 as ICA didn't preserve this
      all_rss <- c()
      for (j in 1:ncol(donor_mat)) {
        rss <- sqrt(sum(donor_mat[,j]**2))
        all_rss <- c(all_rss,rss)
      }
      donor_mat <- sweep(donor_mat,2,all_rss,FUN='/')
    }

    # compute kronecker product
    kron_prod <- kronecker(ctype_by_factors,gene_by_factors,make.dimnames = TRUE)

    # generate rotated core tensor unfolded along donor dimension
    core_new <- t(as.matrix(donor_mat)) %*% rTensor::k_unfold(rTensor::as.tensor(tnsr),1)@data %*% kron_prod

    # compute loadings matrix with rotated core tensor
    ldngs <- core_new %*% t(kron_prod)
  } else if (rotation_type=='ica_lds') {
    # compute kronecker product
    kron_prod <- kronecker(ctype_by_factors,gene_by_factors,make.dimnames = TRUE)

    # generate rotated core tensor unfolded along donor dimension
    core_new <- t(as.matrix(donor_mat)) %*% rTensor::k_unfold(rTensor::as.tensor(tnsr),1)@data %*% kron_prod

    # compute loadings matrix with rotated core tensor
    ldngs <- core_new %*% t(kron_prod)

    ## rotate loadings directly with ICA
    # ica_res <- ica::icafast(t(ldngs),ranks[1],center=FALSE,alg='def',
    #                         maxit = 200,tol = 1e-15)
    ica_res <- ica::icafast(t(ldngs),ranks[1],center=FALSE,alg='def')
    ldngs <- t(ica_res$S)

    # counter rotate the donor scores matrix
    donor_mat <- donor_mat %*% solve(ica_res$W)

    # renormalize donor scores vectors and expand factor matrices by these constants
    all_rss <- c()
    for (j in 1:ncol(donor_mat)) {
      rss <- sqrt(sum(donor_mat[,j]**2))
      all_rss <- c(all_rss,rss)
    }
    donor_mat <- sweep(donor_mat,2,all_rss,FUN='/')

    for (j in 1:nrow(ldngs)) {
      ldngs[j,] <- ldngs[j,] * all_rss[j]
    }
  }

  return(list(donor_mat,ldngs))
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



#' Get the donor scores and loadings matrix for a single-factor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The number corresponding to the factor to extract
#'
#' @return a list with the first element as the donor scores and the second element
#' as the corresponding loadings matrix
#' @export
#' 
#' @examples
#' test_container <- get_one_factor(test_container, factor_select=1)
get_one_factor <- function(container, factor_select) {
  ldngs <- container$tucker_results[[2]]

  # break down a factor from the loadings matrix
  genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
  ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})

  sr_col <- ldngs[factor_select,]

  tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)

  return(list(container$tucker_results[[1]][,factor_select,drop=FALSE],tmp_casted_num))
}


#' Computes singular-value decomposition on the unfolded tensor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ranks numeric The number of factors to extract. Unlike with the Tucker
#' decomposition, this should be a single number.
#'
#' @return container with results of the decomposition in container$tucker_results
#' @export
#' 
#' @examples
#' test_container <- pca_unfolded(test_container, 2)
pca_unfolded <- function(container, ranks) {
  # get tensor data
  tensor_data <- container$tensor_data

  # extract tensor and labels
  donor_nm <- tensor_data[[1]]
  gene_nm  <- tensor_data[[2]]
  ctype_nm  <- tensor_data[[3]]
  tnsr <- tensor_data[[4]]

  d_unfold <- rTensor::k_unfold(rTensor::as.tensor(tnsr),1)@data

  rownames(d_unfold) <- donor_nm
  var_names <- sapply(ctype_nm, function(x) {
    sapply(gene_nm, function(y) {
      paste0(x,':',y)
    })
  })
  colnames(d_unfold) <- var_names

  svd_res <- svd(d_unfold)

  donor_mat <- svd_res$u[,1:ranks]
  rownames(donor_mat) <- rownames(d_unfold)
  ldngs <- t(svd_res$v)[1:ranks,]
  colnames(ldngs) <- colnames(d_unfold)

  # incorporate d values into ldngs
  ldngs <- svd_res$d[1:ranks] * ldngs

  container$tucker_results <- list(donor_mat,ldngs)

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

#' Computes non-negative matrix factorization on the unfolded tensor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ranks numeric The number of factors to extract. Unlike with the Tucker
#' decomposition, this should be a single number.
#'
#' @return container with results of the decomposition in container$tucker_results
#' @export
#' 
#' @examples
#' test_container <- nmf_unfolded(test_container, 2)
nmf_unfolded <- function(container, ranks) {
  # get tensor data
  tensor_data <- container$tensor_data

  # extract tensor and labels
  donor_nm <- tensor_data[[1]]
  gene_nm  <- tensor_data[[2]]
  ctype_nm  <- tensor_data[[3]]
  tnsr <- tensor_data[[4]]

  d_unfold <- rTensor::k_unfold(rTensor::as.tensor(tnsr),1)@data

  rownames(d_unfold) <- donor_nm
  var_names <- sapply(ctype_nm, function(x) {
    sapply(gene_nm, function(y) {
      paste0(x,':',y)
    })
  })
  colnames(d_unfold) <- var_names

  # make data non-negative by adding the minimum value of each gene
  col_m <- matrixStats::colMins(d_unfold)
  d_unfold <- sweep(d_unfold,MARGIN=2,col_m,FUN='-')

  # remove columns that are all 0
  ndx_keep <- which(colSums(d_unfold)!=0)
  d_unfold <- d_unfold[,ndx_keep]

  nmf_res <- NMF::nmf(d_unfold,ranks)
  donor_mat <- nmf_res@fit@W
  ldngs <- nmf_res@fit@H

  # incorporate magnitude values into ldngs
  all_rss <- c()
  for (j in 1:ncol(donor_mat)) {
    rss <- sqrt(sum(donor_mat[,j]**2))
    all_rss <- c(all_rss,rss)
  }
  donor_mat <- sweep(donor_mat,2,all_rss,FUN='/')
  ldngs <- t(sweep(t(ldngs),2,all_rss,FUN='*'))

  container$tucker_results <- list(donor_mat,ldngs)

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






