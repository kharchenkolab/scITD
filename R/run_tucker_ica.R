
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
get_one_factor <- function(container, factor_select) {
  ldngs <- container$tucker_results[[2]]

  # break down a factor from the loadings matrix
  genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
  ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})

  sr_col <- ldngs[factor_select,]

  tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)

  return(list(container$tucker_results[[1]][,factor_select,drop=FALSE],tmp_casted_num))
}




icafast2 <-
  function(X, nc, center = TRUE, maxit = 100, tol = 1e-6,
           Rmat = diag(nc), alg = c("par", "def"),
           fun = c("logcosh", "exp", "kur"), alpha = 1){
    ###### Fast Independent Component Analysis
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: May 23, 2018

    # get sdiag() from ica-internal
    sdiag <- utils::getFromNamespace("sdiag", "ica")


    ### initial checks
    X <- as.matrix(X)
    nobs <- nrow(X)
    nvar <- ncol(X)
    nc <- as.integer(nc[1])
    if(nc<1){ stop("Must set nc>=1 component.") }
    maxit <- as.integer(maxit[1])
    if(maxit<1){ stop("Must set maxit>=1 iteration.") }
    tol <- tol[1]
    if(tol<=0){ stop("Must set ctol>0.") }
    if(nc>min(nobs,nvar)){ stop("Too many components. Set nc<=min(dim(X)).") }
    alpha <- alpha[1]
    if(alpha<1 | alpha>2){ stop("Must set 'alpha' between 1 and 2.") }
    if(nrow(Rmat)!=nc | ncol(Rmat)!=nc){ stop("Input 'Rmat' must be nc-by-nc rotation matrix.") }

    ### center and whiten
    xeig <- myeigen(crossprod(X)/nobs,symmetric=TRUE)
    nze <- sum(xeig$val>xeig$val[1]*.Machine$double.eps)
    if(nze<nc){
      warning("Numerical rank of X is less than requested number of components (nc).\n  Number of components has been redefined as the numerical rank of X.")
      nc <- nze
      Rmat <- diag(nc)
    }
    Dmat <- sdiag(sqrt(xeig$val[1:nc]))
    Mprt <- tcrossprod(Dmat,xeig$vec[,1:nc])
    diag(Dmat) <- 1/diag(Dmat)
    Pmat <- xeig$vec[,1:nc]%*%Dmat
    Xw <- X%*%Pmat   # whitened data

    ### check if nc=1
    if(nc==1L){
      return(list(S=Xw,M=Mprt,W=t(Pmat),Y=Xw,Q=t(Pmat),R=matrix(1),
                  vafs=(sum(Mprt^2)*nobs)/sum(X^2),iter=NA,
                  alg=alg,fun=fun,alpha=alpha))
    }

    ### define contrast function derivatives
    if(fun[1]=="kur"){
      fun1d <- function(x){ x^3 }
      fun2d <- function(x){ 3*(x^2) }
    } else if(fun[1]=="exp"){
      fun1d <- function(x){ x*exp(-(x^2)/2) }
      fun2d <- function(x){ exp(-(x^2)/2)*(1-x^2) }
    } else {
      fun1d <- function(x){ tanh(alpha*x) }
      fun2d <- function(x){ alpha*(1-tanh(alpha*x)^2) }
    }

    ### determine method
    if(alg[1]=="def"){
      myiters <- rep(NA,nc)
      for(j in 1:nc){
        if (j<2){
          # first component
          Rmat[,j] <- Rmat[,j]/sqrt(sum(Rmat[,j]^2))
          iter <- 0
          vtol <- 1
          while(vtol>tol && iter<maxit){
            # update first component
            svec <- Xw%*%Rmat[,j]
            rnew <- colMeans(Xw*matrix(fun1d(svec),nobs,nc))
            rnew <- rnew - mean(fun2d(svec))*Rmat[,j]
            rnew <- rnew / sqrt(sum(rnew^2))
            # check for convergence
            vtol <- 1 - abs(sum(Rmat[,j]*rnew))
            iter <- iter + 1
            Rmat[,j] <- rnew
          }
          myiters[j] <- iter
        } else {
          # decorrelate with previous components
          Rmat[,j] <- Rmat[,j] / sqrt(sum(Rmat[,j]^2))
          svec <- matrix(0,nc,1)
          for(k in 1:(j-1)){
            svec <- svec + sum(Rmat[,k]*Rmat[,j])*Rmat[,k]
          }
          Rmat[,j] <- Rmat[,j] - svec
          Rmat[,j] <- Rmat[,j] / sqrt(sum(Rmat[,j]^2))
          # get j-th component
          iter <- 0
          vtol <- 1
          while(vtol>tol && iter<maxit){
            # update j-th component
            svec <- Xw%*%Rmat[,j]
            rnew <- colMeans(Xw*matrix(fun1d(svec),nobs,nc))
            rnew <- rnew - mean(fun2d(svec))*Rmat[,j]
            rnew <- rnew/sqrt(sum(rnew^2))
            # decorrelate j-th component
            svec <- matrix(0,nc,1)
            for(k in 1:(j-1)){
              svec <- svec + sum(Rmat[,k]*rnew)*Rmat[,k]
            }
            rnew <- rnew - svec
            rnew <- rnew / sqrt(sum(rnew^2))
            # check for convergence
            vtol <- 1 - abs(sum(Rmat[,j]*rnew))
            iter <- iter + 1
            Rmat[,j] <- rnew
          }
          myiters[j] <- iter
        } # if (j<2)
      } # end for(j in 1:nc)
    } else {
      # parallel (symmetric) decorrelation
      rsvd <- svd(Rmat)
      Rmat <- tcrossprod(rsvd$u,rsvd$v)
      iter <- 0
      vtol <- 1
      while(vtol>tol && iter<maxit){
        # update all components
        smat <- Xw%*%Rmat
        rnew <- crossprod(Xw,fun1d(smat))/nobs
        rnew <- rnew - Rmat%*%sdiag(colMeans(fun2d(smat)))
        rsvd <- svd(rnew)
        rnew <- tcrossprod(rsvd$u,rsvd$v)
        # check for convergence
        vtol <- 1 - min(abs(colSums(Rmat*rnew)))
        iter <- iter + 1
        Rmat <- rnew
      }
      myiters <- iter
    } # end if(alg[1]=="def")

    ### sort according to vafs
    M <- crossprod(Rmat,Mprt)
    vafs <- rowSums(M^2)
    ix <- sort(vafs,decreasing=TRUE,index.return=TRUE)$ix
    M <- M[ix,]
    Rmat <- Rmat[,ix]
    vafs <- (vafs[ix]*nobs) / sum(X^2)

    return(list(S=Xw%*%Rmat,M=t(M),W=t(Pmat%*%Rmat),Y=Xw,
                Q=t(Pmat),R=Rmat,vafs=vafs,iter=myiters,
                alg=alg[1],fun=fun,alpha=alpha))

  }




ica_p1 <-
  function(X, nc, center = TRUE, maxit = 100, tol = 1e-6,
           Rmat = diag(nc), alg = c("par", "def"),
           fun = c("logcosh", "exp", "kur"), alpha = 1){
    ###### Fast Independent Component Analysis
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: May 23, 2018

    # get sdiag() from ica-internal
    sdiag <- utils::getFromNamespace("sdiag", "ica")


    ### initial checks
    X <- as.matrix(X)
    nobs <- nrow(X)
    nvar <- ncol(X)
    nc <- as.integer(nc[1])
    if(nc<1){ stop("Must set nc>=1 component.") }
    maxit <- as.integer(maxit[1])
    if(maxit<1){ stop("Must set maxit>=1 iteration.") }
    tol <- tol[1]
    if(tol<=0){ stop("Must set ctol>0.") }
    if(nc>min(nobs,nvar)){ stop("Too many components. Set nc<=min(dim(X)).") }
    alpha <- alpha[1]
    if(alpha<1 | alpha>2){ stop("Must set 'alpha' between 1 and 2.") }
    if(nrow(Rmat)!=nc | ncol(Rmat)!=nc){ stop("Input 'Rmat' must be nc-by-nc rotation matrix.") }

    ### center and whiten
    cp <- crossprod(X)
    # xeig <- eigen(crossprod(X)/nobs,symmetric=TRUE)
    # nze <- sum(xeig$val>xeig$val[1]*.Machine$double.eps)
    # if(nze<nc){
    #   warning("Numerical rank of X is less than requested number of components (nc).\n  Number of components has been redefined as the numerical rank of X.")
    #   nc <- nze
    #   Rmat <- diag(nc)
    # }
    # Dmat <- sdiag(sqrt(xeig$val[1:nc]))
    # Mprt <- tcrossprod(Dmat,xeig$vec[,1:nc])
    # diag(Dmat) <- 1/diag(Dmat)
    # Pmat <- xeig$vec[,1:nc]%*%Dmat
    # Xw <- X%*%Pmat   # whitened data


    return(cp)

  }

ica_p2 <-
  function(X, nc, center = TRUE, maxit = 100, tol = 1e-6,
           Rmat = diag(nc), alg = c("par", "def"),
           fun = c("logcosh", "exp", "kur"), alpha = 1){
    ###### Fast Independent Component Analysis
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: May 23, 2018

    # get sdiag() from ica-internal
    sdiag <- utils::getFromNamespace("sdiag", "ica")


    ### initial checks
    X <- as.matrix(X)
    nobs <- nrow(X)
    nvar <- ncol(X)
    nc <- as.integer(nc[1])
    if(nc<1){ stop("Must set nc>=1 component.") }
    maxit <- as.integer(maxit[1])
    if(maxit<1){ stop("Must set maxit>=1 iteration.") }
    tol <- tol[1]
    if(tol<=0){ stop("Must set ctol>0.") }
    if(nc>min(nobs,nvar)){ stop("Too many components. Set nc<=min(dim(X)).") }
    alpha <- alpha[1]
    if(alpha<1 | alpha>2){ stop("Must set 'alpha' between 1 and 2.") }
    if(nrow(Rmat)!=nc | ncol(Rmat)!=nc){ stop("Input 'Rmat' must be nc-by-nc rotation matrix.") }

    ### center and whiten
    cp <- crossprod(X)
    xeig <- eigen(cp/nobs,symmetric=TRUE)$vectors
    # nze <- sum(xeig$val>xeig$val[1]*.Machine$double.eps)
    # if(nze<nc){
    #   warning("Numerical rank of X is less than requested number of components (nc).\n  Number of components has been redefined as the numerical rank of X.")
    #   nc <- nze
    #   Rmat <- diag(nc)
    # }
    # Dmat <- sdiag(sqrt(xeig$val[1:nc]))
    # Mprt <- tcrossprod(Dmat,xeig$vec[,1:nc])
    # diag(Dmat) <- 1/diag(Dmat)
    # Pmat <- xeig$vec[,1:nc]%*%Dmat
    # Xw <- X%*%Pmat   # whitened data


    return(xeig)

  }



myeigen <- function(x, symmetric, only.values = FALSE, EISPACK = FALSE)
{
  x <- unname(as.matrix(x))
  n <- nrow(x)
  if (!n) stop("0 x 0 matrix")
  if (n != ncol(x)) stop("non-square matrix in 'eigen'")
  n <- as.integer(n)
  if(is.na(n)) stop("invalid nrow(x)")

  complex.x <- is.complex(x)
  if (!all(is.finite(x))) stop("infinite or missing values in 'x'")

  if(missing(symmetric)) symmetric <- isSymmetric.matrix(x)

  if (symmetric) {
    z <- if(!complex.x) .Internal(La_rs(x, only.values))
    else .Internal(La_rs_cmplx(x, only.values))
    ord <- rev(seq_along(z$values))
  } else {
    z <- if(!complex.x) .Internal(La_rg(x, only.values))
    else .Internal(La_rg_cmplx(x, only.values))
    ord <- sort.list(Mod(z$values), decreasing = TRUE)
  }
  return(list(values = z$values[ord],
              vectors = if (!only.values) z$vectors[, ord, drop = FALSE]))
}










