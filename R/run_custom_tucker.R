
#'Tucker Decomposition adapted from rTensor but with sparsity constraints
#'@import rTensor
#'@importFrom sparsepca rspca
#'
#'@param tnsr Tensor with K modes.
#'@param ranks numeric Vector of the modes of the output core Tensor
#'@param max_iter numeric Maximum number of iterations if error stays above tol
#'@param tol numeric Relative Frobenius norm error tolerance
#'
#'@return a list containing the results
#'@export
tucker_sparse <- function(tnsr,ranks=NULL,max_iter=25,tol=1e-5){
  if(is.null(ranks)) stop("ranks must be specified")
  if (sum(ranks>tnsr@modes)!=0) stop("ranks must be smaller than the corresponding mode")
  if (sum(ranks<=0)!=0) stop("ranks must be positive")

  #initialization via truncated hosvd
  num_modes <- tnsr@num_modes
  U_list <- vector("list",num_modes)
  for(m in 1:num_modes){
    temp_mat <- rs_unfold(tnsr,m=m)@data
    # using sparse truncated svd
    U_list[[m]] <- rspca(t(temp_mat),ranks[m],center=FALSE,alpha=.0001)$loadings
    # U_list[[m]] <- svd(temp_mat,nu=ranks[m])$u
  }
  tnsr_norm <- fnorm(tnsr)
  curr_iter <- 1
  converged <- FALSE
  #set up convergence check
  fnorm_resid <- rep(0, max_iter)
  CHECK_CONV <- function(Z,U_list){
    est <- ttl(Z,U_list,ms=1:num_modes)
    curr_resid <- fnorm(tnsr - est)
    fnorm_resid[curr_iter] <<- curr_resid
    if (curr_iter==1) return(FALSE)
    if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
    else{return(FALSE)}
  }

  #main loop (until convergence or max_iter)
  while((curr_iter < max_iter) && (!converged)){
    modes <- tnsr@modes
    modes_seq <- 1:num_modes
    for(m in modes_seq){
      #core Z minus mode m
      X <- ttl(tnsr,lapply(U_list[-m],t),ms=modes_seq[-m])
      #sparse truncated SVD of X
      U_list[[m]] <- rspca(t(rs_unfold(X,m=m)@data),ranks[m],center=FALSE,alpha=.0001)$loadings
      # U_list[[m]] <- svd(rs_unfold(X,m=m)@data,nu=ranks[m])$u
    }
    #compute core tensor Z
    Z <- ttm(X,mat=t(U_list[[num_modes]]),m=num_modes)

    #checks convergence
    if(CHECK_CONV(Z, U_list)){
      converged <- TRUE
      
      # ## trying out adding sparse pca as a last iteration...
      # for(m in modes_seq){
      #   #core Z minus mode m
      #   X <- ttl(tnsr,lapply(U_list[-m],t),ms=modes_seq[-m])
      #   #sparse truncated SVD of X
      #   U_list[[m]] <- rspca(t(rs_unfold(X,m=m)@data),ranks[m],center=FALSE,alpha=.0001)$loadings
      # }
      # #compute core tensor Z
      # Z <- ttm(X,mat=t(U_list[[num_modes]]),m=num_modes)
      
    }else{
      curr_iter <- curr_iter + 1
    }
  }
  #end of main loop
  #put together return list, and returns
  fnorm_resid <- fnorm_resid[fnorm_resid!=0]
  norm_percent<-(1-(tail(fnorm_resid,1)/tnsr_norm))*100
  est <- ttl(Z,U_list,ms=1:num_modes)
  invisible(list(Z=Z, U=U_list, conv=converged, est=est, norm_percent = norm_percent, fnorm_resid=tail(fnorm_resid,1), all_resids=fnorm_resid))
}


