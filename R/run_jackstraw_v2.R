
utils::globalVariables(c("donor_rank", "min_sig"))

#' Run tensor-based jackstraw to get gene_cell type combinations that are significantly
#' associated with donor scores for factors extracted by Tucker decomposition
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ranks numeric The number of donor, gene, and cell type ranks, respectively,
#' to decompose to using Tucker decomposition.
#' @param n_fibers numeric The number of fibers the randomly shuffle in each iteration
#' (default=100)
#' @param n_iter numeric The number of shuffling iterations to complete (default=500)
#' @param tucker_type character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
#' with sparsity constraints (default='regular')
#' @param rotation_type character Set to 'ica' to perform ICA rotation on resulting donor factor
#' matrix and loadings. Otherwise set to 'varimax' to perform varimax rotation. (default='ica')
#'
#' @return the project container with adjusted pvalues in container$gene_score_associations
#' @export
#'
#' @examples
#' test_container <- run_jackstraw(test_container, ranks=c(2,4,2), n_fibers=6,
#' n_iter=50, tucker_type='regular', rotation_type='ica')
run_jackstraw_v2 <- function(container, ranks, n_fibers=100, n_iter=500,
                          tucker_type='regular', rotation_type='ica') {
  # set random seed
  RNGkind("L'Ecuyer-CMRG")
  set.seed(container$experiment_params$rand_seed)
  
  # extract needed inputs from experiment parameters
  ncores <- container$experiment_params$ncores
  
  fstats_shuffled <- mclapply(1:n_iter, function(x) {
    # extract tensor data as we dont want to overwrite container
    tensor_data <- container$tensor_data
    
    # sample fibers and shuffle them across donors
    s_fibers <- sample_fibers_v2(tensor_data, n_fibers)
    tensor_data <- shuffle_fibers_v2(tensor_data, s_fibers)
    
    # compute tucker
    tucker_results <- tucker_ica_helper(tensor_data, ranks, tucker_type,
                                        rotation_type)
    
    
    # compute fiber-factor association F statistics for sampled fibers
    fiber_fstats <- calculate_fiber_fstats_v2(tensor_data, tucker_results, s_fibers)
    return(fiber_fstats)
  }, mc.cores = ncores)
  
  fstats_shuffled_all <- fstats_shuffled[[1]]
  for (i in 2:length(fstats_shuffled)) { # looping through iterations
    # need to loop through factors
    for (j in 1:length(fstats_shuffled_all)) {
      fstats_shuffled_all[[j]] <- c(fstats_shuffled_all[[j]],fstats_shuffled[[i]][[j]])
      
    }
  }
  
  # compute actual F statistics for all real fibers
  fstats_real <- get_real_fstats_v2(container, ncores)
  
  # # temporary testing getting regular pvalues
  # pvals_adj <- fstats_real
  # names(pvals_adj) <- names(fstats_real)
  # names(pvals_adj) <- sapply(names(pvals_adj),function(x){
  #   substr(x,1,nchar(x)-6)
  # })
  
  # calculate p-value by counting how many null F stats greater
  pvals_adj <- get_fstats_pvals_v2(fstats_real, fstats_shuffled_all)
  names(pvals_adj) <- names(fstats_real)
  container$gene_score_associations <- pvals_adj
  return(container)
}


#' Get list of tensor fibers to shuffle
#'
#' @param tensor_data list The tensor data including donor, gene, and cell type labels
#' as well as the tensor array itself
#' @param n_fibers numeric The number of fibers to get
#'
#' @return a list of gene and cell type indices for the randomly selected fibers
#' @export
sample_fibers_v2 <- function(tensor_data, n_fibers) {
  n_genes <- length(tensor_data[[2]])
  n_ctypes <- length(tensor_data[[3]])
  sample_space <- n_genes * n_ctypes
  samp_fibers <- sample(1:sample_space, n_fibers)
  fiber_coords <- lapply(samp_fibers, function(x) {
    ctype_ndx <- ceiling(x / n_genes)
    gene_ndx <- x %% n_genes
    if (gene_ndx == 0) {
      gene_ndx <- n_genes
    }
    return(list(gene_ndx,ctype_ndx))
  })
  return(fiber_coords)
}

#' Shuffle elements within selected fibers
#'
#' @param tensor_data list The tensor data including donor, gene, and cell type labels
#' as well as the tensor array itself
#' @param s_fibers list Gene and cell type indices for the randomly selected fibers
#'
#' @return the tensor_data with the values for the selected fibers shuffled
#' @export
shuffle_fibers_v2 <- function(tensor_data, s_fibers) {
  for (f in s_fibers) {
    gene_ndx <- f[[1]]
    ctype_ndx <- f[[2]]
    shuffled_fiber <- sample(tensor_data[[4]][,gene_ndx,ctype_ndx])
    tensor_data[[4]][,gene_ndx,ctype_ndx] <- shuffled_fiber
  }
  return(tensor_data)
}

#' Calculate F-Statistics for the association between donor scores for each factor
#' donor values of shuffled gene_ctype fibers
#'
#' @param tensor_data list The tensor data including donor, gene, and cell type labels
#' as well as the tensor array itself
#' @param tucker_results list The results from Tucker decomposition. Includes a scores
#' matrix as the first element and the loadings tensor unfolded as the second element.
#' @param s_fibers list Gene and cell type indices for the randomly selected fibers
#'
#' @return the F-Statistics for associations between all shuffled fibers and donor scores
#' @export
calculate_fiber_fstats_v2 <- function(tensor_data, tucker_results, s_fibers) {
  
  all_fstats <- lapply(1:ncol(tucker_results[[1]]), function(i) {
    fiber_fstats <- c()
    t=0
    for (f in s_fibers) { # loop through shuffled fibers
      t=t+1
      # print(t)
      gene_ndx <- f[[1]]
      ctype_ndx <- f[[2]]
      shuffled_fiber <- tensor_data[[4]][,gene_ndx,ctype_ndx]
      
      # calculate fstat for fiber-factor association
      df_test <- as.data.frame(cbind(shuffled_fiber, tucker_results[[1]][,i]))
      colnames(df_test) <- c('shuffled','factor')
      lmres <- lm(factor~shuffled,df_test)
      fstat <- summary(lmres)$fstatistic[[1]]
      fiber_fstats <- c(fiber_fstats,fstat)
    }
    return(fiber_fstats)
  })
  
  return(all_fstats)
}

#' Get F-Statistics for the real (non-shuffled) gene_ctype fibers
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ncores numeric The number of cores to use
#'
#' @return a vector F-Statistics for each gene-celltype-factor combination
#' @export
get_real_fstats_v2 <- function(container, ncores) {
  tensor_data <- container$tensor_data
  tucker_results <- container$tucker_results
  
  n_genes <- length(tensor_data[[2]])
  n_ctypes <- length(tensor_data[[3]])
  fstats_real <- data.frame(matrix(ncol=3,nrow=0))
  
  fstats_real <- mclapply(1:n_genes, function(i) {
    gene <- tensor_data[[2]][i]
    gene_res <- list()
    for (j in 1:n_ctypes) {
      ctype <- tensor_data[[3]][j]
      gene_res[[ctype]] <- list()
      for (k in 1:ncol(tucker_results[[1]])) {
        tmp_fiber <- tensor_data[[4]][,i,j]
        
        # if expression is 0 for all donors make the fstat na
        if (sum(tmp_fiber==0)==length(tmp_fiber)) {
          gene_res[[ctype]][[as.character(k)]] <- NA
        } else {
          df_test <- as.data.frame(cbind(tmp_fiber, tucker_results[[1]][,k]))
          colnames(df_test) <- c('fiber','factor')
          lmres <- lm(factor~fiber,df_test)
          fstat <- summary(lmres)$fstatistic[[1]]
          gene_res[[ctype]][[as.character(k)]] <- fstat
        }
      }
    }
    return(gene_res)
  }, mc.cores = ncores)
  
  names(fstats_real) <- tensor_data[[2]]
  
  # unpack the list
  fstats_real <- unlist(fstats_real)
  
  return(fstats_real)
}

#' Calculate adjusted pvalues for gene celltype fiber-donor score associations
#'
#' @param fstats_real numeric A vector of F-Statistics for gene-cell type-factor combinations
#' @param fstats_shuffled numeric A vector of null F-Statistics
#'
#' @return adjusted pvalues for associations of the unshuffled fibers with factor donor scores
#' @export
get_fstats_pvals_v2 <- function(fstats_real, fstats_shuffled_all) {
  raw_pvals <- c()
  for (i in 1:length(fstats_real)) {
    tmp <- strsplit(names(fstats_real)[i],split = '.', fixed = TRUE)[[1]]
    factor_cur <- as.numeric(tmp[[length(tmp)]])
    num_null_greater <- sum(fstats_shuffled_all[[factor_cur]] > fstats_real[i])
    pval <- num_null_greater / length(fstats_shuffled_all[[factor_cur]])
    raw_pvals <- c(raw_pvals, pval)
  }
  
  adj_pvals <- p.adjust(raw_pvals, method='fdr')
  
  # set the na pvals to 1
  adj_pvals[is.na(adj_pvals)] <- 1
  
  # get smallest non-zero adj pvalue to check if did enough iterations
  adj_pvals_no_zero <- adj_pvals[adj_pvals > 0]
  min_nonzero_padj <- min(adj_pvals_no_zero)
  if (min_nonzero_padj > 0.05) {
    print('Warning: smallest non-zero pvalue > 0.05. Do not trust zero pvalues.')
  }
  
  return(adj_pvals)
}


