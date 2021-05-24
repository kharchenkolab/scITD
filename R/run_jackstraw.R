
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
run_jackstraw <- function(container, ranks, n_fibers=100, n_iter=500,
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
    s_fibers <- sample_fibers(tensor_data, n_fibers)
    tensor_data <- shuffle_fibers(tensor_data, s_fibers)

    # compute tucker
    tucker_results <- tucker_ica_helper(tensor_data, ranks, tucker_type,
                                        rotation_type)


    # compute fiber-factor association F statistics for sampled fibers
    fiber_fstats <- calculate_fiber_fstats(tensor_data, tucker_results, s_fibers)
    return(fiber_fstats)
  }, mc.cores = ncores)

  fstats_shuffled <- unlist(fstats_shuffled)

  # compute actual F statistics for all real fibers
  fstats_real <- get_real_fstats(container, ncores)

  # # temporary testing getting regular pvalues
  # pvals_adj <- fstats_real
  # names(pvals_adj) <- names(fstats_real)
  # names(pvals_adj) <- sapply(names(pvals_adj),function(x){
  #   substr(x,1,nchar(x)-6)
  # })

  # calculate p-value by counting how many null F stats greater
  pvals_adj <- get_fstats_pvals(fstats_real, fstats_shuffled)
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
sample_fibers <- function(tensor_data, n_fibers) {
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
shuffle_fibers <- function(tensor_data, s_fibers) {
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
calculate_fiber_fstats <- function(tensor_data, tucker_results, s_fibers) {
  all_fstats <- c()
  for (f in s_fibers) {
    gene_ndx <- f[[1]]
    ctype_ndx <- f[[2]]
    shuffled_fiber <- tensor_data[[4]][,gene_ndx,ctype_ndx]

    # get F stats for each factor
    for (i in 1:ncol(tucker_results[[1]])) {
      df_test <- as.data.frame(cbind(shuffled_fiber, tucker_results[[1]][,i]))
      colnames(df_test) <- c('shuffled','factor')
      lmres <- lm(factor~shuffled,df_test)
      fstat <- summary(lmres)$fstatistic[[1]]
      all_fstats <- c(all_fstats,fstat)
    }
  }
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
get_real_fstats <- function(container, ncores) {
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
        df_test <- as.data.frame(cbind(tmp_fiber, tucker_results[[1]][,k]))
        colnames(df_test) <- c('fiber','factor')
        lmres <- lm(factor~fiber,df_test)
        fstat <- summary(lmres)$fstatistic[[1]]

        # # testing getting pval instead of fstat
        # x <- summary(lmres)
        # fstat <- stats::pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)

        gene_res[[ctype]][[as.character(k)]] <- fstat
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
get_fstats_pvals <- function(fstats_real, fstats_shuffled) {
  raw_pvals <- c()
  for (i in 1:length(fstats_real)) {
    num_null_greater <- sum(fstats_shuffled > fstats_real[i])
    pval <- num_null_greater / length(fstats_shuffled)
    raw_pvals <- c(raw_pvals, pval)
  }

  adj_pvals <- p.adjust(raw_pvals, method='fdr')

  # get smallest non-zero adj pvalue to check if did enough iterations
  adj_pvals_no_zero <- adj_pvals[adj_pvals > 0]
  min_nonzero_padj <- min(adj_pvals_no_zero)
  if (min_nonzero_padj > 0.05) {
    print('Warning: smallest non-zero pvalue > 0.05. Do not trust zero pvalues.')
  }

  return(adj_pvals)
}



#' Evaluate the minimum number for significant genes in any factor for a given number of
#' factors extracted by tucker.
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses. Should have
#' @param donor_rank_range numeric Range of possible number of donor factors to use.
#' @param gene_ranks numeric The number of gene ranks to use in the decomposition
#' @param ctype_ranks numeric The number of cell type ranks to use in the decomposition
#' @param tucker_type character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
#' with sparsity constraints (default='regular')
#' @param rotation_type character Set to 'ica' to perform ICA rotation on resulting donor factor
#' matrix and loadings. Otherwise set to 'varimax' to perform varimax rotation. (default='ica')
#' @param n_fibers numeric The number of fibers the randomly shuffle in each iteration
#' (default=100)
#' @param n_iter numeric The number of shuffling iterations to complete (default=500)
#' @param thresh numeric Pvalue threshold for significant genes in calculating the
#' number of significant genes identified per factor. (default=0.05)
#'
#' @return the project container with a plot of the minimum significant genes for
#' each decomposition with varying number of donor factors located in
#' container$plots$min_sig_genes
#' @export
get_min_sig_genes <- function(container, donor_rank_range, gene_ranks,
                              ctype_ranks, tucker_type='regular',
                              rotation_type='ica', n_fibers=100, n_iter=500,
                              thresh=0.05) {

  min_per_decomp <- data.frame(matrix(ncol=2,nrow=0))
  colnames(min_per_decomp) <- c('donor_rank','min_sig')
  for (i in donor_rank_range) {
    container <- run_tucker_ica(container, c(i,gene_ranks,ctype_ranks),
                                tucker_type=tucker_type, rotation_type=rotation_type)
    container <- run_jackstraw(container, ranks=c(i,gene_ranks,ctype_ranks),
                               n_fibers=n_fibers, n_iter=n_iter,
                               tucker_type=tucker_type, rotation_type=rotation_type)

    padj <- container$gene_score_associations
    padj_factors <- sapply(names(padj),function(x) {
      strsplit(x,split = '.', fixed = TRUE)[[1]][[3]]
    })

    # loop through factors to get the min number of significant genes out of any factor
    num_sig_genes <- c()
    for (j in 1:i) {
      padj_use <- padj[which(padj_factors == as.character(j))]
      num_sig_genes <- c(num_sig_genes, sum(padj_use < thresh))
    }
    print(num_sig_genes)
    tmp <- as.data.frame(t(c(i,min(num_sig_genes))))
    colnames(tmp) <- colnames(min_per_decomp)
    min_per_decomp <- rbind(min_per_decomp,tmp)
  }

  # plot results
  p <- ggplot(min_per_decomp, aes(x=donor_rank,y=min_sig)) +
    geom_line() +
    xlab("Number of Donor Factors") +
    ylab("Minimum Significant Genes of Any Factor")

  container$plots$min_sig_genes <- p

  return(container)
}

























