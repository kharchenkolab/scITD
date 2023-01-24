
utils::globalVariables(c("donor_rank", "min_sig"))

#' Run jackstraw to get genes that are significantly associated with donor scores
#' for factors extracted by Tucker decomposition
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ranks numeric The number of donor ranks and gene ranks to decompose to
#' using Tucker decomposition
#' @param n_fibers numeric The number of fibers the randomly shuffle in each iteration
#' (default=100)
#' @param n_iter numeric The number of shuffling iterations to complete (default=500)
#' @param tucker_type character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
#' with sparsity constraints (default='regular')
#' @param rotation_type character Set to 'hybrid' to perform hybrid rotation on resulting donor factor
#' matrix and loadings. Otherwise set to 'ica_lds' to perform ica rotation on loadings or
#' ica_dsc to perform ica on donor scores. (default='hybrid')
#' @param seed numeric Seed passed to set.seed() (default=container$experiment_params$rand_seed)
#' @param ncores numeric The number of cores to use (default=container$experiment_params$ncores)
#'
#' @return The project container with a vector of adjusted pvalues in container$gene_score_associations.
#' @export
#' 
#' @examples
#' test_container <- run_jackstraw(test_container, ranks=c(2,4), n_fibers=2, n_iter=10,
#' tucker_type='regular', rotation_type='hybrid', ncores=1)
run_jackstraw <- function(container, ranks, n_fibers=100, n_iter=500,
                          tucker_type='regular', rotation_type='hybrid', 
                          seed=container$experiment_params$rand_seed, ncores=container$experiment_params$ncores) {
  # set random seed
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)


  fstats_shuffled <- sccore::plapply(1:n_iter, function(x) {
    # extract tensor data as we dont want to overwrite container
    tensor_data <- container$tensor_data

    # set third element of ranks to the number of cell types
    ranks[3] <- length(tensor_data[[3]])

    # sample fibers and shuffle them across donors
    s_fibers <- sample_fibers(tensor_data, n_fibers)
    tensor_data <- shuffle_fibers(tensor_data, s_fibers)

    # compute tucker
    tucker_results <- tucker_ica_helper(tensor_data, ranks, tucker_type,
                                        rotation_type)


    # compute fiber-factor association F statistics for sampled fibers
    fiber_fstats <- calculate_fiber_fstats(tensor_data, tucker_results, s_fibers)
    return(fiber_fstats)
  }, n.cores = ncores, mc.preschedule=TRUE, progress = TRUE)

  fstats_shuffled <- unlist(fstats_shuffled)

  # compute actual F statistics for all real fibers
  fstats_real <- get_real_fstats(container, ncores)

  # calculate p-value by counting how many null F stats greater
  pvals_adj <- get_fstats_pvals(fstats_real, fstats_shuffled)
  names(pvals_adj) <- names(fstats_real)
  container$gene_score_associations <- pvals_adj
  return(container)
}


#' Get a list of tensor fibers to shuffle
#'
#' @param tensor_data list The tensor data including donor, gene, and cell type labels
#' as well as the tensor array itself
#' @param n_fibers numeric The number of fibers to get
#'
#' @return A list of gene and cell type indices for the randomly selected fibers
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

#' Shuffle elements within the selected fibers
#'
#' @param tensor_data list The tensor data including donor, gene, and cell type labels
#' as well as the tensor array itself
#' @param s_fibers list Gene and cell type indices for the randomly selected fibers
#'
#' @return The tensor_data object with the values for the selected fibers shuffled.
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
#' @return A numeric vector of F-statistics for associations between all shuffled fibers 
#' and donor scores.
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
#' @return A vector F-statistics for each gene_celltype-factor association of the
#' unshuffled data.
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
        if (all(tmp_fiber==0)) {
          fstat <- 0
        } else {
          df_test <- as.data.frame(cbind(tmp_fiber, tucker_results[[1]][,k]))
          colnames(df_test) <- c('fiber','factor')
          lmres <- lm(factor~fiber,df_test)
          fstat <- summary(lmres)$fstatistic[[1]]
        }

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

#' Calculate adjusted p-values for gene_celltype fiber-donor score associations
#'
#' @param fstats_real numeric A vector of F-Statistics for gene-cell type-factor combinations
#' @param fstats_shuffled numeric A vector of null F-Statistics
#'
#' @return A vector of adjusted p-values for associations of the unshuffled fibers with
#' factor donor scores.
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
    warning('Warning: smallest non-zero pvalue > 0.05. Do not trust zero pvalues.')
  }

  return(adj_pvals)
}



#' Evaluate the minimum number for significant genes in any factor for a given number of
#' factors extracted by the decomposition
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses. Should have
#' @param donor_rank_range numeric Range of possible number of donor factors to use.
#' @param gene_ranks numeric The number of gene ranks to use in the decomposition
#' @param use_lm logical Set to true to use get_lm_pvals otherwise uses jackstraw (default=TRUE)
#' @param tucker_type character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
#' with sparsity constraints (default='regular')
#' @param rotation_type character Set to 'hybrid' to perform hybrid rotation on resulting donor factor
#' matrix and loadings. Otherwise set to 'ica_lds' to perform ica rotation on loadings or
#' ica_dsc to perform ica on donor scores. (default='hybrid')
#' @param n_fibers numeric The number of fibers the randomly shuffle in each jackstraw iteration
#' (default=100)
#' @param n_iter numeric The number of jackstraw shuffling iterations to complete (default=500)
#' @param n.cores Number of cores to use in get_lm_pvals() (default = container$experiment_params$ncores)
#'
#' @param thresh numeric Pvalue threshold for significant genes in calculating the
#' number of significant genes identified per factor. (default=0.05)
#'
#' @return The project container with a plot of the minimum significant genes for
#' each decomposition with varying number of donor factors located in
#' container$plots$min_sig_genes.
#' @export
#' 
#' @examples
#' test_container <- get_min_sig_genes(test_container, donor_rank_range=c(2:4),
#' gene_ranks=4, tucker_type='regular', rotation_type='hybrid', n.cores=1)
get_min_sig_genes <- function(container, donor_rank_range, gene_ranks,
                              use_lm=TRUE, tucker_type='regular',
                              rotation_type='hybrid', n_fibers=100, n_iter=500,
                              n.cores = container$experiment_params$ncores, thresh=0.05) {

  min_per_decomp <- data.frame(matrix(ncol=2,nrow=0))
  colnames(min_per_decomp) <- c('donor_rank','min_sig')
  for (i in donor_rank_range) {
    container <- run_tucker_ica(container, c(i,gene_ranks),
                                tucker_type=tucker_type, rotation_type=rotation_type)
    if (use_lm) {
      container <- get_lm_pvals(container, n.cores=n.cores)
    } else {
      container <- run_jackstraw(container, ranks=c(i,gene_ranks),
                                 n_fibers=n_fibers, n_iter=n_iter,
                                 tucker_type=tucker_type, rotation_type=rotation_type,
                                 ncores=n.cores)
    }

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


#' Compute gene-factor associations using univariate linear models
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param n.cores Number of cores to use (default = container$experiment_params$ncores)
#'
#' @return The project container with a vector of adjusted p-values for the gene-factor
#' associations in container$gene_score_associations.
#' @export
#'
#' @examples
#' test_container <- get_lm_pvals(test_container, n.cores=1)
get_lm_pvals <- function(container, n.cores = container$experiment_params$ncores) {
  tensor_data <- container$tensor_data
  tucker_results <- container$tucker_results

  if (is.null(tucker_results)) {
    stop('Need to run run_tucker_ica() first')
  }

  n_genes <- length(tensor_data[[2]])
  n_ctypes <- length(tensor_data[[3]])
  all_pvals <- data.frame(matrix(ncol=3,nrow=0))
  if (is.null(container$experiment_params$ncores)){
    n.cores = 1
  }
  all_pvals <- mclapply(1:n_genes, function(i) {
    gene <- tensor_data[[2]][i]
    gene_res <- list()
    for (j in 1:n_ctypes) {
      ctype <- tensor_data[[3]][j]
      gene_res[[ctype]] <- list()
      for (k in 1:ncol(tucker_results[[1]])) {
        tmp_fiber <- tensor_data[[4]][,i,j]

        # if expression is 0 for all donors just skip
        if (sum(tmp_fiber==0)==length(tmp_fiber)) {
          gene_res[[ctype]][[as.character(k)]]$value <- NA
        } else {
          df_test <- as.data.frame(cbind(tmp_fiber, tucker_results[[1]][,k]))
          colnames(df_test) <- c('fiber','factor')
          lmres <- lm(factor~fiber,df_test)

          x <- summary(lmres)
          pval <- stats::pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)

          gene_res[[ctype]][[as.character(k)]] <- pval
        }
      }
    }
    return(gene_res)
  }, mc.cores = n.cores)

  names(all_pvals) <- tensor_data[[2]]

  # unpack the list
  all_pvals <- unlist(all_pvals)

  all_pvals <- p.adjust(all_pvals,method='fdr')

  # set NA values to 1
  all_pvals[is.na(all_pvals)] <- 1

  new_names <- sapply(names(all_pvals),function(x) {
    tmp <- strsplit(x,split = '.', fixed = TRUE)[[1]]
    if (length(tmp)==4) {
      return(paste0(tmp[[1]],'.',tmp[[2]],'.',tmp[[3]]))
    } else if (length(tmp)==5) {
      return(paste0(tmp[[1]],'.',tmp[[2]],'.',tmp[[3]],'.',tmp[[4]]))
    } else if (length(tmp)==6) {
      return(paste0(tmp[[1]],'.',tmp[[2]],'.',tmp[[3]],'.',tmp[[4]],'.',tmp[[5]]))
    }
  })
  names(new_names) <- NULL
  names(all_pvals) <- new_names
  container[["gene_score_associations"]] <- all_pvals
  return(container)
}





#' Get significant genes for a factor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param factor_select numeric The number corresponding to the factor to extract
#'
#' @return A gene by cell type matrix of gene significance p-values for a factor
#' @export
get_one_factor_gene_pvals <- function(container, factor_select) {
  if (is.null(container$gene_score_associations)) {
    stop('Run get_lm_pvals() first to compute the gene p-values')
  }
  
  ldngs <- container$tucker_results[[2]]
  
  # break down a factor from the loadings matrix
  genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
  ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})
  
  sr_col <- ldngs[factor_select,]
  
  tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)
  
  sig_vectors <- get_significance_vectors(container,
                                          factor_select, colnames(tmp_casted_num))
  # convert list to df
  sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
  
  # order df same way as in tmp_casted_num
  sig_df <- sig_df[rownames(tmp_casted_num),colnames(tmp_casted_num)]
  
  return(sig_df)
}











