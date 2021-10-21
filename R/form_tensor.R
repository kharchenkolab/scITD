
#' Form the tensor and scale the data
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param donor_min_cells numeric Minimum threshold for number of cells per
#' donor (default=5)
#' @param norm_method character The normalization method to use on the pseudobulked
#' count data. Set to 'regular' to do standard normalization of dividing by
#' library size. Set to 'trim' to use edgeR trim-mean normalization, whereby counts
#' are divided by library size times a normalization factor. (default='trim')
#' @param scale_factor numeric The number that gets multiplied by fractional counts
#' during normalization of the pseudobulked data (default=10000)
#' @param vargenes_method character The method by which to select highly variable
#' genes from each cell type. Set to 'anova' or 'norm_var' (default='norm_var')
#' @param vargenes_thresh numeric The threshold to use in variable gene selection.
#' For 'anova' and 'empir' this should be a p-value threshold. For 'norm_var' this
#' should be the number of most variably expressed genes to select from each cell
#' type (default=500)
#' @param batch_var character A batch variable from metadata to remove (default=NULL)
#' @param scale_var logical TRUE to scale the gene expression variance across donors
#' for each cell type. If FALSE then all genes are scaled to unit variance across
#' donors for each cell type. (default=TRUE)
#' @param var_scale_power numeric Exponent of normalized variance that is
#' used for variance scaling. Variance for each gene
#' is initially set to unit variance across donors (for a given cell type).
#' Variance for each gene is then scaled by multiplying the unit scaled values
#' by each gene's normalized variance (where the effect of the mean-variance
#' dependence is taken into account) to the exponent specified here.
#' If NULL, uses var_scale_power from container$experiment_params. (default=.5)
#' @param custom_genes character A vector of genes to include in the tensor.
#' Overrides the default gene selection if not NULL. (default=NULL)
#' @param verbose logical Set to TRUE to print out progress (default=TRUE)
#'
#' @return the project container with tensor data added in the
#' container$tensor_data slot
#' @export
#'
#' @examples
#' test_container <- form_tensor(test_container, donor_min_cells=0,
#' norm_method='trim', scale_factor=10000, vargenes_method='norm_var', vargenes_thresh=500,
#' scale_var = TRUE, var_scale_power = 1.5)
form_tensor <- function(container, donor_min_cells=5, norm_method='trim',
                        scale_factor=10000, vargenes_method='norm_var',
                        vargenes_thresh=500, batch_var=NULL, scale_var=TRUE,
                        var_scale_power=.5, custom_genes=NULL, verbose=TRUE) {
  # parse data by cell type
  if (verbose) {
    print('parsing data matrix by cell/tissue type...')
  }
  container <- parse_data_by_ctypes(container)

  # clean counts
  if (verbose) {
    print('cleaning data...')
  }
  container <- clean_data(container, donor_min_cells=donor_min_cells)

  # collapse data to donor-level
  if (verbose) {
    print('collapsing count matrices from cells to donors (aka pseudobulk operation)...')
  }
  container <- get_pseudobulk(container)

  # normalize data
  if (verbose) {
    print('normalizing data...')
  }
  container <- normalize_pseudobulk(container, method=norm_method, scale_factor=scale_factor)

  # get normalized variances
  if (verbose) {
    print('calculating gene overdispersion factors...')
  }
  container <- get_normalized_variance(container)

  # reduce number of genes to use in the tensor
  if (!is.null(custom_genes)) {
    if (verbose) {
      print('reducing tensor to selected genes...')
    }
    # check custom genes all in tensor
    all_genes <- colnames(container$scMinimal_ctype[[1]]$pseudobulk)
    if (any(!(custom_genes %in% all_genes))) {
      stop('some of custom_genes are not present in the count data')
    }
    # set them as "vargenes" and reduce data to just these
    container$all_vargenes <- custom_genes
    container <- reduce_to_vargenes(container)
  } else {
    # select highly variable genes
    if (verbose) {
      print('selecting highly variable genes from each cell type...')
    }
    container <- get_ctype_vargenes(container, method=vargenes_method, thresh=vargenes_thresh)
  }

  if (scale_var) {
    # scale gene expression
    if (verbose) {
      print('scaling variance...')
    }
    container <- scale_variance(container,var_scale_power=var_scale_power)
  }

  # apply batch correction if specified
  if (!is.null(batch_var)) {
    if (verbose) {
      print('applying ComBat for batch correction...')
    }
    container <- apply_combat(container,batch_var=batch_var)
  }

  # build the tensor
  if (verbose) {
    print('forming tensor...')
  }
  container <- stack_tensor(container)

  if (verbose) {
    print('Complete!')
  }

  return(container)
}

#' Parse main counts matrix into per-celltype-matrices
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#'
#' @return the project container with scMinimal objects per cell type
#' @export
parse_data_by_ctypes <- function(container) {
  # check that ctypes_use param has been set
  if (is.null(container$experiment_params$ctypes_use)) {
    stop("ctypes_use parameter from container$experiment_params is NULL. Use set_experiment_params()")
  }

  for (ct in container$experiment_params$ctypes_use) {
    ctype_sub <- subset_scMinimal(container$scMinimal_full, ctypes_use=ct, in_place=FALSE)
    container$scMinimal_ctype[[ct]] <- ctype_sub
  }

  return(container)
}

#' Clean data to remove genes only expressed in a few cells and donors
#' with very few cells
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param donor_min_cells numeric Minimum threshold for number of cells per
#' donor (default=5)
#'
#' @return the project container with cleaned counts matrices
#' @export
clean_data <- function(container, donor_min_cells=5) {
  for (ct in container$experiment_params$ctypes_use) {
    ctype_sub <- container$scMinimal_ctype[[ct]]

    # identify donors with few cells
    donor_counts <- table(ctype_sub$metadata$donors)
    donors_keep <- names(donor_counts)[donor_counts > donor_min_cells]

    # subset on donors
    ctype_sub <- subset_scMinimal(ctype_sub, donors_use = donors_keep)
  }

  # get donors present in all ctype matrices
  donors_in_all <- unique(container$scMinimal_ctype[[1]]$metadata$donors)
  for (ct in container$experiment_params$ctypes_use) {
    ctype_donors <- unique(container$scMinimal_ctype[[ct]]$metadata$donors)
    donors_in_all <- intersect(donors_in_all,ctype_donors)
  }

  # reduce data to only the intersection of donors that have all ctypes
  for (ct in container$experiment_params$ctypes_use) {
    ctype_sub <- container$scMinimal_ctype[[ct]]
    ctype_sub <- subset_scMinimal(ctype_sub, donors_use=donors_in_all)
  }

  print(paste0('Keeping ',length(donors_in_all),' donors. All donors have at least ',donor_min_cells,' cells in each cell type included.'))

  # get total num donors
  total_num_donors <- length(unique(container$scMinimal_full$metadata$donors))

  if (length(donors_in_all) < (.5*total_num_donors)) {
    print('Consider using fewer cell types or reducing the donor_min_cells parameter to include more donors.')
  }

  return(container)
}


#' Collapse data from cell-level to donor-level via summing counts
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param shuffle logical Set to TRUE to shuffle cell-donor linkages (default=FALSE)
#' @param shuffle_within character A metadata variable to shuffle cell-donor linkages
#' within (default=NULL)
#'
#' @return the project container with pseudobulked count matrices in
#' scMinimal_ctype$ctype$pseudobulk slots for each cell type
#' @export
get_pseudobulk <- function(container,shuffle=FALSE,shuffle_within=NULL) {
  for (ct in container$experiment_params$ctypes_use) {
    scMinimal <- container$scMinimal_ctype[[ct]]

    cdata <- t(scMinimal$count_data)

    if (shuffle) {
      if (is.null(shuffle_within)) {
        donor_meta <- as.factor(sample(scMinimal$metadata$donors))
      } else {
        bvars <- unique(scMinimal$metadata[[shuffle_within]])
        donor_meta <- scMinimal$metadata$donors
        for (bv in bvars) {
          batch_ndx <- which(scMinimal$metadata[[shuffle_within]]==bv)
          donor_meta[batch_ndx] <- as.factor(sample(donor_meta[batch_ndx]))
        }
      }
    } else {
      donor_meta <- as.factor(scMinimal$metadata$donors)
    }
    donor_sum_counts <- get_sums(cdata,donor_meta)

    # remove first row because it's all NaN
    donor_sum_counts <- donor_sum_counts[2:nrow(donor_sum_counts),]

    scMinimal$pseudobulk <- t(donor_sum_counts)
  }
  return(container)
}

#' Normalize the pseudobulked counts matrices
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param method character The normalization method to use on the pseudobulked
#' count data. Set to 'regular' to do standard normalization of dividing by
#' library size. Set to 'trim' to use edgeR trim-mean normalization, whereby counts
#' are divided by library size times a normalization factor. (default='trim')
#' @param scale_factor numeric The number that gets multiplied by fractional counts
#' during normalization of the pseudobulked data (default=10000)
#'
#' @return the project container with normalized matrices
#' @export
normalize_pseudobulk <- function(container,method='trim',scale_factor=10000) {
  for (ct in container$experiment_params$ctypes_use) {
    scMinimal <- container$scMinimal_ctype[[ct]]

    if (method=='regular') {
      scMinimal$pseudobulk <- normalize_counts(scMinimal$pseudobulk,scale_factor=scale_factor)
    } else if (method=='trim') {
      all_nf <- edgeR::calcNormFactors(scMinimal$pseudobulk)

      # divide by lib size and multiply by scale factor
      lib_sizes <- Matrix::colSums(scMinimal$pseudobulk)
      scMinimal$pseudobulk <- sweep(scMinimal$pseudobulk,MARGIN=2,lib_sizes*all_nf,FUN='/')

      # log transform result
      scMinimal$pseudobulk <- log1p(scMinimal$pseudobulk * scale_factor)
    }

    # transpose for downstream analyses
    scMinimal$pseudobulk <- t(scMinimal$pseudobulk)

    # convert to sparse matrix for downstream operations
    scMinimal$pseudobulk <- Matrix(scMinimal$pseudobulk, sparse = TRUE)
  }
  return(container)
}


#' Helper function to normalize and log-transform count data
#'
#' @param count_data matrix or sparse matrix Gene by cell matrix of counts
#' @param scale_factor numeric The number that gets multiplied by fractional counts
#' during normalization of the pseudobulked data (default=10000)
#'
#' @return the normalized, log-transformed data
#' @export
normalize_counts <- function(count_data, scale_factor=10000) {
  # divide by lib size and multiply by scale factor
  lib_sizes <- Matrix::colSums(count_data)
  tmp <- sweep(count_data,MARGIN=2,lib_sizes,FUN='/') * scale_factor

  # log transform result
  tmp <- log1p(tmp)

  return(tmp)
}


#' Get normalized variance for each gene, taking into account mean-variance trend
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#'
#' @return the project container with normalized variances values in scMinimal objects
#' @export
get_normalized_variance <- function(container) {
  # get normalized variance for each cell type
  for (ct in container$experiment_params$ctypes_use) {
    norm_variances <- norm_var_helper(container$scMinimal_ctype[[ct]])
    container$scMinimal_ctype[[ct]]$norm_variances <- norm_variances[[1]]
    container$scMinimal_ctype[[ct]]$var_pvals <- norm_variances[[2]]
  }
  return(container)
}


#' Calculates the normalized variance for each gene. This is based on a function
#' from Pagoda2
#'
#' @param scMinimal environment A sub-container for the project typically
#' consisting of gene expression data in its raw and processed forms as well
#' as metadata
#'
#' @return a vector of the normalized variance for each gene
#' @export
norm_var_helper <- function(scMinimal) {

  donor_sum_counts <- scMinimal$pseudobulk

  df <- colMeanVars(donor_sum_counts, rowSel = NULL)
  df$m <- log(df$m); df$v <- log(df$v);
  rownames(df) <- colnames(donor_sum_counts);

  # min.gene.cells <- round(nrow(scMinimal$pseudobulk)*.02)
  min.gene.cells <- 0
  vi <- which(is.finite(df$v) & df$nobs>=min.gene.cells);
  gam.k <- 5
  m <- mgcv::gam(stats::as.formula(paste0('v ~ s(m, k = ',gam.k,')')), data = df[vi,])

  df$res <- -Inf
  df$res[vi] <- stats::resid(m,type='response')
  n.obs <- df$nobs
  suppressWarnings(df$lp <- as.numeric(stats::pf(exp(df$res),n.obs,n.obs,lower.tail=FALSE,log.p=FALSE)))
  var_pvals <- log(p.adjust(df$lp,method='fdr'))
  df$lp <- log(df$lp)
  n.cells <- nrow(donor_sum_counts)
  scaled_var <- as.numeric(stats::qchisq(df$lp, n.cells-1, lower.tail = FALSE,log.p=TRUE)/n.cells) #gives same answer as inputting non-log pval and log.p=F
  names(scaled_var) <- colnames(donor_sum_counts)
  names(var_pvals) <- colnames(donor_sum_counts)


  # make sure no scaled_var values == 0 as I use it to scale the variance later
  scaled_var[is.nan(scaled_var)] <- 0  # first make any nan to 0
  min_non_zero <- min(scaled_var[scaled_var!=0])
  ndx_zero <- which(scaled_var==0)
  scaled_var[ndx_zero] <- min_non_zero

  # make log(pvals) for nan elements to be 0
  var_pvals[is.nan(var_pvals)] <- 0

  return(list(scaled_var,var_pvals))
}



#' Partition main gene by cell matrix into per cell type matrices with significantly
#' variable genes only
#' @import Matrix
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param method character The method used to select significantly variable
#' genes across donors within a cell type. Can be either "anova" to use basic
#' anova with cells grouped by donor or "norm_var" to get the top overdispersed
#' genes by normalized variance. Set to "norm_var_pvals" to use normalized variance
#' p-values as calculated in pagoda2.
#' @param thresh numeric A pvalue threshold to use for gene significance when
#' method is set to "anova" or "empir". For the method "norm_var" thresh is the
#' number of top overdispersed genes from each cell type to include.
#'
#' @return the project container with scMinimal environments added for each
#' cell type. The scMinimal environments contain expression count matrices
#' with columns as the cells belonging to the respective cell type and rows
#' as the genes which were identified as significantly variable across donors
#' in any cell type. Metadata is also present in these sub-containers for the
#' corresponding included cells and is in the "metadata" slot.
#' @export
get_ctype_vargenes <- function(container, method, thresh) {

  # set random seed to work with mclapply
  RNGkind("L'Ecuyer-CMRG")
  set.seed(container$experiment_params$rand_seed)

  ncores <- container$experiment_params$ncores

  if (method == "norm_var") {
    all_vargenes <- c()
    for (ct in container$experiment_params$ctypes_use) {
      norm_variances <- container$scMinimal_ctype[[ct]]$norm_variances
      norm_variances <- norm_variances[order(norm_variances,decreasing=TRUE)]

      # limit to overdispersed genes
      norm_variances <- norm_variances[norm_variances > 1]

      # limit to the top overdispersed genes
      if (thresh < length(norm_variances)) {
        norm_variances <- norm_variances[1:thresh]
      }
      all_vargenes <- c(all_vargenes,names(norm_variances))

      container$scMinimal_ctype[[ct]]$vargenes <- names(norm_variances)
    }

    container$all_vargenes <- unique(all_vargenes)

  } else if (method == "norm_var_pvals") {
    all_vargenes <- c()
    for (ct in container$experiment_params$ctypes_use) {
      var_pvals <- container$scMinimal_ctype[[ct]]$var_pvals

      if (thresh==1) {
        sig_var <- var_pvals[var_pvals <= log(thresh)]
      } else {
        sig_var <- var_pvals[var_pvals < log(thresh)]
      }

      all_vargenes <- c(all_vargenes,names(sig_var))

      container$scMinimal_ctype[[ct]]$vargenes <- names(sig_var)
    }

    container$all_vargenes <- unique(all_vargenes)
  } else if (method == 'anova') {
    var_res <- data.frame(matrix(ncol=3, nrow=0))
    for (ct in container$experiment_params$ctypes_use) {
      pvals <- vargenes_anova(container$scMinimal_ctype[[ct]], ncores)
      var_res <- rbind(var_res,cbind(names(pvals), pvals, rep(ct,length(pvals))))
    }
    colnames(var_res) <- c('genes', 'pvalues', 'ctypes')
    var_res <- as.data.frame(var_res)
    var_res$pvalues <- as.numeric(as.character(var_res$pvalues))
    var_res$genes <- as.character(var_res$genes)
    var_res$ctypes <- as.character(var_res$ctypes)
    var_res$padj <- p.adjust(var_res$pvalues, "fdr")
    vargenes_all <- var_res$genes[var_res$padj < thresh]

    container$all_vargenes <- unique(vargenes_all)

    # put cell type specific vargenes in slots in case want to access them later
    for (ct in container$experiment_params$ctypes_use) {
      ctype_res <- var_res[var_res$ctypes == ct,]
      ctype_vargenes <- ctype_res[ctype_res$padj < thresh, 'genes']
      container$scMinimal_ctype[[ct]]$vargenes <- ctype_vargenes
    }
  } else {
    stop('need to select one of the available options for vargenes_method parameter')
  }

  # reduce ctype data to only significantly variable genes
  container <- reduce_to_vargenes(container)

  return(container)
}

#' Compute significantly variable genes via anova
#' @importFrom stats aov cor lm p.adjust sd
#'
#' @param scMinimal environment A sub-container for the project typically
#' consisting of gene expression data in its raw and processed forms
#' @param ncores numeric Number of cores to use
#'
#' @return the raw pvalues for each gene
#' @export
vargenes_anova <- function(scMinimal, ncores) {

  dge_sparse <- normalize_counts(scMinimal$count_data)
  dge_sparse <- t(dge_sparse)

  # calculate anova for each gene
  pvals <- mclapply(1:ncol(dge_sparse),function(x) {
    tmp <- cbind(dge_sparse[,x],scMinimal$metadata$donors)

    anova_res <- aov(tmp[,1]~tmp[,2])
    pval <- summary(anova_res)[[1]][["Pr(>F)"]][[1]]
  }, mc.cores = ncores)

  names(pvals) <- colnames(dge_sparse)
  return(pvals)
}

#' Reduce each cell type's expression matrix to just the significantly variable genes
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#'
#' @return the project container with pseudobulked matrices reduced to only
#' the most variable genes
#' @export
reduce_to_vargenes <- function(container) {
  vargenes <- container$all_vargenes
  for (ct in container$experiment_params$ctypes_use) {
    scMinimal <- container$scMinimal_ctype[[ct]]
    scMinimal$pseudobulk <- scMinimal$pseudobulk[,vargenes]
  }
  return(container)
}

#' Apply ComBat batch correction
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param batch_var character A batch variable from metadata to remove
#'
#' @return the project container with pseudobulked matrices corrected for batch
#' @export
apply_combat <- function(container,batch_var) {
  for (ct in container$experiment_params$ctypes_use) {
    scMinimal <- container$scMinimal_ctype[[ct]]

    # need metadata at donor level
    metadata <- unique(scMinimal$metadata)
    rownames(metadata) <- metadata$donors
    metadata <- metadata[rownames(scMinimal$pseudobulk),]

    modcombat <- stats::model.matrix(~1, data=metadata)
    tmp <- sva::ComBat(dat=t(scMinimal$pseudobulk),
                       batch=metadata[,batch_var],
                       mod=modcombat, par.prior=TRUE,
                       prior.plots=FALSE)

    scMinimal$pseudobulk <- t(tmp)
  }
  return(container)
}


#' Scale variance across donors for each gene within each cell type
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param var_scale_power numeric Exponent of normalized variance that is
#' used for variance scaling. Variance for each gene
#' is initially set to unit variance across donors (for a given cell type).
#' Variance for each gene is then scaled by multiplying the unit scaled values
#' by each gene's normalized variance (where the effect of the mean-variance
#' dependence is taken into account) to the exponent specified here.
#' If NULL, uses var_scale_power from container$experiment_params.
#'
#' @return the project container with the variance altered for each gene within
#' the pseudobulked matrices for each cell type
#' @export
scale_variance <- function(container, var_scale_power) {

  for (ct in container$experiment_params$ctypes_use) {
    scMinimal <- container$scMinimal_ctype[[ct]]

    pb <- scMinimal$pseudobulk

    # center with unit variance
    pb <- scale(pb, center=TRUE)

    # if gene was all 0's it is now NaN, so need to change the values back
    pb[is.nan(pb)] <- 0

    norm_variances <- container$scMinimal_ctype[[ct]]$norm_variances
    scale_factor <- norm_variances[colnames(pb)]

    pb <- apply(pb,MARGIN=1,function(x) {
      x * (scale_factor ** var_scale_power)
    })
    scMinimal$pseudobulk <- t(pb)
  }

  return(container)
}

#' Create the tensor object by stacking each bulk cell type matrix
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#'
#' @return the project container with the tensor data in container$tensor_data
#' @export
stack_tensor <- function(container) {

  ctypes_use <- container$experiment_params$ctypes_use
  donors_in_all <- rownames(container$scMinimal_ctype[[1]]$pseudobulk)
  gene_order <- colnames(container$scMinimal_ctype[[1]]$pseudobulk)

  # make empty tensor of correct dimensions
  tnsr <- array(NA, dim = c(length(donors_in_all),
                            length(gene_order),
                            length(ctypes_use)))
  # fill the tensor
  for (i in 1:length(ctypes_use)) {
    ct <- ctypes_use[i]
    pb <- container$scMinimal_ctype[[ct]]$pseudobulk
    pb <- pb[donors_in_all,gene_order]
    tnsr[, ,i] <- as.matrix(pb)
  }

  container$tensor_data <- list(donors_in_all, gene_order, ctypes_use, tnsr)

  return(container)
}















