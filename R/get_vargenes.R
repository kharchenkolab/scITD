
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
  # calculate anova for each gene
  dge_sparse <- methods::as(scMinimal$data_sparse,'sparseMatrix')
  dge_sparse <- Matrix::t(dge_sparse)
  pvals <- mclapply(as.data.frame(dge_sparse),function(x) {
    tmp <- as.data.frame(cbind(x,scMinimal$metadata$donors))
    colnames(tmp) <- c('expres','donors')
    tmp$expres <- as.numeric(as.character(tmp$expres))
    anova_res <- aov(expres~donors,tmp)
    pval <- summary(anova_res)[[1]][["Pr(>F)"]][[1]]
  }, mc.cores = ncores)

  names(pvals) <- colnames(dge_sparse)
  return(pvals)
}


#' Compute significant variable genes via empirical shuffling
#' @useDynLib scITD
#' @importFrom Rcpp sourceCpp
#'
#' @param scMinimal environment A sub-container for the project typically
#' consisting of gene expression data in its raw and processed forms
#' @param num_iter numeric Number of null statistics to generate
#' @param ncores numeric Number of cores to use
#'
#' @return the raw pvalues for each gene
#' @export
vargenes_shuffle <- function(scMinimal,num_iter,ncores) {

  dge_sparse <- methods::as(scMinimal$data_sparse,'sparseMatrix')
  dge_sparse <- Matrix::t(dge_sparse)

  all_null_dists <- parallel::mclapply(1:num_iter,function(x) {
    donor_meta <- as.factor(sample(scMinimal$metadata$donors))
    donor_vars <- get_vars(dge_sparse,donor_meta,table(donor_meta))
  }, mc.cores = ncores)

  # need to convert from list to df
  all_null_dists <- do.call(rbind, all_null_dists)

  # get actual variances for each gene
  donor_meta <- as.factor(scMinimal$metadata$donors)
  donor_vars_act <- get_vars(dge_sparse,donor_meta,table(donor_meta))

  # calculate pvals
  act_vals <- log10(donor_vars_act)
  all_null_dists <- log10(all_null_dists)

  tmp <- rbind(act_vals,all_null_dists)
  pvals <- apply(tmp,MARGIN = 2,function(x) {
    y <- x[1]
    x <- x[2:length(x)]
    mask <- y <= x #number of null dist vals greater than the act val
    pval <- sum(mask)/length(x)
    return(pval)
  })

  names(pvals) <- colnames(dge_sparse)
  return(pvals)
}



























