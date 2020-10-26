
#' Get normalized variance for each gene. This takes into account the
#' mean-variance relation. Based on a Pagoda2 method...
#'
#' @param scMinimal environment A sub-container for the project typically
#' consisting of gene expression data in its raw and processed forms as well
#' as metadata
#'
#' @return a vector of the normalized variance for each gene
#' @export
get_normalized_variance <- function(scMinimal) {

  dge_sparse <- methods::as(scMinimal$data_sparse,'sparseMatrix')
  dge_sparse <- Matrix::t(dge_sparse)

  donor_meta <- as.factor(scMinimal$metadata$donors)
  donor_sum_counts <- get_sums(dge_sparse,donor_meta)

  # remove first row because it's all NaN
  donor_sum_counts <- donor_sum_counts[2:nrow(donor_sum_counts),]
  donor_sum_counts <- methods::as(donor_sum_counts,'sparseMatrix')

  # trying to normalize same way as in pagoda...
  depth <- Matrix::rowSums(donor_sum_counts)
  depthScale <- 1e3
  donor_sum_counts <- donor_sum_counts/as.numeric(depth/depthScale);

  df <- colMeanVars(donor_sum_counts, rowSel = NULL)
  df$m <- log(df$m); df$v <- log(df$v);
  rownames(df) <- colnames(donor_sum_counts);

  gam.k <- 5
  min.gene.cells <- 0
  vi <- which(is.finite(df$v) & df$nobs>=min.gene.cells);
  if(length(vi)<gam.k*1.5) { gam.k=1 };# too few genes
  if(gam.k<2) {
    m <- lm(v ~ m, data = df[vi,])
  } else {
    m <- mgcv::gam(stats::as.formula(paste0('v ~ s(m, k = ',gam.k,')')), data = df[vi,])
  }

  df$res <- -Inf;  df$res <- stats::resid(m,type='response')
  n.obs <- df$nobs;
  suppressWarnings(df$lp <- as.numeric(stats::pf(exp(df$res),n.obs,n.obs,lower.tail=FALSE,log.p=TRUE)))
  n.cells <- nrow(donor_sum_counts)
  scaled_var <- as.numeric(stats::qchisq(df$lp, n.cells-1, lower.tail = FALSE,log.p=TRUE)/n.cells)
  names(scaled_var) <- colnames(donor_sum_counts)
  return(scaled_var)
}








