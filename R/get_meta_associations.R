
#' Get metadata associations with all factors
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param vars_test character The names of meta variables to get associations for
#' @param stat_use character Set to either 'rsq' to get r-squared values or 'pval'
#' to get adjusted pvalues (default='rsq)
#'
#' @return the project container with the metadata associations in container$meta_associations
#' @export
#' 
#' @examples
#' test_container <- get_meta_associations(test_container, vars_test='lanes', stat_use='pval')
get_meta_associations <- function(container, vars_test, stat_use='rsq') {
  # check that tucker has already been run
  if (is.null(container$tucker_results)) {
    stop("Need to run tucker first")
  }

  meta <- container$scMinimal_full$metadata
  # meta <- meta[,colnames(meta)!='ctypes']
  meta <- meta[,c('donors',vars_test)]
  meta <- unique(meta)
  rownames(meta) <- meta$donors
  meta$donors <- NULL

  # vars_test <- colnames(meta)

  donor_scores <- container$tucker_results[[1]]

  res_df <- data.frame(matrix(nrow=length(vars_test),ncol=ncol(donor_scores)))
  colnames(res_df) <- sapply(1:ncol(donor_scores),function(x){
    paste0('Factor',x)
  })
  rownames(res_df) <- vars_test

  for (vtest in vars_test) {
    for (i in 1:ncol(donor_scores)) {
      factor_name <- paste0('Factor',i)

      # ensure class of the meta variable is preserved (factor -> numeric conversion messes up F-stats)
      meta_type <- class(meta[[vtest]])

      # put donor scores and meta variable in same df
      tmp <- as.data.frame(cbind(donor_scores[,i],meta[rownames(donor_scores),vtest]))
      colnames(tmp) <- c("dscore",vtest)
      if (meta_type=='factor' || meta_type=='character') {
        tmp[[vtest]] <- as.factor(tmp[[vtest]])
      } else if (meta_type=='numeric') {
        tmp[[vtest]] <- as.numeric(tmp[[vtest]])
      }

      # construct the model
      mymodel <- stats::as.formula(paste0("dscore ~ ",vtest))
      lmres <- lm(mymodel,data=tmp)
      lmres <- summary(lmres)
      rsq <- lmres$adj.r.squared
      pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],
                            lmres$fstatistic[3],lower.tail=FALSE)

      # store results
      if (stat_use=='rsq') {
        res_df[vtest,factor_name] <- rsq
      } else if (stat_use=='pval') {
        res_df[vtest,factor_name] <- pval
      }
    }
  }

  res_df <- as.matrix(res_df)
  res_df[res_df<0] <- 0

  # adjust pvalues if using them instead of rsq
  if (stat_use=='pval') {
    res_df2 <- c(res_df)
    res_df2 <- p.adjust(res_df2,method='fdr')
    res_df2 <- matrix(res_df2,ncol=ncol(res_df),nrow=nrow(res_df))
    colnames(res_df2) <- colnames(res_df)
    rownames(res_df2) <- rownames(res_df)
    res_df <- res_df2
  }

  container$meta_associations <- res_df
  return(container)
}
























