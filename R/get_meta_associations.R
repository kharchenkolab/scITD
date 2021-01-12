
#' Plot meta data associations in a heatmap
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param vars_test character The names of meta variables to get associations for
#'
#' @return the project container with the heatmap in the slot
#' container$plots$meta_associations
#' @export
get_meta_associations <- function(container, vars_test) {
  # check that tucker has already been run
  if (is.null(container$tucker_results)) {
    stop("need to run tucker first")
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
      res_df[vtest,factor_name] <- rsq
    }
  }
  
  res_df <- as.matrix(res_df)
  res_df[res_df<0] <- 0
  
  container$meta_associations <- res_df
  return(container)
}
























