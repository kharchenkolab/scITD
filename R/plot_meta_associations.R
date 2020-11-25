


#' Plot meta data associations in a heatmap
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#'
#' @return the project container with the heatmap in the slot
#' container$plots$meta_associations
#' @export
plot_meta_associations <- function(container) {
  # check that tucker has already been run
  if (is.null(container$tucker_results)) {
    stop("need to run tucker first")
  }
  
  meta <- container$scMinimal_full$metadata
  meta <- meta[,colnames(meta)!='ctypes']
  meta <- unique(meta)
  rownames(meta) <- meta$donors
  meta$donors <- NULL
  
  vars_test <- colnames(meta)
  
  donor_scores <- container$tucker_results[[1]]
  
  res_df <- data.frame(matrix(nrow=length(vars_test),ncol=ncol(donor_scores)))
  colnames(res_df) <- sapply(1:ncol(donor_scores),function(x){
    paste0('Factor',x)
  })
  rownames(res_df) <- vars_test
  
  res_df2 <- res_df
  
  for (vtest in vars_test) {
    for (i in 1:ncol(donor_scores)) {
      factor_name <- paste0('Factor',i)
      
      # put donor scores and meta variable in same df
      tmp <- cbind(donor_scores[,i],meta[rownames(donor_scores),vtest])
      colnames(tmp) <- c("dscore",vtest)
      
      # construct the model
      mymodel <- stats::as.formula(paste0("dscore ~ ",vtest))
      lmres <- lm(mymodel,data=data.frame(tmp))
      lmres <- summary(lmres)
      rsq <- lmres$adj.r.squared
      pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],
                            lmres$fstatistic[3],lower.tail=FALSE)
      
      # store results
      res_df[vtest,factor_name] <- rsq
      res_df2[vtest,factor_name] <- pval
    }
  }
  
  # adjust p values
  res_df2 <- stats::p.adjust(as.matrix(res_df2),method = 'fdr')
  res_df2 <- matrix(res_df2,ncol=ncol(donor_scores),nrow=length(vars_test))
  res_df2 <- as.data.frame(res_df2)
  colnames(res_df2) <- colnames(res_df)
  rownames(res_df2) <- rownames(res_df)
  res_df <- as.matrix(res_df)
  res_df[res_df<0] <- 0
  res_df2 <- as.matrix(res_df2)
  
  
  # create heatmap
  col_fun = colorRamp2(c(0, 1), c("white", "red"))
  myhmap <- Heatmap(res_df, name = "adj rsq",  col = col_fun,
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    row_names_side = "left",
                    border=TRUE,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      if(res_df2[i, j] < 0.001) {
                        grid::grid.text(sprintf("%.2f***", res_df[i, j]), x, y, gp = gpar(fontsize = 10))
                      } else if (res_df2[i, j] < 0.01) {
                        grid::grid.text(sprintf("%.2f**", res_df[i, j]), x, y, gp = gpar(fontsize = 10))
                      } else if (res_df2[i, j] < 0.05) {
                        grid::grid.text(sprintf("%.2f*", res_df[i, j]), x, y, gp = gpar(fontsize = 10))
                      } else {
                        grid::grid.text(sprintf("%.2f", res_df[i, j]), x, y, gp = gpar(fontsize = 10))
                      }
                    })
  container$plots$meta_associations <- myhmap
  return(container)
}
























