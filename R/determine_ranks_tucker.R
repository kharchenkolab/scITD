
utils::globalVariables(c("num_ranks", "rec_error", "num_iter", "run_type", "error_diff", "total_ranks", "myrsq", "n_factors", "factor_type", "total_n_factors","Status"))

#' Rank determination by svd on the tensor unfolded along each mode.
#'
#' @import ggplot2
#' @importFrom parallel mclapply
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param max_ranks_test numeric Vector of length 3 with maximum number of
#' ranks to test for donor, gene, and cell type modes in that order
#' @param shuffle_level character Either "cells" to shuffle cell-donor linkages or
#' "tensor" to shuffle values within the tensor. Currently "tensor" only works with
#' the svd method (default="cells")
#' @param shuffle_within character A metadata variable to shuffle cell-donor linkages
#' within (default=NULL)
#' @param num_iter numeric Number of null iterations (default=100)
#' @param batch_var character A batch variable from metadata to remove. No batch
#' correction applied if NULL. (default=NULL)
#' @param norm_method character The normalization method to use on the pseudobulked
#' count data. Set to 'regular' to do standard normalization of dividing by
#' library size. Set to 'trim' to use edgeR trim-mean normalization, whereby counts
#' are divided by library size times a normalization factor. (default='trim')
#' @param scale_factor numeric The number that gets multiplied by fractional counts
#' during normalization of the pseudobulked data (default=10000)
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
#'
#' @return the project container with rank determination plot in
#' container$plots$rank_determination_plot
#' @export
determine_ranks_tucker <- function(container, max_ranks_test,
                                   shuffle_level='cells', shuffle_within=NULL,
                                   num_iter=100, batch_var=NULL,
                                   norm_method='trim',
                                   scale_factor=10000,
                                   scale_var=TRUE,
                                   var_scale_power=.5) {

  # check if run tensor formation yet...
  if (is.null(container$tensor_data)) {
    stop("need to run form_tensor() first")
  }

  # set random seed
  RNGkind("L'Ecuyer-CMRG")
  set.seed(container$experiment_params$rand_seed)

  # extract needed inputs from experiment parameters
  ncores <- container$experiment_params$ncores

  # generate reconstruction errors under the null condition
  null_res <- lapply(1:num_iter, function(x) {

    if (shuffle_level == "cells") {

      # pseudobulk with shuffle donor-cell linkages
      container <- get_pseudobulk(container, shuffle=TRUE, shuffle_within=shuffle_within)

      # normalize data
      container <- normalize_pseudobulk(container, method=norm_method, scale_factor=scale_factor)

      # reduce to previously identified vargenes
      container <- reduce_to_vargenes(container)

      # apply batch correction if specified
      if (!is.null(batch_var)) {
        container <- apply_combat(container,batch_var=batch_var)
      }

      if (scale_var) {
        # scale gene expression
        container <- scale_variance(container,
                                    var_scale_power=var_scale_power)
      }

      # build the tensor
      container <- stack_tensor(container)
    }

    tmp_tnsr <- container$tensor_data[[4]]

    # get and store reconstruction errors
    if (shuffle_level == "cells") {
      rec_errors <- get_reconstruct_errors_svd(tmp_tnsr,max_ranks_test,shuffle_tensor=FALSE)
    } else {
      rec_errors <- get_reconstruct_errors_svd(tmp_tnsr,max_ranks_test,shuffle_tensor=TRUE)
    }
    return(rec_errors)
  })

  # recompute unshuffled tensor if shuffling was done at cell level
  if (shuffle_level == "cells") {
    # get actual reconstruction errors
    container <- get_pseudobulk(container, shuffle=FALSE)

    # normalize data
    container <- normalize_pseudobulk(container, method=norm_method, scale_factor=scale_factor)

    # reduce to previously identified vargenes
    container <- reduce_to_vargenes(container)

    # apply batch correction if specified
    if (!is.null(batch_var)) {
      container <- apply_combat(container,batch_var=batch_var)
    }

    if (scale_var) {
      # scale gene expression
      container <- scale_variance(container,
                                  var_scale_power=var_scale_power)
    }

    # build the tensor
    container <- stack_tensor(container)
  }

  tmp_tnsr <- container$tensor_data[[4]]

  # get and store reconstruction errors
  rec_errors_real <- get_reconstruct_errors_svd(tmp_tnsr,max_ranks_test,shuffle_tensor=FALSE)

  # plot the results
  all_line_plots <- list()
  all_bar_plots <- list()
  for (i in 1:length(max_ranks_test)) {
    line_plot <- plot_rec_errors_line_svd(rec_errors_real, null_res, i)
    bar_plot <- plot_rec_errors_bar_svd(rec_errors_real, null_res, i)
    all_line_plots[[i]] <- line_plot
    all_bar_plots[[i]] <- bar_plot
  }
  p <- ggpubr::ggarrange(all_line_plots[[1]], all_bar_plots[[1]], all_line_plots[[2]],
                         all_bar_plots[[2]], all_line_plots[[3]], all_bar_plots[[3]],
                         ncol = 2, nrow = 3)

  # save data and plot
  container$rank_determination_results <- list(rec_errors_real,null_res)
  container$plots$rank_determination_plot <- p

  return(container)
}


#' Calculate reconstruction errors using svd approach
#'
#' @param tnsr array A 3-dimensional array with dimensions of
#' donors, genes, and cell types in that order
#' @param max_ranks_test numeric Vector of length 3 with maximum number of
#' ranks to test for donor, gene, and cell type modes in that order
#' @param shuffle_tensor logical Set to TRUE to shuffle values within the tensor
#'
#' @return reconstruction errors
#' @export
get_reconstruct_errors_svd <- function(tnsr, max_ranks_test, shuffle_tensor) {
  tnsr <- rTensor::as.tensor(tnsr)

  if (shuffle_tensor) {
    # unfold along donor mode
    d_unfold <- rTensor::k_unfold(tnsr,1)@data

    # shuffle values in each column (preserves gene and ct structure)
    for (i in 1:ncol(d_unfold)) {
      d_unfold[,i] <- sample(d_unfold[,i])
    }

    # refold tensor
    tnsr <- rTensor::k_fold(d_unfold,m=1,modes=tnsr@modes)
  }

  mode_rank_errors <- list()

  for (m in 1:length(max_ranks_test)) {
    rnk_errors <- c()
    d_unfold <- rTensor::k_unfold(tnsr,m)@data

    svd_res <- svd(d_unfold)
    d <- diag(svd_res$d)
    for (rnk in 1:max_ranks_test[m]) {
      rec <- as.matrix(svd_res$u[,1:rnk]) %*% as.matrix(d[1:rnk,1:rnk]) %*% as.matrix(t(svd_res$v[,1:rnk]))
      fnorm_relative <- norm(rec - d_unfold, "F")**2 / norm(d_unfold, "F")**2

      rnk_errors <- c(rnk_errors,fnorm_relative)
    }
    mode_rank_errors[[m]] <- rnk_errors
  }

  return(mode_rank_errors)
}

#' Evaluate ability to extract sex-linked factor over varying the
#' variance scaling exponent parameter
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param min_ranks_test numeric Vector of length 3 with minimum number of
#' ranks to test for donor, gene, and cell type modes in that order
#' @param max_ranks_test numeric Vector of length 3 with maximum number of
#' ranks to test for donor, gene, and cell type modes in that order
#' @param min_power_test numeric The minimum exponent parameter value to test
#' @param max_power_test numeric The maximum parameter value to test
#' @param tucker_type character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
#' with sparsity constraints (default='regular')
#' @param rotation_type character Set to 'ica' to perform ICA rotation on resulting donor factor
#' matrix and loadings. Otherwise set to 'varimax' to perform varimax rotation. (default='ica')
#' @param norm_method character The normalization method to use on the pseudobulked
#' count data. Set to 'regular' to do standard normalization of dividing by
#' library size. Set to 'trim' to use edgeR trim-mean normalization, whereby counts
#' are divided by library size times a normalization factor. (default='trim')
#' @param scale_factor numeric The number that gets multiplied by fractional counts
#' during normalization of the pseudobulked data (default=10000)
#' @param batch_var character A batch variable from metadata to remove. No batch
#' correction applied if NULL. (default=NULL)
#'
#' @return the experiment container with plot in container$plots$var_scale_plot
#' @export
optimize_var_scale_power <- function(container, min_ranks_test, max_ranks_test,
                                     min_power_test, max_power_test, tucker_type='regular',
                                     rotation_type='ica', norm_method='trim',
                                     scale_factor=10000, batch_var=NULL) {

  # check if run tensor formation yet...
  if (is.null(container$tensor_data)) {
    stop("need to run form_tensor() first")
  }

  ncores <- container$experiment_params$ncores

  if (ncores > 10) {
    ncores <- 10
  }

  mycombos <- get_factor_combos(max_ranks_test)

  # only test those combos with max ctype factors
  mycombos <- mycombos[mycombos[,3]==max_ranks_test[3],]

  # only test combos above min ranks to test
  mycombos <- mycombos[mycombos[,1]>=min_ranks_test[1],]
  mycombos <- mycombos[mycombos[,2]>=min_ranks_test[2],]

  powers_test <- seq(min_power_test,max_power_test,by=.25)
  power_results <- list()
  # iterate through different scaling powers
  for (scale_power in powers_test) {
    print(scale_power)

    # get pseudobulk
    container <- get_pseudobulk(container, shuffle=FALSE)

    # normalize data
    container <- normalize_pseudobulk(container, method=norm_method, scale_factor=scale_factor)

    # reduce to previously identified vargenes
    container <- reduce_to_vargenes(container)

    # apply batch correction if specified
    if (!is.null(batch_var)) {
      container <- apply_combat(container,batch_var=batch_var)
    }

    # scale gene expression
    container <- scale_variance(container,
                                var_scale_power=scale_power)

    # build the tensor
    container <- stack_tensor(container)

    tensor_data <- container$tensor_data

    combo_max_f <- mclapply(1:nrow(mycombos),function(i,container,tensor_data) {
      tucker_result <- tucker_ica_helper(tensor_data, ranks=unlist(mycombos[i,1:3]),
                                         tucker_type=tucker_type,
                                         rotation_type=rotation_type)
      donor_mat <- as.data.frame(as.matrix(tucker_result[[1]]))

      meta <- container$scMinimal_full$metadata[,c('donors','sex')]
      meta <- unique(meta)
      rownames(meta) <- meta$donors
      meta$donors <- NULL

      # limit rows of meta to those of donor_mat
      meta <- meta[rownames(donor_mat),,drop=FALSE]

      # compute F statistic of association between donor scores and sex
      factor_fstats <- c()
      for (j in 1:ncol(donor_mat)) {
        tmp <- cbind(donor_mat[,j], meta)
        colnames(tmp) <- c('scores', 'sex')
        lmres <- lm(scores~sex, tmp)
        fstat <- summary(lmres)$fstatistic[[1]]
        factor_fstats <- c(factor_fstats,fstat)
      }
      best_fstat <- max(factor_fstats)
      return(best_fstat)
    }, container=container, tensor_data=tensor_data, mc.cores = ncores)

    # get names of combos
    combo_names <- sapply(1:nrow(mycombos), function(i) {
      paste(as.character(mycombos[i,1]), as.character(mycombos[i,2]),
                                  as.character(mycombos[i,3]), sep=',')
    })
    names(combo_max_f) <- combo_names

    power_results[[as.character(scale_power)]] <- combo_max_f
  }

  # plot the result
  var_scale_plot <- plot_var_scale_power(power_results)

  # save raw result and plot
  container$var_scale_results <- power_results
  container$plots$var_scale_plot <- var_scale_plot

  return(container)
}


#' Gets valid combinations of numbers of ranks to test with tucker
#'
#' @param max_ranks_test numeric Vector of length 3 with maximum number of
#' ranks to test for donor, gene, and cell type modes in that order
#'
#' @return a dataframe containing donor, gene, cell type rank combinations
#' @export
get_factor_combos <- function(max_ranks_test) {
  combos <- expand.grid(1:max_ranks_test[1], 1:max_ranks_test[2],
                        1:max_ranks_test[3])

  # remove combos that do not meet criteria of Timmerman, Kiers et al
  keeps <- c()
  for (i in 1:nrow(combos)) {
    mycombo <- combos[i,]
    t1 <- (mycombo[1]*mycombo[2]) >= mycombo[3]
    t2 <- (mycombo[2]*mycombo[3]) >= mycombo[1]
    t3 <- (mycombo[1]*mycombo[3]) >= mycombo[2]
    if (sum(t1,t2,t3)==3) {
      keeps <- c(keeps,i)
    }
  }
  combos <- combos[keeps,]
  combos$totalrank <- rowSums(combos)
  combos$fit <- NA

  return(combos)
}


#' Plot reconstruction errors as line plot for svd method
#'
#' @param real list The real reconstruction errors
#' @param shuffled list The reconstruction errors under null model
#' @param mode_to_show numeric The mode to plot the results for
#'
#' @return plot
#' @export
plot_rec_errors_line_svd <- function(real,shuffled,mode_to_show) {
  plot_res <- data.frame(matrix(ncol=2,nrow=0))
  for (i in 1:length(shuffled)) {
    shuffle_iter <- shuffled[[i]][[mode_to_show]]
    tmp <- cbind(shuffle_iter,c(1:length(shuffle_iter)),
                 rep(as.character(i),length(shuffle_iter)),rep('shuffled',length(shuffle_iter)))
    plot_res <- rbind(plot_res,tmp)
  }
  colnames(plot_res) <- c("rec_error","num_ranks","num_iter","run_type")

  # append real data
  tmp <- cbind(real[[mode_to_show]],c(1:length(shuffle_iter)),
               rep('real',length(shuffle_iter)),rep('real',length(shuffle_iter)))
  colnames(tmp) <- c("rec_error","num_ranks","num_iter","run_type")
  plot_res <- rbind(plot_res,tmp)
  plot_res$rec_error <- as.numeric(as.character(plot_res$rec_error))
  plot_res$num_ranks <- as.numeric(as.character(plot_res$num_ranks))

  if (mode_to_show==1) {
    by_interval <- 2
  } else if (mode_to_show==2) {
    by_interval <- 5
  } else {
    by_interval <- 1
  }


  p <- ggplot(plot_res,aes(x=num_ranks,y=rec_error,group=num_iter,color=run_type)) +
    geom_line() +
    geom_point() +
    ylim(0,1) +
    xlab("Number of Factors") +
    ylab("Relative Error") +
    labs(color = "Type") +
    ggtitle(paste0('Mode ', as.character(mode_to_show), ' Error')) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = seq(1, max(plot_res$num_ranks), by = by_interval))

  return(p)
}

#' Plot reconstruction errors as bar plot for svd method
#'
#' @param real list The real reconstruction errors
#' @param shuffled list The reconstruction errors under null model
#' @param mode_to_show numeric The mode to plot the results for
#'
#' @return plot
#' @export
plot_rec_errors_bar_svd <- function(real,shuffled,mode_to_show) {
  plot_res <- data.frame(matrix(ncol=2,nrow=0))
  for (i in 1:length(shuffled)) {
    shuffle_iter <- shuffled[[i]][[mode_to_show]]
    diffs <- rev(diff(rev(c(1,shuffle_iter))))
    tmp <- cbind(diffs,c(1:length(diffs)),
                 rep(as.character(i),length(diffs)),rep('shuffled',length(diffs)))
    plot_res <- rbind(plot_res,tmp)
  }

  colnames(plot_res) <- c("error_diff","num_ranks","num_iter","run_type")

  # append real data
  diff_real <- rev(diff(rev(c(1,real[[mode_to_show]]))))
  tmp <- cbind(diff_real,c(1:length(diff_real)),
               rep('real',length(diff_real)),rep('real',length(diff_real)))
  colnames(tmp) <- c("error_diff","num_ranks","num_iter","run_type")
  plot_res <- rbind(plot_res,tmp)
  plot_res$error_diff <- as.numeric(as.character(plot_res$error_diff))
  plot_res$num_ranks <- as.numeric(as.character(plot_res$num_ranks))

  # calculate summary stats
  suppressWarnings(tgc <- Rmisc::summarySE(plot_res, measurevar="error_diff",
                                           groupvars=c("num_ranks","run_type")))

  if (mode_to_show==1) {
    by_interval <- 2
  } else if (mode_to_show==2) {
    if (max(tgc$num_ranks) < 20) {
      by_interval <- 2
    } else {
      by_interval <- 5
    }
  } else {
    by_interval <- 1
  }

  p <- ggplot(tgc, aes(x=num_ranks, y=error_diff, fill=run_type)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=error_diff-sd, ymax=error_diff+sd), width=1,position=position_dodge()) +
    xlab("Number of Factors") +
    ylab("Error(n-1) - Error(n)") +
    labs(fill = "Type") +
    ggtitle(paste0('Mode ', as.character(mode_to_show), ' Error Differences')) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = seq(1, max(tgc$num_ranks), by = by_interval))


  return(p)
}


#' Plot the variance scaling effect on extraction of sex-linked factor
#'
#' @param power_results list F-Statistics for max factor-sex association
#' for Tucker decompositions over a range of variance scaling parameters
#'
#' @return plot
#' @export
plot_var_scale_power <- function(power_results) {
  power_results <- data.frame(do.call(rbind.data.frame, power_results))

  # sort columns by increasing donor factor then increasing gene factors
  ranks <- colnames(power_results)
  donor_rank <- sapply(ranks,function(x) {
    tmp <- strsplit(x,split='[.]')[[1]][[1]]
    return(as.numeric(substr(tmp,2,nchar(tmp))))
  })

  power_results <- power_results[,order(donor_rank)]

  colnames(power_results) <- sapply(colnames(power_results),function(x) {
    substr(x,2,nchar(x))
  })

  col_fun = colorRamp2(c(0, 150), c("lightgray", "red"))

  myhmap <- Heatmap(power_results, name = "F-Statistic",
                    cluster_columns = FALSE,
                    cluster_rows = FALSE, col = col_fun,
                    column_names_gp = gpar(fontsize = 10),
                    rect_gp = gpar(col = "white", lwd = 2),
                    column_title = "Rank Combination (donor.gene.celltype)",
                    row_title = "Variance Scaling Exponent",
                    row_names_side = "left", column_title_side = 'bottom')
  return(myhmap)
}


#' Get plot of r-squared values for each factor's association with a batch variable
#' to varying numbers of donor factors
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param donor_ranks_test numeric The number of donor rank values to test
#' @param gene_ranks numeric The number of gene ranks to use throughout
#' @param batch_var character The name of the batch meta variable
#' @param thresh numeric The threshold r-squared cutoff for considering a
#' factor to be a batch factor. Can be a vector of multiple values to get plots
#' at varying thresholds. (default=0.5)
#' @param tucker_type character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
#' with sparsity constraints (default='regular')
#' @param rotation_type character Set to 'ica' to perform ICA rotation on resulting donor factor
#' matrix and loadings. Otherwise set to 'varimax' to perform varimax rotation. (default='ica')
#'
#' @return plots placed in container$plots$num_batch_factors slot
#' @export
get_num_batch_ranks <- function(container, donor_ranks_test, gene_ranks, batch_var, thresh=0.5,
                                tucker_type='regular',rotation_type='ica') {
  n_ctypes <- length(container$experiment_params$ctypes_use)
  res <- data.frame(matrix(nrow=0,ncol=2))
  for (r in donor_ranks_test) {
    container <- run_tucker_ica(container, ranks=c(r,gene_ranks,n_ctypes),
                                tucker_type=tucker_type, rotation_type=rotation_type)
    var_exp <- sum(container$exp_var)
    container <- get_meta_associations(container,vars_test=c(batch_var))
    tmp <- cbind(t(container[["meta_associations"]]),rep(r,r),rep(var_exp,r))
    res <- rbind(res,tmp)
  }
  rownames(res) <- NULL
  colnames(res) <- c('myrsq','total_n_factors','var_exp')

  raw_plot <- ggplot(res,aes(x=total_n_factors,y=myrsq)) +
    geom_point() +
    xlab('Total Number of Factors') +
    ylab('Batch Association (r^2)') +
    ggtitle('Factor-Batch Associations') +
    theme(plot.title = element_text(hjust = 0.5))

  var_plot <- ggplot(res,aes(x=total_n_factors,y=var_exp)) +
    geom_line() +
    xlab('Total Number of Factors') +
    ylab('Explained Variance %') +
    ggtitle('Explained Variance') +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(0,100)

  all_plots <- list(raw_plot)
  for (i in 1:length(thresh)) {
    thr <- thresh[i]
    res2 <- data.frame(matrix(nrow=0,ncol=3))
    for (r in donor_ranks_test) {
      res_sub <- res[res$total_n_factors==r,]
      num_less <- sum(res_sub$myrsq < thr)
      num_greater <- sum(res_sub$myrsq >= thr)
      tmp <- cbind(c(num_less,num_greater),rep(r,2),c('non-batch','batch'))
      res2 <- rbind(res2,tmp)
    }
    colnames(res2) <- c('n_factors','total_n_factors','factor_type')
    res2$n_factors <- as.numeric(res2$n_factors)
    res2$total_n_factors <- as.numeric(res2$total_n_factors)

    thresh_plot <- ggplot(res2,aes(x=total_n_factors,y=n_factors,color=factor_type)) +
      geom_line() +
      xlab('Total Number of Factors') +
      ylab('Number of Factors') +
      ggtitle(paste0('Threshold r^2 = ',as.character(thr))) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = c(.15, 0.925)) +
      labs(color='Factor Type')

    all_plots[[i+1]] <- thresh_plot
  }

  all_plots[[length(all_plots)+1]] <- var_plot

  p <- ggpubr::ggarrange(plotlist=all_plots, nrow = 1)

  container$plots$num_batch_factors <- p

  return(container)
}


















