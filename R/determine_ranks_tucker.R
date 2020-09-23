
#' Rank determination method
#' @import ggplot2
#' @importFrom parallel mclapply
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param max_ranks_test numeric Vector of length 3 with maximum number of
#' ranks to test for donor, gene, and cell type modes in that order
#' @param method character The method to use for rank determination. For 'svd',
#' reconstruction errors are computed using svd on the unfolded tensor along each
#' mode. For 'tucker', reconstruction errors are computed using tucker decomposition
#' (default='svd')
#' @param num_iter numeric Number of null iterations (default=100)
#'
#' @return the project container with rank determination plot and results
#' added into slots
#' @export
determine_ranks_tucker <- function(container, max_ranks_test, method='svd', num_iter=100) {

  # check that var_scale_power has been set if scale_var is TRUE
  if (container$experiment_params$scale_var && is.null(container$experiment_params$var_scale_power)) {
    stop("Need to set variance scaling power parameter, var_scale_power. Use set_experiment_params()")
  }

  # set random seed
  RNGkind("L'Ecuyer-CMRG")
  set.seed(container$experiment_params$rand_seed)

  # extract needed inputs from experiment parameters
  ncores <- container$experiment_params$ncores

  # generate reconstruction errors under the null condition
  null_res <- lapply(1:num_iter, function(x) {

    # first get donor means shuffled
    container <- collapse_by_donors(container, shuffle=T)

    # form a tensor
    container <- form_tensor(container)
    tmp_tnsr <- container$tensor_data[[4]]

    # get and store reconstruction errors
    if (method == 'svd') {
      rec_errors <- get_reconstruct_errors_svd(tmp_tnsr,max_ranks_test)
    } else if (method == 'tucker') {
      rec_errors <- get_reconstruct_errors_tucker(tmp_tnsr,max_ranks_test,ncores)
    }
    return(rec_errors)
  })

  # get actual reconstruction errors
  container <- collapse_by_donors(container, shuffle=F)

  # form a tensor
  container <- form_tensor(container)
  tmp_tnsr <- container$tensor_data[[4]]

  # get and store reconstruction errors
  if (method == 'svd') {
    rec_errors_real <- get_reconstruct_errors_svd(tmp_tnsr,max_ranks_test)
  } else if (method == 'tucker') {
    rec_errors_real <- get_reconstruct_errors_tucker(tmp_tnsr,max_ranks_test,ncores)
    print(rec_errors_real)
  }

  # plot the results
  if (method == 'svd') {
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
  } else if (method == 'tucker') {
    line_plot <- plot_rec_errors_line_tucker(rec_errors_real, null_res)
    bar_plot <- plot_rec_errors_bar_tucker(rec_errors_real, null_res)
    p <- ggpubr::ggarrange(line_plot, bar_plot,
              ncol = 2, nrow = 1)
  }

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
#' @return reconstruction errors
#' @export
get_reconstruct_errors_svd <- function(tnsr, max_ranks_test) {
  # do svd to rnk ranks
  mode_rank_errors <- list()
  for (m in 1:length(max_ranks_test)) {
    rnk_errors <- c()
    d_unfold <- rTensor::k_unfold(rTensor::as.tensor(tnsr),m)@data
    svd_res <- svd(d_unfold)
    d <- diag(svd_res$d)
    for (rnk in 1:max_ranks_test[m]) {
      rec <- as.matrix(svd_res$u[,1:rnk]) %*% as.matrix(d[1:rnk,1:rnk]) %*% as.matrix(t(svd_res$v[,1:rnk]))
      fnorm_relative <- norm(rec - d_unfold, "F") / norm(d_unfold, "F")

      rnk_errors <- c(rnk_errors,fnorm_relative)
    }
    mode_rank_errors[[m]] <- rnk_errors
  }

  return(mode_rank_errors)
}

#' Calculate reconstruction errors using Tucker approach
#'
#' @param tnsr array A 3-dimensional array with dimensions of
#' donors, genes, and cell types in that order
#' @param max_ranks_test numeric Vector of length 3 with maximum number of
#' ranks to test for donor, gene, and cell type modes in that order
#' @param ncores numeric The number of cores to use
#'
#' @return the reconstruction errors
#' @export
get_reconstruct_errors_tucker <- function(tnsr,max_ranks_test,ncores) {
  mycombos <- get_factor_combos(max_ranks_test)

  fits <- mclapply(1:nrow(mycombos),function(i) {
    invisible(utils::capture.output(
      tucker_decomp <- rTensor::tucker(rTensor::as.tensor(tnsr), ranks=unlist(mycombos[i,1:3]))
    ))
    return(tucker_decomp$norm_percent)
  },mc.cores = ncores)
  mycombos$fit <- unlist(fits)

  # keep only max fit for a given totalrank
  totalrank_vals <- unique(mycombos$totalrank)
  totalrank_vals <- totalrank_vals[order(totalrank_vals,decreasing=F)]
  new_combos <- data.frame(matrix(nrow=0,ncol=5))
  for (i in totalrank_vals) {
    sets_for_totalrank <- mycombos[mycombos$totalrank==i,]
    set_keep <- which(sets_for_totalrank$fit==max(sets_for_totalrank$fit))
    new_combos <- rbind(new_combos,sets_for_totalrank[set_keep,])
  }

  return(new_combos)
}

#' Evaluate ability to extract sex-linked factor over varying the
#' variance scaling exponent parameter
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param max_ranks_test numeric Vector of length 3 with maximum number of
#' ranks to test for donor, gene, and cell type modes in that order
#'
#' @return the experiment container with added plot and results
#' @export
optimize_var_scale_power <- function(container, max_ranks_test) {
  # check that scale_var parameter is true
  if (!container$experiment_params$scale_var) {
    stop("scale_var parameter from container$experiment_params is NULL. This
         function is only needed when variance scaling is used")
  }

  # extract needed inputs from experiment parameters
  rotate_modes <- container$experiment_params$rotate_modes
  ncores <- container$experiment_params$ncores

  if (ncores > 10) {
    ncores <- 10
  }

  mycombos <- get_factor_combos(max_ranks_test)

  # only test those combos with max ctype factors
  mycombos <- mycombos[mycombos[,3]==max_ranks_test[3],]

  # get donor means
  container <- collapse_by_donors(container, shuffle=F)

  powers_test <- seq(.5,2,by=.25)
  power_results <- list()
  # iterate through different scaling powers
  for (scale_power in powers_test) {
    print(scale_power)

    # form tensor with appropriate scaling
    container <- form_tensor(container, var_scale_power=scale_power)

    tensor_data <- container$tensor_data

    combo_max_f <- mclapply(1:nrow(mycombos),function(i,container,tensor_data) {
      tucker_result <- tucker_ica_helper(tensor_data, ranks=unlist(mycombos[i,1:3]),
                                         rotate_modes = rotate_modes)
      donor_mat <- as.data.frame(as.matrix(tucker_result[[1]]))

      meta <- container$scMinimal_full$metadata[,c('donors','sex')]
      meta <- unique(meta)
      rownames(meta) <- meta$donors
      meta$donors <- NULL

      # limit rows of meta to those of donor_mat
      meta <- meta[rownames(donor_mat),,drop=F]

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


  p <- ggplot(plot_res,aes(x=num_ranks,y=rec_error,group=num_iter,color=run_type)) +
    geom_line() +
    geom_point() +
    ylim(0,1) +
    xlab("Number of Factors") +
    ylab("Relative Error") +
    labs(color = "Type") +
    scale_x_continuous(breaks=c(1:length(real[[mode_to_show]]))) +
    ggtitle(paste0('Mode ', as.character(mode_to_show), ' Error')) +
    theme(plot.title = element_text(hjust = 0.5))

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
  tgc <- Rmisc::summarySE(plot_res, measurevar="error_diff",
                   groupvars=c("num_ranks","run_type"))

  p <- ggplot(tgc, aes(x=num_ranks, y=error_diff, fill=run_type)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=error_diff-sd, ymax=error_diff+sd), width=1,position=position_dodge()) +
    xlab("Number of Factors") +
    ylab("Error(n-1) - Error(n)") +
    labs(fill = "Type") +
    scale_x_continuous(breaks=c(1:length(real[[mode_to_show]]))) +
    ggtitle(paste0('Mode ', as.character(mode_to_show), ' Error Differences')) +
    theme(plot.title = element_text(hjust = 0.5))

  return(p)
}


#' Plot reconstruction errors as line plot for tucker method
#'
#' @param real list The real reconstruction errors
#' @param shuffled list The reconstruction errors under null model
#'
#' @return plot
#' @export
plot_rec_errors_line_tucker <- function(real,shuffled) {
  plot_res <- data.frame(matrix(ncol=4,nrow=0))
  for (i in 1:length(shuffled)) {
    shuffle_iter <- shuffled[[i]][,4:5]
    tmp <- cbind(shuffle_iter,
                 rep(as.character(i),nrow(shuffle_iter)),rep('shuffled',nrow(shuffle_iter)))
    plot_res <- rbind(plot_res,tmp)
  }
  colnames(plot_res) <- c("total_ranks","rec_error","num_iter","run_type")

  # append real data
  tmp <- cbind(real[,4:5],
               rep('real',nrow(shuffle_iter)),rep('real',nrow(shuffle_iter)))
  colnames(tmp) <- c("total_ranks","rec_error","num_iter","run_type")
  plot_res <- rbind(plot_res,tmp)
  plot_res$rec_error <- as.numeric(as.character(plot_res$rec_error))
  plot_res$total_ranks <- as.numeric(as.character(plot_res$total_ranks))

  plot_res$rec_error <- (100 - plot_res$rec_error) / 100

  p <- ggplot(plot_res,aes(x=total_ranks,y=rec_error,group=num_iter,color=run_type)) +
    geom_line() +
    geom_point() +
    ylim(0,1) +
    xlab("Total Ranks") +
    ylab("Relative Error") +
    labs(color = "Type") +
    scale_x_continuous(breaks=shuffle_iter[,1]) +
    ggtitle('Error') +
    theme(plot.title = element_text(hjust = 0.5))

  return(p)
}

#' Plot reconstruction errors as bar plot for tucker method
#'
#' @param real list The real reconstruction errors
#' @param shuffled list The reconstruction errors under null model
#'
#' @return plot
#' @export
plot_rec_errors_bar_tucker <- function(real,shuffled) {
  plot_res <- data.frame(matrix(ncol=4,nrow=0))
  for (i in 1:length(shuffled)) {
    shuffle_iter <- shuffled[[i]][,4:5]
    shuffle_iter$dif <- NA
    shuffle_iter$dif[1] <- shuffle_iter$fit[1]
    for (j in 2:nrow(shuffle_iter)) {
      shuffle_iter$dif[j] <- shuffle_iter$fit[j] - shuffle_iter$fit[j-1]
    }
    tmp <- cbind(shuffle_iter[,c(1,3)],
                 rep(as.character(i),nrow(shuffle_iter)),rep('shuffled',nrow(shuffle_iter)))
    plot_res <- rbind(plot_res,tmp)
  }
  colnames(plot_res) <- c("total_ranks","error_diff","num_iter","run_type")

  # append real data
  real <- real[,4:5]
  real$dif <- NA
  real$dif[1] <- real$fit[1]
  for (i in 2:nrow(real)) {
    real$dif[i] <- real$fit[i] - real$fit[i-1]
  }
  tmp <- cbind(real[,c(1,3)],
               rep('real',nrow(real)),rep('real',nrow(real)))
  colnames(tmp) <- c("total_ranks","error_diff","num_iter","run_type")
  plot_res <- rbind(plot_res,tmp)
  plot_res$error_diff <- as.numeric(as.character(plot_res$error_diff))
  plot_res$total_ranks <- as.numeric(as.character(plot_res$total_ranks))

  plot_res$error_diff <- plot_res$error_diff / 100

  # calculate summary stats
  tgc <- Rmisc::summarySE(plot_res, measurevar="error_diff",
                   groupvars=c("total_ranks","run_type"))

  p <- ggplot(tgc, aes(x=total_ranks, y=error_diff, fill=run_type)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=error_diff-sd, ymax=error_diff+sd), width=1,position=position_dodge()) +
    xlab("Total Ranks") +
    ylab("Error(n-1) - Error(n)") +
    labs(fill = "Type") +
    scale_x_continuous(breaks=shuffle_iter[,1]) +
    ggtitle('Error Differences') +
    theme(plot.title = element_text(hjust = 0.5))

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
    substr(x,2,2)
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





















