
### PROBABLY WILL NEED TO FIX THIS AND CHECK CHECK() AGAIN AFTER
utils::globalVariables(c("dscore", "donor_proportion", "ctypes"))

#' Compute associations between donor factor scores and donor proportions of cell subtypes
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param max_num_k numeric The largest number of subclusters to get for any cell type
#' @param stat_type character Either "fstat" to get F-Statistics, "adj_rsq" to get adjusted
#' R-squared values, or "adj_pval" to get adjusted pvalues.
#'
#' @return the project container with a plot of results in
#' container$plots$subtype_prop_factor_associations
#' @export
get_subtype_prop_associations <- function(container,max_num_k,stat_type) {
  if (is.null(container$scMinimal_full$umap)) {
    umap_all <- reduce_dimensions(container$scMinimal_full)
  } else {
    umap_all <- container$scMinimal_full$umap
  }

  # create dataframe to store results
  res <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(res) <- c(stat_type,'k','factor','ctype')

  # loop through cell types
  for (ct in container$experiment_params$ctypes_use) {
    scMinimal <- container[["scMinimal_ctype"]][[ct]]

    # get umap coordinates for all cells of current cell type
    umap_ctype <- umap_all[colnames(scMinimal$data_sparse),]

    # loop through increasing k values
    for (k in 2:max_num_k) {
      # run kmeans
      cell_clusters <- stats::kmeans(umap_ctype,k,nstart=20,iter.max=100,
                                     algorithm = "MacQueen")$cluster
      # get donor proportions of subclusters
      donor_props <- compute_donor_props(cell_clusters,scMinimal$metadata)

      # transform from proportions to balances
      donor_balances <- coda.base::coordinates(donor_props)
      rownames(donor_balances) <- rownames(donor_props)

      # compute regression statistics
      reg_stats <- compute_associations(donor_balances,container$tucker_results[[1]],stat_type)

      # store association results
      for (i in 1:length(reg_stats)) {
        new_row <- as.data.frame(list(reg_stats[i], k, paste0("Factor ", as.character(i)), ct),stringsAsFactors = F)
        colnames(new_row) <- colnames(res)
        res <- rbind(res,new_row)
      }
    }
  }

  # generate plot
  reg_stat_plots <- plot_subclust_associations(res)

  # save results
  container$plots$subtype_prop_factor_associations <- reg_stat_plots

  return(container)
}

#' Compute associations between donor factor scores and donor proportions of major cell types
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param stat_type character Either "fstat" to get F-Statistics, "adj_rsq" to get adjusted
#' R-squared values, or "adj_pval" to get adjusted pvalues.
#'
#' @return the project container with the results plot in container$plots$ctype_prop_factor_associations
#' @export
get_ctype_prop_associations <- function(container,stat_type) {
  scMinimal <- container$scMinimal_full
  donor_scores <- container$tucker_results[[1]]
  metadata <- scMinimal$metadata
  umap_all <- scMinimal$umap

  # map cell types to numbers temporarily
  ctypes <- unique(as.character(metadata$ctypes)) # index of this is the mapping
  cell_clusters <- sapply(as.character(metadata$ctypes),function(x){
    return(which(ctypes %in% x))
  })
  names(cell_clusters) <- rownames(metadata)

  # get matrix of donor proportions of different cell types
  donor_props <- compute_donor_props(cell_clusters,metadata)

  # transform from proportions to balances
  donor_balances <- coda.base::coordinates(donor_props)
  rownames(donor_balances) <- rownames(donor_props)

  # compute associations
  sig_res <- compute_associations(donor_balances,donor_scores,stat_type)

  # plot results
  prop_figure <- plot_donor_props(donor_props,donor_scores,sig_res,ctypes,stat_type)

  # save results
  container$plots$ctype_prop_factor_associations <- prop_figure

  return(container)
}

#' Gets umap coordinates if not already provided in container$scMinimal_full$umap
#'
#' @param scMinimal environment A sub-container for the project typically
#' consisting of gene expression data in its raw and processed forms as well
#' as metadata
#'
#' @return a dataframe with umap coordinates of each cell in the dataset
#' @export
reduce_dimensions <- function(scMinimal) {
  data_use <- t(scMinimal$data_sparse)
  data_use <- scale(data_use, center = TRUE)
  data_reduced <- gmodels::fast.prcomp(data_use, scale. = FALSE, center = FALSE)
  data_reduced <- umap::umap(data_reduced$x[,1:20])$layout
  return(data_reduced)
}

#' Get donor proportions of each cell type or subtype
#'
#' @param clusts integer Cluster assignments for each cell with names as cell barcodes
#' @param metadata data.frame The $metadata field for the given scMinimal
#'
#' @return a data.frame of cluster proportions for each donor
#' @export
compute_donor_props <- function(clusts,metadata) {
  names(clusts) <- metadata[names(clusts),"donors"]
  all_donors <- unique(as.character(metadata$donors))

  # store results in df
  donor_props <- data.frame(matrix(0,ncol=length(unique(clusts)),nrow = length(all_donors)))
  colnames(donor_props) <- sapply(1:ncol(donor_props),function(x) {
    paste0('K',as.character(x))
  })
  rownames(donor_props) <- all_donors
  for (d in all_donors) {
    tmp_clusts <- clusts[names(clusts)==d]
    counts <- table(tmp_clusts)
    names(counts) <- sapply(names(counts),function(x) {
      paste0('K',as.character(x))
    })
    for (j in 1:length(counts)) {
      donor_props[d,names(counts)[j]] <- counts[j]
    }
  }
  donor_props <- donor_props + 1 #adding pseudocount to avoid infinities when make balances
  donor_props <- t(apply(donor_props, 1, function(i) i/sum(i))) # counts -> props
  return(donor_props)
}


#' Compute associations between donor proportions
#'
#' @param donor_balances matrx The balances computed from donor cell type proportions
#' @param donor_scores data.frame The donor scores matrix from tucker results
#' @param stat_type character Either "fstat" to get F-Statistics, "adj_rsq" to get adjusted
#' R-squared values, or "adj_pval" to get adjusted pvalues.
#'
#' @return a numeric vector of F-Statistics (one for each factor)
#' @export
compute_associations <- function(donor_balances, donor_scores, stat_type) {

  all_reg_stats <- c()
  # loop through factors
  for (f in 1:ncol(donor_scores)) {
    # compute association between donor proportions and donor scores
    tmp <- as.data.frame(cbind(donor_scores[,f],donor_balances[rownames(donor_scores),]))

    if (ncol(tmp)==2) {
      colnames(tmp) <- c('dscore','ilr1')
    } else {
      colnames(tmp)[1] <- "dscore"
    }

    # construct the model
    prop_model <- paste0("dscore ~ ", paste(colnames(donor_balances),collapse=" + "))

    # run lm (need to figure out how to specify multiple explanatory vars)
    lmres <- stats::lm(prop_model, data=tmp)

    # extract regression statistic
    if (stat_type == 'fstat') {
      reg_stat <- summary(lmres)$fstatistic[[1]]
    } else if (stat_type == 'adj_rsq') {
      reg_stat <- summary(lmres)$adj.r.squared
    } else if (stat_type == 'adj_pval') {
      x <- summary(lmres)
      reg_stat <- stats::pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)
    }

    all_reg_stats <- c(all_reg_stats,reg_stat)
  }
  return(all_reg_stats)
}

#' Plot donor proportions for each factor
#'
#' @param donor_props data.frame Donor proportions as output from compute_donor_props()
#' @param donor_scores data.frame Donor scores from tucker results
#' @param significance numeric F-Statistics as output from compute_associations()
#' @param ctype_mapping character The cell types corresponding with columns of donor_props
#' @param stat_type character Either "fstat" to get F-Statistics, "adj_rsq" to get adjusted
#' R-squared values, or "adj_pval" to get adjusted pvalues.
#'
#' @return plots of donor proportions for each cell type vs donor factor scores for each factor
#' @export
plot_donor_props <- function(donor_props,donor_scores,significance,ctype_mapping=NULL,stat_type) {
  if (stat_type == 'adj_pval') {
    significance <- stats::p.adjust(significance)
  }

  all_plots <- list()
  # loop through factors
  for (f in 1:ncol(donor_scores)) {
    # compute association between donor proportions and donor scores
    tmp <- cbind(donor_scores[,f],as.data.frame(donor_props[rownames(donor_scores),]))
    colnames(tmp)[1] <- "dscore"

    # need to reshape the matrix
    tmp2 <- reshape2::melt(data=tmp, id.vars = 'dscore',
                           variable.name = 'ctypes', value.name = 'donor_proportion')

    if (!is.null(ctype_mapping)) {
      tmp2$ctypes <- sapply(tmp2$ctypes,function(x){
        return(ctype_mapping[x])
      })
    }

    if (stat_type=='fstat') {
      plot_stat_name <- 'F-Statistic'
    } else if (stat_type=='adj_rsq') {
      plot_stat_name <- 'Adjusted R-Squared'
    } else if (stat_type == 'adj_pval') {
      plot_stat_name <- 'Adjusted P-Value'
    }

    p <- ggplot(tmp2, aes(x=dscore,y=donor_proportion,color=ctypes)) +
      geom_line() +
      ggtitle(paste0("Factor ",as.character(f))) +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(color = "Cell Type") +
      xlab("") +
      ylab("") +
      annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
               label=paste0(plot_stat_name,': ',round(significance[f],digits=2)))

    all_plots[[f]] <- p
  }

  p_total <- ggpubr::ggarrange(plotlist = all_plots, ncol=1,
                               common.legend = T, legend = 'right')
  p_total <- ggpubr::annotate_figure(p_total,
                                     bottom = ggpubr::text_grob("Donor Factor Score",
                                                                size = 15, hjust = .7),
                                     left = ggpubr::text_grob("Cell Type Proportion", rot = 90, size = 15, hjust = .375))
  return(p_total)
}

#' Title
#'
#' @param res data.frame Regression statistics for each subcluster analysis
#'
#' @return plots of regression statistics for each subtypes at varying k values and
#' for each factor
#' @export
plot_subclust_associations <- function(res) {

  stat_type <- colnames(res)[1]

  # if plotting pvalues, fdr adjust and transform to -log10(pval)
  if (stat_type == 'adj_pval') {
    res[,stat_type] <- stats::p.adjust(res[,stat_type], method = 'fdr')
    res[,stat_type] <- -log10(res[,stat_type])
  }

  num_factors <- length(unique(res$factor))
  ctypes <- unique(res$ctype)
  plot_list <- list()

  for (f in 1:num_factors) {
    factor_name <- paste0("Factor ",as.character(f))
    res_factor <- res[res$factor==factor_name,]

    p <- ggplot(res_factor,aes_string(x='k',y=stat_type,color='ctype')) +
      geom_line() +
      xlab("") +
      ylab("") +
      labs(color = "Broad Cell Type") +
      ggtitle(factor_name) +
      theme(plot.title = element_text(hjust = 0.5))

    # if plotting r-squared change y-limits to 0-1
    if (stat_type == 'adj_rsq') {
      p <- p + ylim(c(-.1,1))
    }

    # if plotting -log10 pvals draw significance line
    if (stat_type == 'adj_pval') {
      p <- p + geom_hline(yintercept=-log10(.05), linetype="dashed", color = "red")
    }

    plot_list[[factor_name]] <- p

  }
  f_plots <- ggpubr::ggarrange(plotlist = plot_list, ncol=1,
                               common.legend = T, legend = 'right')
  if (stat_type=='fstat') {
    y_axis_name <- 'F-Statistic'
  } else if (stat_type=='adj_rsq') {
    y_axis_name <- 'Adjusted R-Squared'
  } else if (stat_type == 'adj_pval') {
    y_axis_name <- '-log10(Adjusted P-Value)'
  }

  f_plots <- ggpubr::annotate_figure(f_plots,
                  bottom = ggpubr::text_grob("Number of Subclusters",
                                     size = 15, hjust = .7),
                  left = ggpubr::text_grob(y_axis_name, rot = 90, size = 15, hjust = .375))
  return(f_plots)
}















