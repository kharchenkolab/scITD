
utils::globalVariables(c("dscore", "donor_proportion", "ctypes", "AUC", "Specificity", "Precision"))

#' Compute associations between donor factor scores and donor proportions of cell subtypes
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param max_res numeric The maximum clustering resolution to use
#' @param stat_type character Either "fstat" to get F-Statistics, "adj_rsq" to get adjusted
#' R-squared values, or "adj_pval" to get adjusted pvalues.
#' @param integration_var character The meta data variable to use for creating
#' the joint embedding with Conos if not already provided in container$embedding (default=NULL)
#'
#' @return the project container with a plot of association results in
#' container$plots$subtype_prop_factor_associations and proportion plots in
#' container$plots$prop_plots_all
#' @export
get_subtype_prop_associations <- function(container,max_res,stat_type,integration_var=NULL) {
  if (!(stat_type %in% c("fstat","adj_rsq","adj_pval"))) {
    stop("stat_type parameter is not one of the three options")
  }
  
  if (is.null(integration_var)) {
    if (is.null(container$embedding)) {
      stop("need to set integration_var parameter to get an embedding")
    }
  } else {
    container <- reduce_dimensions(container,integration_var)
  }
  
  donor_scores <- container$tucker_results[[1]]

  # create dataframe to store association results
  res <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(res) <- c(stat_type,'resolution','factor','ctype')
  
  # make list to store subclustering results
  subc_all <- list()
  
  # store proportion plots
  prop_plots_all <- list()
  
  # loop through cell types
  for (ct in container$experiment_params$ctypes_use) {
    scMinimal <- container[["scMinimal_ctype"]][[ct]]

    # loop through increasing clustering resolutions
    cluster_res <- seq(.5,max_res,by=.1)
    for (r in cluster_res) {
      # run clustering
      subclusts <- get_subclusters(container,ct,r,min_cells_group=50,small_clust_action='merge')
      subclusts <- subclusts + 1 # moves subcluster index from 0 to 1
      subc_all[[ct]][[paste0('res:',as.character(r))]] <- subclusts
      
      num_subclusts <- length(unique(subclusts))
      if (num_subclusts > 1) {
        sub_meta_tmp <- scMinimal$metadata[names(subclusts),]
        
        # get donor proportions of subclusters
        donor_props <- compute_donor_props(subclusts,sub_meta_tmp)
        
        # transform from proportions to balances
        donor_balances <- coda.base::coordinates(donor_props)
        rownames(donor_balances) <- rownames(donor_props)
        
        # compute regression statistics
        reg_stats <- compute_associations(donor_balances,donor_scores,stat_type)
        
        # generate plot of donor proportions and scores
        colnames(donor_props) <- sapply(1:ncol(donor_props),function(x){paste0(ct,'_',x)})
        prop_plot <- plot_donor_props(donor_props,donor_scores,reg_stats,ctype_mapping=NULL,stat_type)
        prop_plots_all[[ct]][[paste0('res:',as.character(r))]] <- prop_plot
        
      } else {
        if (stat_type=='fstat' || stat_type=='adj_rsq') {
          reg_stats <- rep(0,ncol(container$tucker_results[[1]]))
        } else if (stat_type=='adj_pval') {
          reg_stats <- rep(1,ncol(container$tucker_results[[1]]))
        }
      }
      
      # store association results
      for (i in 1:length(reg_stats)) {
        new_row <- as.data.frame(list(reg_stats[i], r, paste0("Factor ", as.character(i)), ct),stringsAsFactors = F)
        colnames(new_row) <- colnames(res)
        res <- rbind(res,new_row)
      }
    }
  }

  # generate plot of associations
  reg_stat_plots <- plot_subclust_associations(res)

  # save results
  container$plots$subtype_prop_factor_associations <- reg_stat_plots
  container$plots$prop_plots_all <- prop_plots_all
  container$subclusters <- subc_all
  
  return(container)
}

#' Do leiden subclustering to get cell subtypes
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ctype character The cell type to do subclustering for
#' @param resolution numeric The leiden resolution to use
#' @param min_cells_group numeric The minimum allowable cluster size (default=50)
#' @param small_clust_action character Either 'remove' to remove subclusters or
#' 'merge' to merge clusters below min_cells_group threshold to the nearest cluster
#' above the size threshold (default='merge')
#'
#' @return a vector of cell subclusters
#' @export
get_subclusters <- function(container,ctype,resolution,min_cells_group=50,small_clust_action='merge') {
  con <- container$embedding
  
  # using leiden community detection
  clusts <- conos::findSubcommunities(con,method=conos::leiden.community, resolution=resolution, target.clusters=ctype)
  
  # limit clusts to just cells of the cell type
  ctype_bcodes <- rownames(container$scMinimal_ctype[[ctype]]$metadata)
  clusts <- clusts[names(clusts) %in% ctype_bcodes]
  
  if (small_clust_action=='remove') {
    # remove subclusters with less than n cells
    clust_sizes <- table(clusts)
    clusts_keep <- names(clust_sizes)[clust_sizes > min_cells_group]
    large_clusts <- clusts[clusts %in% clusts_keep]
  } else if (small_clust_action=='merge') {
    large_clusts <- merge_small_clusts(con,clusts,min_cells_group)
  }

  # change cluster names to numbers
  large_clusts <- sapply(large_clusts,function(x) {
    return(as.numeric(strsplit(x,split='_')[[1]][2]))
  })
  
  return(large_clusts)
}

#' Merge small subclusters into larger ones
#'
#' @param con conos Object for the dataset with umap projection and groups as cell types
#' @param clusts character The initially assigned subclusters by leiden clustering
#' @param min_cells_group numeric The minimum allowable cluster size
#'
#' @return the subclusters with small clusters below the size threshold merged into
#' the nearest larger cluster
merge_small_clusts <- function(con,clusts,min_cells_group) {
  # get names of large cluster
  clust_sizes <- table(clusts)
  clusts_keep <- names(clust_sizes)[clust_sizes > min_cells_group]
  clusts_merge <- names(clust_sizes)[clust_sizes <= min_cells_group]
  
  coords <- con[["embedding"]]
  
  # get centroids of large clusters
  get_centroid <- function(clust_name) {
    ndx <- names(clusts)[clusts==clust_name]
    x_y <- coords[ndx,]
    if (length(ndx)>1) {
      x_med <- stats::median(x_y[,1])
      y_med <- stats::median(x_y[,2])
      return(c(x_med,y_med))
    } else {
      return(x_y)
    }
  }
  
  main_centroids <- lapply(clusts_keep,get_centroid)
  names(main_centroids) <- clusts_keep
  small_centroids <- lapply(clusts_merge,get_centroid)
  names(small_centroids) <- clusts_merge
  
  # for each small cluster, find its nearest large cluster and assigns it's subtypes accordingly
  get_nearest_large_clust <- function(clust_name) {
    cent <- small_centroids[[clust_name]]
    c_distances <- c()
    # calculate euclidean distance to each big cluster's centroid
    for (big_clust in names(main_centroids)) {
      clust_dist <- sqrt(sum((main_centroids[[big_clust]] - cent)**2))
      c_distances[big_clust] <- clust_dist
    }
    nearest_big_clust <- names(main_centroids)[which(c_distances == min(c_distances))]
    return(nearest_big_clust)
  }
  
  
  for (cmerge in clusts_merge) {
    merge_partner <- get_nearest_large_clust(cmerge)
    clusts[clusts==cmerge] <- merge_partner 
  }
  
  return(clusts)
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
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param integration_var character The meta data variable to use for creating
#' the joint embedding with Conos.
#' @return a dataframe with umap coordinates of each cell in the dataset
#' @export
reduce_dimensions <- function(container, integration_var) {
  ncores <- container$experiment_params$ncores
  
  # some cells have been removed because donors had too few cells per ctype
  # need to make sure the full data is limited to the cells used in analysis
  all_cells <- c()
  for (ct in container$experiment_params$ctypes_use) {
    cells_in_ctype <- rownames(container$scMinimal_ctype[[ct]]$metadata)
    all_cells <- c(all_cells,cells_in_ctype)
  }
  # ### for testing only
  # metadata_copy <- container$scMinimal_full$metadata
  # data_copy <- container$scMinimal_full$data_sparse
  # metadata_copy <- metadata_copy[all_cells,]
  # data_copy <- data_copy[,all_cells]
  # ###
  
  container$scMinimal_full$metadata <- container$scMinimal_full$metadata[all_cells,]
  container$scMinimal_full$data_sparse <- container$scMinimal_full$data_sparse[,all_cells]
  
  # create a list of subsetted data matrices (one per var value)
  panel <- list()
  meta <- as.character(container$scMinimal_full$metadata[,integration_var])
  var_vals <- unique(meta)
  for (v in var_vals) {
    cell_ndx <- which(meta == v)
    panel[[v]] <- container$scMinimal_full$data_sparse[,cell_ndx]
  }
  
  # turn the list of matrices to list of pagoda2 objects
  panel.preprocessed <- lapply(panel, pagoda2::basicP2proc, n.cores=ncores,
                               min.cells.per.gene=0, n.odgenes=2e3,
                               get.largevis=FALSE, make.geneknn=FALSE)
  
  con <- conos::Conos$new(panel.preprocessed, n.cores=ncores)
  
  # build graph
  con$buildGraph()
  
  # make umap embedding
  con$embedGraph(method="UMAP", min.dist=0.01, spread=15, n.cores=ncores, min.prob.lower=1e-3)
  
  # assign ctype names to the cells
  con$findCommunities(method=conos::leiden.community, resolution=1)
  cell_assigns <- container$scMinimal_full$metadata[,"ctypes"]
  names(cell_assigns) <- rownames(container$scMinimal_full$metadata)
  con$clusters$leiden$groups <- cell_assigns[names(con$clusters$leiden$groups)]
  
  container$embedding <- con
  
  return(container)
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
    prop_model <- stats::as.formula(paste0("dscore ~ ",
                                    paste(colnames(donor_balances),collapse=" + ")))

    # # run lm
    # lmres <- stats::lm(prop_model, data=tmp)
    # 
    # # extract regression statistic
    # if (stat_type == 'fstat') {
    #   reg_stat <- summary(lmres)$fstatistic[[1]]
    # } else if (stat_type == 'adj_rsq') {
    #   reg_stat <- summary(lmres)$adj.r.squared
    # } else if (stat_type == 'adj_pval') {
    #   x <- summary(lmres)
    #   reg_stat <- stats::pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)
    # }
    
    # run robust regression 
    lmres <- MASS::rlm(prop_model, data=tmp, maxit = 100)
    lmres <- sfsmisc::f.robftest(lmres)
    
    # extract regression statistic
    if (stat_type == 'fstat') {
      reg_stat <- lmres$statistic
    } else if (stat_type == 'adj_pval') {
      reg_stat <- lmres$p.value
    }
    
    all_reg_stats <- c(all_reg_stats,reg_stat)
  }
  return(all_reg_stats)
}


#' Gets cell subtype plots including an embedding, a factor association plot, and
#' a heatmap of differentially expressed genes between subtypes
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ctype character The cell type for which subtypes are to be investigated
#' @param res numeric The clustering resolution that was used to generate
#' the clustering
#' @param factor_use numeric The factor to plot scores for
#'
#' @return The embedding plot for the cell type
#' @export
get_subclust_plots <- function(container,ctype,res,factor_use) {
  
  con <- container[["embedding"]]
  resolution_name <- paste0('res:',as.character(res))
  subclusts <- container$subclusters[[ctype]][[resolution_name]]
  
  # append large cell type name to subclusters
  subclusts <- sapply(subclusts,function(x){paste0(ctype,'_',x)})
  
  # save original embedding
  orig_embed <- con[["embedding"]]

  # save original cluster labels
  orig_clusts <- con$clusters$leiden$groups

  # limit cells in subclusts to those that we actually have scores for
  donor_scores <- container$tucker_results[[1]]
  donor_vec <- container$scMinimal_full$metadata[names(subclusts),'donors']
  subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]
  con$clusters$leiden$groups <- as.factor(subclusts)
  con[["embedding"]] <- orig_embed[names(subclusts),]
  
  # get IQR so can remove outliers
  qt_x <- stats::quantile(con[["embedding"]][,1], c(.25,.75)) 
  qt_y <- stats::quantile(con[["embedding"]][,2], c(.25,.75)) 
  iqr_x <- qt_x[2] - qt_x[1]
  iqr_y <- qt_y[2] - qt_y[1]
  outlier_up_lim_x <- qt_x[2] + 2 * iqr_x
  outlier_down_lim_x <- qt_x[1] - 2 * iqr_x 
  outlier_up_lim_y <- qt_y[2] + 2 * iqr_y
  outlier_down_lim_y <- qt_y[1] - 2 * iqr_y 
  
  subc_embed_plot <- con$plotGraph()
  subc_embed_plot <- subc_embed_plot + 
    ggtitle(paste0(ctype,' Subclusters')) + 
    xlab('UMAP 1') +
    ylab('UMAP 2') +
    xlim(outlier_down_lim_x,outlier_up_lim_x) +
    ylim(outlier_down_lim_y,outlier_up_lim_y) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title.y = element_text(size = rel(.8)),
          axis.title.x = element_text(size = rel(.8)))
  
  # make subtype association plot
  subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
  scMinimal <- container$scMinimal_ctype[[ctype]]
  sub_meta_tmp <- scMinimal$metadata[names(subclusts),]

  # get donor proportions of subclusters
  donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)

  subtype_associations <- get_indv_subtype_associations(container,donor_props,factor_use)
  
  # get directionality of associations
  for (i in 1:length(subtype_associations)) {
    subc_name <- names(subtype_associations)[i]
    subc_name <- strsplit(subc_name,split="_")[[1]][1]
    
    # top top and bottom percentile of donor score
    scores_eval <- donor_scores[,factor_use]
    cutoffs <- stats::quantile(scores_eval, c(.25, .75)) 
    donors_low <- names(scores_eval)[scores_eval < cutoffs[1]]
    donors_high <- names(scores_eval)[scores_eval > cutoffs[2]]
    
    donors_high_props <- donor_props[donors_high,subc_name]
    donors_low_props <- donor_props[donors_low,subc_name]
    
    donors_high_props_mean <- mean(donors_high_props)
    donors_low_props_mean <- mean(donors_low_props)
    
    subtype_associations[i] <- -log10(subtype_associations[i])
    
    if (donors_high_props_mean < donors_low_props_mean) {
      subtype_associations[i] <- subtype_associations[i] * -1
    }
    
  }
  
  # plot enrichment results - use pval cutoff line, make up red and down blue
  subtype_names <- sapply(1:length(subtype_associations),function(x){
    paste0(ctype,"_",x)})
  tmp <- data.frame(cbind(subtype_names,subtype_associations))
  subc_assoc_plot <- ggplot(tmp,aes(x=subtype_names,y=as.numeric(subtype_associations), 
                      fill = as.numeric(subtype_associations)>0)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values = c("lightblue", "firebrick")) +
    ylab('-log10(Adj. P-Value)') +
    xlab('') +
    geom_hline(yintercept=0, linetype="solid", color = "black") +
    geom_hline(yintercept=-log10(.05), linetype="dashed", color = "red") +
    geom_hline(yintercept=log10(.05), linetype="dashed", color = "red") +
    ggtitle(paste0("Subcluster-Factor ",factor_use," Associations")) +
    theme_bw() +
    theme(legend.position = "none", axis.title.y = element_text(size = rel(.8)),
          plot.title = element_text(hjust = 0.5)) 
    
  # get subtype DE results heamap
  myde <- con$getDifferentialGenes(groups=as.factor(subclusts),append.auc=TRUE,z.threshold=0,upregulated.only=TRUE)
  subc_de_hmap <- plotDEheatmap_conos(con, groups=as.factor(subclusts), de=myde, container,
                                 row.label.font.size=8)
  
  # make heatmap into a grob
  subc_hmap_grob <- grid::grid.grabExpr(draw(subc_de_hmap,annotation_legend_side = "bottom"))
  
  # store results in container
  container$plots$subc_plots[[paste0(ctype,"_",resolution_name)]] <- list(subc_embed_plot,
                                                                 subc_assoc_plot,
                                                                 subc_hmap_grob)

  # reset the embedding and clusters
  con$clusters$leiden$groups <- orig_clusts
  con[["embedding"]] <- orig_embed
  container$embedding <- con
  
  return(container)
}


#' Compute single subtype associations with a factor
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param donor_props matrix Donor proportions of subtypes
#' @param factor_select numeric The factor to get associations for
#'
#' @return the association statistics for each factor
#' @export
get_indv_subtype_associations <- function(container, donor_props, factor_select) {
  reg_stats_all <- list()
  for (j in 1:ncol(donor_props)) {
    # choose a column (subtype)
    subtype <- donor_props[,j,drop=FALSE]
    
    # add a second column with value of 1 - first column
    subtype <- cbind(subtype,1-subtype)
    
    # get balances
    donor_balances <- coda.base::coordinates(subtype)
    rownames(donor_balances) <- rownames(subtype)
    
    # compute regression statistics
    reg_stats <- compute_associations(donor_balances,container$tucker_results[[1]],"adj_pval")
    reg_stats_all[[paste0("K",j,"_")]] <- reg_stats
  }
  
  reg_stats_all <- unlist(reg_stats_all)
  reg_stats_all <- stats::p.adjust(reg_stats_all, method = 'fdr')
  
  parsed_name <- sapply(names(reg_stats_all),function(x){
    return(as.numeric(strsplit(x,split="_")[[1]][2]))
  })
  reg_stats_all <- reg_stats_all[parsed_name==factor_select]
  
  return(reg_stats_all)
}

#' Generate subcluster plots for multiple subclusterings at different resolutsions
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ctypes character The major cell types for which to get subcluster plots for
#' @param res numeric The subcluster resolution corresponding to the ctypes vector
#' @param factors numeric A vector of factors to get compute subtype associations
#' with. Should be same length as ctypes vector.
#'
#' @return the subtype plots for cell types at the specified resolutions
#' @export
get_all_subclust_plots <- function(container,ctypes,res,factors) {
  for (i in 1:length(ctypes)) {
    ct <- ctypes[i]
    r <- res[i]
    f <- factors[i]
    container <- get_subclust_plots(container=container,ctype=ct,res=r,factor_use=f)
  }
  return(container)
}

#' Render a figure of all subcluster plots
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#'
#' @export
render_subtype_plots <- function(container) {
  
  num_subfig_cols <- length(container$plots$subc_plots)

  select_grobs <- function(lay) {
    id <- unique(c(t(lay)))
    id[!is.na(id)]
  }

  hlay <- c()
  gs <- list()
  for (j in 1:num_subfig_cols) {
    # create single list of plots in right order
    # start by aligning first two plots for the cell type
    subc_embed_plot <- container$plots$subc_plots[[j]][[1]]
    subc_assoc_plot <- container$plots$subc_plots[[j]][[2]]
    
    # align these two plots
    ggplots_combined <- ggpubr::ggarrange(subc_embed_plot,subc_assoc_plot,
                                          ncol = 1,align='v')
    
    # set up layout for figure positioning
    start_ndx <- ((j-1)*2) + 1
    col_lay <- rbind(c(rep(start_ndx,26),NA),c(NA,NA,rep(start_ndx+1,24),NA))
    hlay <- cbind(hlay,col_lay)
    
    # store plots for layout
    gs[[start_ndx]] <- ggplots_combined
    gs[[start_ndx+1]] <- container$plots$subc_plots[[j]][[3]]
    
  }

  gridExtra::grid.arrange(grobs=gs, layout_matrix=hlay)

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
      round_digits <- 3
    } else if (stat_type=='adj_rsq') {
      plot_stat_name <- 'Adjusted R-Squared'
      round_digits <- 3
    } else if (stat_type == 'adj_pval') {
      plot_stat_name <- 'Adjusted P-Value'
      round_digits <- 5
    }

    p <- ggplot(tmp2, aes(x=dscore,y=donor_proportion,color=ctypes)) +
      geom_line() +
      ggtitle(paste0("Factor ",as.character(f))) +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(color = "Cell Type") +
      xlab("") +
      ylab("") +
      annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
               label=paste0(plot_stat_name,': ',round(significance[f],digits=round_digits)))

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

#' Get plot for associations between subcluster proportions for each major cell
#' type and each factor
#'
#' @param res data.frame Regression statistics for each subcluster analysis
#'
#' @return plots of regression statistics for each subtypes at varying clustering
#' resolutions and for each factor
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

    p <- ggplot(res_factor,aes_string(x='resolution',y=stat_type,color='ctype')) +
      geom_line() +
      xlab("") +
      ylab("") +
      labs(color = "Cell Type") +
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
                  bottom = ggpubr::text_grob("Leiden Resolution",
                                     size = 15, hjust = .7),
                  left = ggpubr::text_grob(y_axis_name, rot = 90, size = 15, hjust = .375))
  return(f_plots)
}






#' Plot a heatmap of differential genes
#' @importFrom dplyr %>%
#'
#' @param con conos (or p2) object
#' @param groups groups in which the DE genes were determined (so that the cells can be ordered correctly)
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param de differential expression result (list of data frames)
#' @param min.auc optional minimum AUC threshold
#' @param min.specificity optional minimum specificity threshold
#' @param min.precision optional minimum precision threshold
#' @param n.genes.per.cluster number of genes to show for each cluster
#' @param additional.genes optional additional genes to include (the genes will be assigned to the closest cluster)
#' @param exclude.genes an optional list of genes to exclude from the heatmap
#' @param labeled.gene.subset a subset of gene names to show (instead of all genes). Can be a vector of gene names, or a number of top genes (in each cluster) to show the names for.
#' @param expression.quantile expression quantile to show (0.98 by default)
#' @param pal palette to use for the main heatmap
#' @param ordering order by which the top DE genes (to be shown) are determined (default "-AUC")
#' @param column.metadata additional column metadata, passed either as a data.frame with rows named as cells, or as a list of named cell factors.
#' @param show.gene.clusters whether to show gene cluster color codes
#' @param remove.duplicates remove duplicated genes (leaving them in just one of the clusters)
#' @param column.metadata.colors a list of color specifications for additional column metadata, specified according to the HeatmapMetadata format. Use "clusters" slot to specify cluster colors.
#' @param show.cluster.legend whether to show the cluster legend
#' @param show_heatmap_legend whether to show the expression heatmap legend
#' @param border show borders around the heatmap and annotations
#' @param return.details if TRUE will return a list containing the heatmap (ha), but also raw matrix (x), expression list (expl) and other info to produce the heatmap on your own.
#' @param row.label.font.size font size for the row labels
#' @param order.clusters whether to re-order the clusters according to the similarity of the expression patterns (of the genes being shown)
#' @param split logical If TRUE splits the heatmap by cell type (default=FALSE)
#' @param split.gap numeric The distance to put in the gaps between split parts of the heatmap if split=TRUE (default=0)
#' @param cell.order explicitly supply cell order
#' @param averaging.window optional window averaging between neighboring cells within each group (turned off by default) - useful when very large number of cells shown (requires zoo package)
#' @param ... extra parameters are passed to pheatmap
#' @return ComplexHeatmap::Heatmap object (see return.details param for other output)
#' @export
plotDEheatmap_conos <- function(con,groups,container,de=NULL,min.auc=NULL,min.specificity=NULL,min.precision=NULL,n.genes.per.cluster=10,additional.genes=NULL,exclude.genes=NULL, labeled.gene.subset=NULL, expression.quantile=0.99,pal=grDevices::colorRampPalette(c('dodgerblue1','grey95','indianred1'))(1024),ordering='-AUC',column.metadata=NULL,show.gene.clusters=TRUE, remove.duplicates=TRUE, column.metadata.colors=NULL, show.cluster.legend=TRUE, show_heatmap_legend=FALSE, border=TRUE, return.details=FALSE, row.label.font.size=10, order.clusters=FALSE, split=FALSE, split.gap=0, cell.order=NULL, averaging.window=0, ...) {
  
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || utils::packageVersion("ComplexHeatmap") < "2.4") {
    stop("ComplexHeatmap >= 2.4 package needs to be installed to use plotDEheatmap. Please run \"devtools::install_github('jokergoo/ComplexHeatmap')\".")
  }
  
  getGeneExpression <- utils::getFromNamespace("getGeneExpression", "conos")
  
  groups <- as.factor(groups)
  
  if(is.null(de)) { # run DE
    de <- con$getDifferentialGenes(groups=groups,append.auc=TRUE,z.threshold=0,upregulated.only=TRUE)
  }
  
  # drop empty results
  de <- de[unlist(lapply(de,nrow))>0]
  
  # drop results that are not in the factor levels
  de <- de[names(de) %in% levels(groups)]
  
  # order de list to match groups order
  de <- de[order(match(names(de),levels(groups)))]
  
  
  # apply filters
  if(!is.null(min.auc)) {
    if(!is.null(de[[1]]$AUC)) {
      de <- lapply(de,function(x) x %>% dplyr::filter(AUC>min.auc))
    } else {
      warning("AUC column lacking in the DE results - recalculate with append.auc=TRUE")
    }
  }
  if(!is.null(min.specificity)) {
    if(!is.null(de[[1]]$Specificity)) {
      de <- lapply(de,function(x) x %>% dplyr::filter(Specificity>min.specificity))
    } else {
      warning("Specificity column lacking in the DE results - recalculate append.specificity.metrics=TRUE")
    }
  }
  
  if(!is.null(min.precision)) {
    if(!is.null(de[[1]]$Precision)) {
      de <- lapply(de,function(x) x %>% dplyr::filter(Precision>min.precision))
    } else {
      warning("Precision column lacking in the DE results - recalculate append.specificity.metrics=TRUE")
    }
  }
  
  #de <- lapply(de,function(x) x%>%arrange(-Precision)%>%head(n.genes.per.cluster))
  if(n.genes.per.cluster==0) { # want to show only expliclty specified genes
    if(is.null(additional.genes)) stop("if n.genes.per.cluster is 0, additional.genes must be specified")
    additional.genes.only <- TRUE;
    n.genes.per.cluster <- 30; # leave some genes to establish cluster association for the additional genes
  } else {
    additional.genes.only <- FALSE;
  }
  
  de <- lapply(de,function(x) x%>%dplyr::arrange(!!rlang::parse_expr(ordering))%>%head(n.genes.per.cluster))
  de <- de[unlist(lapply(de, nrow))>0]
  
  gns <- lapply(de,function(x) as.character(x$Gene)) %>% unlist
  sn <- function(x) stats::setNames(x,x)
  expl <- lapply(de,function(d) do.call(rbind,lapply(sn(as.character(d$Gene)),function(gene) getGeneExpression(con,gene))))
  
  # place additional genes
  if(!is.null(additional.genes)) {
    genes.to.add <- setdiff(additional.genes,unlist(lapply(expl,rownames)))
    if(length(genes.to.add)>0) {
      x <- setdiff(genes.to.add,conos::getGenes(con)); if(length(x)>0) warning('the following genes are not found in the dataset: ',paste(x,collapse=' '))
      
      age <- do.call(rbind,lapply(sn(genes.to.add),function(gene) getGeneExpression(con,gene)))
      
      # for each gene, measure average correlation with genes of each cluster
      acc <- do.call(rbind,lapply(expl,function(og) rowMeans(cor(t(age),t(og)),na.rm=TRUE)))
      acc <- acc[,apply(acc,2,function(x) any(is.finite(x))),drop=FALSE]
      acc.best <- stats::na.omit(apply(acc,2,which.max))
      
      for(i in 1:length(acc.best)) {
        gn <- names(acc.best)[i];
        expl[[acc.best[i]]] <- rbind(expl[[acc.best[i]]],age[gn,,drop=FALSE])
      }
      if(additional.genes.only) { # leave only genes that were explictly specified
        expl <- lapply(expl,function(d) d[rownames(d) %in% additional.genes,,drop=FALSE])
        expl <- expl[unlist(lapply(expl,nrow))>0]
        
      }
    }
  }
  
  # omit genes that should be excluded
  if(!is.null(exclude.genes)) {
    expl <- lapply(expl,function(x) {
      x[!rownames(x) %in% exclude.genes,,drop=FALSE]
    })
  }
  
  
  exp <- do.call(rbind,expl)
  # limit to cells that were participating in the de
  exp <- stats::na.omit(exp[,colnames(exp) %in% names(stats::na.omit(groups))])
  
  if(order.clusters) {
    # group clusters based on expression similarity (of the genes shown)
    xc <- do.call(cbind,tapply(1:ncol(exp),groups[colnames(exp)],function(ii) rowMeans(exp[,ii,drop=FALSE])))
    hc <- stats::hclust(stats::as.dist(2-cor(xc)),method='ward.D2')
    groups <- factor(groups,levels=hc$labels[hc$order])
    expl <- expl[levels(groups)]
    # re-create exp (could just reorder it)
    exp <- do.call(rbind,expl)
    exp <- stats::na.omit(exp[,colnames(exp) %in% names(stats::na.omit(groups))])
  }
  
  if(averaging.window>0) {
    # check if zoo is installed
    if(requireNamespace("zoo", quietly = TRUE)) {
      exp <- do.call(cbind,tapply(1:ncol(exp),as.factor(groups[colnames(exp)]),function(ii) {
        xa <- t(zoo::rollapply(as.matrix(t(exp[,ii,drop=FALSE])),averaging.window,mean,align='left',partial=TRUE))
        colnames(xa) <- colnames(exp)[ii]
        xa
      }))
    } else {
      warning("window averaging requires zoo package to be installed. skipping.")
    }
  }
  
  # transform expression values
  x <- t(apply(as.matrix(exp), 1, function(xp) {
    if(expression.quantile<1) {
      qs <- stats::quantile(xp,c(1-expression.quantile,expression.quantile))
      if(diff(qs)==0) { # too much, set to adjacent values
        xps <- unique(xp)
        if(length(xps)<3) { qs <- range(xp) } # only two values, just take the extremes
        xpm <- stats::median(xp)
        if(sum(xp<xpm) > sum(xp>xpm)) { # more common to have values below the median
          qs[1] <- max(xp[xp<xpm])
        } else { # more common to have values above the median
          qs[2] <- min(xps[xps>xpm]) # take the next one higher
        }
      }
      xp[xp<qs[1]] <- qs[1]
      xp[xp>qs[2]] <- qs[2]
    }
    xp <- xp-min(xp);
    if(max(xp)>0) xp <- xp/max(xp);
    xp
  }))
  
  
  
  
  if(!is.null(cell.order)) {
    o <- cell.order[cell.order %in% colnames(x)]
  } else { 
    o <- order(groups[colnames(x)])
  }
  x=x[,o]
  
  annot <- data.frame(clusters=groups[colnames(x)],row.names = colnames(x))
  
  if(!is.null(column.metadata)) {
    if(is.data.frame(column.metadata)) { # data frame
      annot <- cbind(annot,column.metadata[colnames(x),])
    } else if(is.list(column.metadata)) { # a list of factors
      annot <- cbind(annot,data.frame(do.call(cbind.data.frame,lapply(column.metadata,'[',rownames(annot)))))
    } else {
      warning('column.metadata must be either a data.frame or a list of cell-named factors')
    }
  }
  annot <- annot[,rev(1:ncol(annot)),drop=FALSE]
  
  if(is.null(column.metadata.colors))  {
    column.metadata.colors <- list();
  } else {
    if(!is.list(column.metadata.colors)) stop("column.metadata.colors must be a list in a format accepted by HeatmapAnnotation col argument")
    # reorder pallete to match the ordering in groups
    if(!is.null(column.metadata.colors[['clusters']])) {
      if(!all(levels(groups) %in% names(column.metadata.colors[['clusters']]))) {
        stop("column.metadata.colors[['clusters']] must be a named vector of colors containing all levels of the specified cell groups")
      }
      column.metadata.colors[['clusters']] <- column.metadata.colors[['clusters']][levels(groups)]
    }
  }
  
  # make sure cluster colors are defined
  if(is.null(column.metadata.colors[['clusters']])) {
    uc <- unique(annot$clusters);
    column.metadata.colors$clusters <- stats::setNames(grDevices::rainbow(length(uc)),uc)
  }
  
  tt <- unlist(lapply(expl,nrow));
  rannot <- stats::setNames(rep(names(tt),tt),unlist(lapply(expl,rownames)))
  #names(rannot) <- rownames(x);
  rannot <- rannot[!duplicated(names(rannot))]
  rannot <- rannot[names(rannot) %in% rownames(x)]
  rannot <- data.frame(clusters=factor(rannot,levels=names(expl)))
  
  if(remove.duplicates) { x <- x[!duplicated(rownames(x)),] }
  
  # draw heatmap
  ha <- ComplexHeatmap::HeatmapAnnotation(df=annot,border=border,
                                          col=column.metadata.colors,
                                          show_legend=TRUE,
                                          show_annotation_name=FALSE,
                                          annotation_legend_param = list(nrow=1))

  if(show.gene.clusters) {
    ra <- ComplexHeatmap::HeatmapAnnotation(df=rannot,which='row',show_annotation_name=FALSE, show_legend=FALSE, border=border,col=column.metadata.colors)
  } else { ra <- NULL }
  
  ## turns off ComplexHeatmap warning:
  ## `use_raster` is automatically set to TRUE for a matrix with more than
  ## 2000 columns. You can control `use_raster` argument by explicitly
  ## setting TRUE/FALSE to it.
  ## Set `ht_opt$message = FALSE` to turn off this message.
  ## 
  ht_opt$message = FALSE
  
  #ComplexHeatmap::Heatmap(x, col=pal, cluster_rows=FALSE, cluster_columns=FALSE, show_column_names=FALSE, top_annotation=ha , left_annotation=ra, column_split=groups[colnames(x)], row_split=rannot[,1], row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border=TRUE,  ...);
  if(split) {
    ha <- ComplexHeatmap::Heatmap(x, name='expression', col=pal, row_labels=convert_gn(container,rownames(x)), cluster_rows=FALSE, cluster_columns=FALSE, show_row_names=is.null(labeled.gene.subset), show_column_names=FALSE, top_annotation=ha , left_annotation=ra, border=border,  show_heatmap_legend=show_heatmap_legend, row_names_gp = grid::gpar(fontsize = row.label.font.size), column_split=groups[colnames(x)], row_split=rannot[,1], row_gap = unit(split.gap, "mm"), column_gap = unit(split.gap, "mm"), ...);
  } else {
    ha <- ComplexHeatmap::Heatmap(x, name='expression', col=pal,
                                  row_labels=convert_gn(container,rownames(x)), cluster_rows=FALSE, cluster_columns=FALSE, show_row_names=is.null(labeled.gene.subset), show_column_names=FALSE, top_annotation=ha , left_annotation=ra, border=border,  show_heatmap_legend=show_heatmap_legend, row_names_gp = grid::gpar(fontsize = row.label.font.size), ...);
  }
  if(!is.null(labeled.gene.subset)) {
    if(is.numeric(labeled.gene.subset)) {
      # select top n genes to show
      labeled.gene.subset <- unique(unlist(lapply(de,function(x) x$Gene[1:min(labeled.gene.subset,nrow(x))])))
    }
    gene.subset <- which(rownames(x) %in% labeled.gene.subset)
    labels <- rownames(x)[gene.subset];
    ha <- ha + ComplexHeatmap::rowAnnotation(link = ComplexHeatmap::anno_mark(at = gene.subset, labels = labels, labels_gp = grid::gpar(fontsize = row.label.font.size)))
    
  }
  
  if(return.details) {
    return(list(ha=ha,x=x,annot=annot,rannot=rannot,expl=expl,pal=pal,labeled.gene.subset=labeled.gene.subset))
  }
  
  return(ha)
}
















