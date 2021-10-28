
utils::globalVariables(c("dscore", "donor_proportion", "ctypes", "AUC", "Specificity",
                         "Precision", "subtype_names","subtype_associations","dsc",
                         "prop", "cell_types", "myx", "myy"))

#' Compute associations between donor factor scores and donor proportions of cell subtypes
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param max_res numeric The maximum clustering resolution to use. Minimum is 0.5.
#' @param stat_type character Either "fstat" to get F-Statistics, "adj_rsq" to get adjusted
#' R-squared values, or "adj_pval" to get adjusted pvalues.
#' @param integration_var character The meta data variable to use for creating
#' the joint embedding with Conos if not already provided in container$embedding (default=NULL)
#' @param min_cells_group numeric The minimum allowable size for cell subpopulations
#' (default=50)
#' @param use_existing_subc logical Set to TRUE to use existing subcluster annotations
#' (default=FALSE)
#' @param alt_ct_names character Cell type names used in clustering if different from those
#' used in the main analysis. Should match the order of container$experiment_params$ctypes_use.
#'(default=NULL)
#' @param n_col numeric The number of columns to organize the plots into (default=2)
#'
#' @return the project container with a plot of association results in
#' container$plots$subtype_prop_factor_associations
#' @export
get_subtype_prop_associations <- function(container, max_res, stat_type,
                                          integration_var=NULL, min_cells_group=50,
                                          use_existing_subc=FALSE,
                                          alt_ct_names=NULL,n_col=2) {
  if (!(stat_type %in% c("fstat","adj_rsq","adj_pval"))) {
    stop("stat_type parameter is not one of the three options")
  }

  if (is.null(integration_var)) {
    if (!use_existing_subc) {
      if (is.null(container$embedding)) {
        stop("need to set integration_var parameter to get an embedding")
      }
    }
  } else {
    container <- reduce_dimensions(container,integration_var)
  }

  # make sure that groups doesn't contain cell types not present
  container$embedding$clusters$leiden$groups <- factor(container$embedding$clusters$leiden$groups,
                                                       levels=unique(container$embedding$clusters$leiden$groups))

  donor_scores <- container$tucker_results[[1]]

  # create dataframe to store association results
  res <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(res) <- c(stat_type,'resolution','factor','ctype')

  # make list to store subclustering results
  if (use_existing_subc) {
    subc_all <- container$subclusters
  } else {
    subc_all <- list()
  }
  # loop through cell types
  for (ct in container$experiment_params$ctypes_use) {
    scMinimal <- container[["scMinimal_ctype"]][[ct]]

    # loop through increasing clustering resolutions
    cluster_res <- seq(.5,max_res,by=.1)
    for (r in cluster_res) {
      if (!use_existing_subc) {
        # run clustering
        # subclusts <- get_subclusters(container,ct,r,min_cells_group=min_cells_group,
        #                              small_clust_action='merge')
        subclusts <- get_subclusters(container,ct,r,min_cells_group=min_cells_group,
                                     small_clust_action='remove')
        subclusts <- subclusts + 1 # moves subcluster index from 0 to 1
        subc_all[[ct]][[paste0('res:',as.character(r))]] <- subclusts
      } else {
        if (!is.null(alt_ct_names)) {
          ct_ndx <- which(container$experiment_params$ctypes_use==ct)
          ct_new <- alt_ct_names[ct_ndx]
          subclusts <- container$subclusters[[ct_new]][[paste0('res:',as.character(r))]]
        } else {
          subclusts <- container$subclusters[[ct]][[paste0('res:',as.character(r))]]
        }
      }

      num_subclusts <- length(unique(subclusts))

      if (num_subclusts > 1) {
        # get cells in both metadata and subclusts
        cell_intersect <- intersect(names(subclusts),rownames(scMinimal$metadata))

        sub_meta_tmp <- scMinimal$metadata[cell_intersect,]

        # get donor proportions of subclusters
        donor_props <- compute_donor_props(subclusts,sub_meta_tmp)

        # transform from proportions to balances
        donor_balances <- coda.base::coordinates(donor_props)
        rownames(donor_balances) <- rownames(donor_props)

        # compute regression statistics
        reg_stats <- compute_associations(donor_balances,donor_scores,stat_type)

        # rename donor_props columns for generating plot of donor proportions and scores
        colnames(donor_props) <- sapply(1:ncol(donor_props),function(x){paste0(ct,'_',x)})

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

  # adjust p-values if using adj_pval stat_type
  if (stat_type=='adj_pval') {
    res$adj_pval <- p.adjust(res$adj_pval,method = 'fdr')
  }

  # generate plot of associations
  reg_stat_plots <- plot_subclust_associations(res,n_col=n_col)

  # save results
  container$plots$subtype_prop_factor_associations <- reg_stat_plots
  container$subclusters <- subc_all
  container$subc_factor_association_res <- res

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
#' @param n_col numeric The number of columns to organize the plots into (default=2)
#'
#' @return the project container with the results plot in container$plots$ctype_prop_factor_associations
#' @export
get_ctype_prop_associations <- function(container,stat_type,n_col=2) {

  # need to make sure the full data is limited to the cells used in analysis
  all_cells <- c()
  for (ct in container$experiment_params$ctypes_use) {
    cells_in_ctype <- rownames(container$scMinimal_ctype[[ct]]$metadata)
    all_cells <- c(all_cells,cells_in_ctype)
  }

  container$scMinimal_full$metadata <- container$scMinimal_full$metadata[all_cells,]
  container$scMinimal_full$count_data <- container$scMinimal_full$count_data[,all_cells]

  scMinimal <- container$scMinimal_full
  donor_scores <- container$tucker_results[[1]]
  metadata <- scMinimal$metadata

  # map cell types to numbers temporarily
  all_ctypes <- unique(as.character(metadata$ctypes)) # index of this is the mapping
  cell_clusters <- sapply(as.character(metadata$ctypes),function(x){
    return(which(all_ctypes %in% x))
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
  prop_figure <- plot_donor_props(donor_props, donor_scores, sig_res, all_ctypes,
                                  stat_type, n_col=n_col)

  # save results
  container$plots$ctype_prop_factor_associations <- prop_figure

  return(container)
}

#' Compute associations between donor factor scores and donor proportions of cell subtypes
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ctype character The cell type to get results for
#' @param res numeric The clustering resolution to retrieve
#' @param n_col numeric The number of columns to organize the plots into (default=2)
#' @param alt_name character Alternate name for the cell type used in clustering (default=NULL)
#'
#' @return the project container with the results plot in container$plots$ctype_prop_factor_associations
#' @export
get_ctype_subc_prop_associations <- function(container,ctype,res,n_col=2,alt_name=NULL) {

  scMinimal <- container$scMinimal_ctype[[ctype]]
  donor_scores <- container$tucker_results[[1]]
  metadata <- scMinimal$metadata

  if (!is.null(alt_name)) {
    cell_clusters <- container[["subclusters"]][[alt_name]][[paste0('res:',as.character(res))]]
  } else {
    cell_clusters <- container[["subclusters"]][[ctype]][[paste0('res:',as.character(res))]]
  }

  # make sure same cells are in clusters as in metadata
  cells_both <- intersect(names(cell_clusters),rownames(metadata))
  cell_clusters <- cell_clusters[cells_both]
  metadata <- metadata[cells_both,]

  # get matrix of donor proportions of different cell types
  donor_props <- compute_donor_props(cell_clusters,metadata)

  # transform from proportions to balances
  donor_balances <- coda.base::coordinates(donor_props)
  rownames(donor_balances) <- rownames(donor_props)

  # compute associations
  sig_res <- compute_associations(donor_balances,donor_scores,'adj_pval')

  # plot results
  all_ctypes <- sapply(1:ncol(donor_props), function(x) {
    paste0(ctype,"_",x)
  })
  prop_figure <- plot_donor_props(donor_props, donor_scores, sig_res, all_ctypes,
                                  'adj_pval', n_col=n_col)

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
#' @param ncores numeric The number of cores to use (default=container$experiment_params$ncores)
#'
#' @return a dataframe with umap coordinates of each cell in the dataset
#' @export
reduce_dimensions <- function(container, integration_var, ncores =container$experiment_params$ncores) {

  # some cells have been removed because donors had too few cells per ctype
  # need to make sure the full data is limited to the cells used in analysis
  all_cells <- c()
  for (ct in container$experiment_params$ctypes_use) {
    cells_in_ctype <- rownames(container$scMinimal_ctype[[ct]]$metadata)
    all_cells <- c(all_cells,cells_in_ctype)
  }

  container$scMinimal_full$metadata <- container$scMinimal_full$metadata[all_cells,]
  container$scMinimal_full$count_data <- container$scMinimal_full$count_data[,all_cells]

  # create a list of subsetted data matrices (one per var value)
  panel <- list()
  meta <- as.character(container$scMinimal_full$metadata[,integration_var])
  var_vals <- unique(meta)
  for (v in var_vals) {
    cell_ndx <- which(meta == v)
    panel[[v]] <- container$scMinimal_full$count_data[,cell_ndx]
  }

  # turn the list of matrices to list of pagoda2 objects
  panel.preprocessed <- lapply(panel, pagoda2::basicP2proc, n.cores=ncores,
                               min.cells.per.gene=0, n.odgenes=2e3,
                               get.largevis=FALSE, make.geneknn=FALSE)

  con <- conos::Conos$new(panel.preprocessed, n.cores=ncores)

  # build graph
  con$buildGraph()

  # make umap embedding
  con$embedGraph(method="UMAP", min.dist=0.01, spread=15, min.prob.lower=1e-3)

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
    if (ncol(donor_balances)==1) {
      prop_model <- stats::as.formula('ilr1 ~ dscore')
    } else {
      prop_model <- stats::as.formula(paste0("dscore ~ ",
                                             paste(colnames(donor_balances),collapse=" + ")))
    }

    if (rowSums(donor_balances)[1]==1) { # tests if table has proportions
      # testing out beta regression
      breg <- betareg::betareg(prop_model, data = tmp)
      tmp <- summary(breg)
      reg_stat <- tmp$coefficients$mean['dscore','Pr(>|z|)']
    } else { # if no proportions, then table has balances instead
      # use lm
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
    }

    all_reg_stats <- c(all_reg_stats,reg_stat)
  }
  return(all_reg_stats)
}


#' Get a figure showing cell subtype proportion associations with each factor. Combines
#' this plot with subtype UMAPs and differential expression heatmaps.
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param all_ctypes character A vector of the cell types to include
#' @param all_res numeric A vector of resolutions matching the all_ctypes parameter
#'
#' @return the figure placed in the slot container$plots$subc_fig. Note that this
#' function runs better if the number of cores in the conos object in
#' container$embedding has n.cores set to a relatively small value < 10.
#' @export
get_subclust_enr_fig <- function(container,all_ctypes,all_res) {

  # make heatmap of enrichment significance pvalues
  container <- get_subclust_enr_hmap(container,all_ctypes,all_res,1:ncol(container$tucker_results[[1]]))
  enr_hmap <- container$plots$subc_enr_hmap
  enr_hmap <- grid::grid.grabExpr(draw(enr_hmap))

  # make fig panel of umaps and heatmaps
  de_hmaps <- get_subclust_de_hmaps(container,all_ctypes,all_res)

  # generate UMAPs
  container <- get_subclust_umap(container,all_ctypes=all_ctypes,all_res=all_res)
  all_umaps <- list()
  for (j in 1:length(all_ctypes)) {
    ctype <- all_ctypes[j]
    res <- all_res[j]
    ct_res <- paste0(ctype,':',as.character(res))
    all_umaps[[j]] <- container$plots$subc_umaps[[ct_res]]
  }

  r1 <- cowplot::plot_grid(plotlist=all_umaps,nrow=1,scale = 0.97)
  r2 <- cowplot::plot_grid(plotlist=de_hmaps,nrow=1)

  fig <- cowplot::plot_grid(r1,r2,enr_hmap,ncol=1,rel_heights=c(1,1.65,1))

  container$plots$subc_fig <- fig

  return(container)

}

#' Get heatmap of subtype proportion associations for each cell type and factor combo
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param all_ctypes character A vector of the cell types to include
#' @param all_res numeric A vector of resolutions matching the all_ctypes parameter
#' @param all_factors numerc A vector of the factors to compute associations for
#'
#' @return the association heatmap object in container$plots$subc_enr_hmap
get_subclust_enr_hmap <- function(container,all_ctypes,all_res,all_factors) {

  res_df <- data.frame(matrix(ncol=length(all_factors),nrow=0))
  hmap_groupings <- c()

  for (j in 1:length(all_ctypes)) {
    ctype <- all_ctypes[j]
    res <- all_res[j]
    resolution_name <- paste0('res:',as.character(res))
    subclusts <- container$subclusters[[ctype]][[resolution_name]]

    # append large cell type name to subclusters
    subclusts <- sapply(subclusts,function(x){paste0(ctype,'_',x)})

    # limit cells in subclusts to those that we actually have scores for
    donor_scores <- container$tucker_results[[1]]
    donor_vec <- container$scMinimal_full$metadata[names(subclusts),'donors']
    subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]

    # make subtype association plot
    subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
    scMinimal <- container$scMinimal_ctype[[ctype]]
    sub_meta_tmp <- scMinimal$metadata[names(subclusts),]

    # get donor proportions of subclusters
    donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)

    tmp_df <- data.frame(matrix(ncol=length(all_factors),nrow=length(unique(subclusts))))
    rownames(tmp_df) <- rownames(tmp_df) <- sapply(1:length(unique(subclusts)),function(x){
      paste0(ctype,"_",x)})

    hmap_groupings <- c(hmap_groupings, rep(ctype,length(unique(subclusts))))

    for (factor_use in all_factors) {
      subtype_associations <- get_indv_subtype_associations(container,donor_props,factor_use)

      # get directionality of associations
      for (i in 1:length(subtype_associations)) {
        subc_name <- names(subtype_associations)[i]
        subc_name <- strsplit(subc_name,split="_")[[1]][1]

        # get top and bottom percentile of donor score
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

      tmp_df[,factor_use] <- subtype_associations
    }

    # add to the all cell types results...
    res_df <- rbind(res_df,tmp_df)
  }

  hmap_groupings <- factor(hmap_groupings,levels=all_ctypes)

  # get mask of the signs
  neg_vals <- res_df < 0

  # unsign, undo log10, adjust p-values, re log10, re sign
  res_df <- abs(res_df)
  res_df <- 10**-res_df
  res_vec <- unlist(res_df)
  res_vec <- stats::p.adjust(res_vec, method = 'fdr')
  res_df_adj <- matrix(res_vec, nrow = nrow(res_df), ncol = ncol(res_df))
  colnames(res_df_adj) <- colnames(res_df)
  rownames(res_df_adj) <- rownames(res_df)
  res_df_adj <- -log10(res_df_adj)
  res_df_adj[neg_vals] <- res_df_adj[neg_vals] * -1

  # make heatmap
  res_df_adj <- t(res_df_adj)
  rownames(res_df_adj) <- sapply(all_factors,function(x) {
    paste0('Factor',x)
  })

  col_fun = colorRamp2(c(-8, log10(.05), 0, -log10(.05), 8), c("blue",  "white", "white", "white", "red"))

  res_df_adj <- as.matrix(res_df_adj)

  p <- Heatmap(res_df_adj, name='enr',
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          column_names_gp = gpar(fontsize = 8),
          row_names_gp = gpar(fontsize = 10),
          col = col_fun, column_split = hmap_groupings,
          border=TRUE, row_names_side='left',
          cluster_column_slices=FALSE, column_gap = unit(8, "mm"))
  container$subc_associations <- res_df_adj
  container$plots$subc_enr_hmap <- p
  return(container)
}


#' Get scatter plot showing significance of associations for cell subtypes
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ctype character The cell type to plot
#' @param res numeric The subcluster resolution to use
#' @param subtype numeric The number corresponding with the subtype of the major
#' cell type to plot
#' @param factor_use numeric The factor to plot
#' @param ctype_cur character The name of the major cell type used in the main analysis
#'
#' @return the plot
#' @export
get_subclust_enr_dotplot <- function(container,ctype,res,subtype,factor_use,ctype_cur=ctype) {
  resolution_name <- paste0('res:',as.character(res))
  subclusts <- container$subclusters[[ctype]][[resolution_name]]

  # append large cell type name to subclusters
  subclusts <- sapply(subclusts,function(x){paste0(ctype,'_',x)})

  # limit cells in subclusts to those that we actually have scores for
  donor_scores <- container$tucker_results[[1]]
  cell_intersect <- intersect(names(subclusts),rownames(container$scMinimal_full$metadata))
  donor_vec <- container$scMinimal_full$metadata[cell_intersect,'donors']
  subclusts <- subclusts[cell_intersect]
  subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]

  # make subtype association plot
  subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
  scMinimal <- container$scMinimal_ctype[[ctype_cur]]
  sub_meta_tmp <- scMinimal$metadata[names(subclusts),]

  # get donor proportions of subclusters
  donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)
  donor_props <- donor_props[,subtype,drop=FALSE]
  colnames(donor_props) <- 'prop'

  # append dscores for factor 4
  donor_props2 <- cbind(donor_props,donor_scores[rownames(donor_props),factor_use])
  colnames(donor_props2)[ncol(donor_props2)] <- 'dsc'

  donor_props2 <- as.data.frame(donor_props2)
  donor_props2$dsc <- as.numeric(donor_props2$dsc)
  donor_props2$prop <- as.numeric(donor_props2$prop)

  lmres <- lm(prop~dsc,data=donor_props2)
  line_range <- seq(min(donor_props2$dsc),max(donor_props2$dsc),.001)
  line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
  line_df <- cbind.data.frame(line_range,line_dat)
  line_df <- cbind.data.frame(line_df,rep('1',nrow(line_df)))
  colnames(line_df) <- c('myx','myy')

  p <- ggplot(donor_props2,aes(x=dsc,y=prop)) +
    geom_point(alpha = 0.5,pch=19,size=2) +
    geom_line(data=line_df,aes(x=myx,y=myy)) +
    xlab(paste0('Factor ',as.character(factor_use),' Donor Score')) +
    ylab(paste0('Proportion of All ',ctype)) +
    ylim(0,1) +
    ggtitle(paste0(ctype,'_',as.character(subtype),' Proportions')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))

  return(p)
}


#' Get list of cell subtype differential expression heatmaps
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param all_ctypes character A vector of the cell types to include
#' @param all_res numeric A vector of resolutions matching the all_ctypes parameter
#'
#' @return a list of the DE heatmaps as grob objects
get_subclust_de_hmaps <- function(container,all_ctypes,all_res) {
  all_plots <- list()
  con <- container$embedding

  for (j in 1:length(all_ctypes)) {
    ctype <- all_ctypes[j]
    res <- all_res[j]
    ct_res <- paste0(ctype,':',as.character(res))
    resolution_name <- paste0('res:',as.character(res))
    if (is.null(container$plots$subtype_de[[ct_res]])) {
      subclusts <- container$subclusters[[ctype]][[resolution_name]]

      # append large cell type name to subclusters
      subclusts <- sapply(subclusts,function(x){paste0(ctype,'_',x)})

      # limit cells in subclusts to those that we actually have scores for
      donor_scores <- container$tucker_results[[1]]
      donor_vec <- container$scMinimal_full$metadata[names(subclusts),'donors']
      subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]

      # save original embedding
      orig_embed <- con[["embedding"]]

      # save original cluster labels
      orig_clusts <- con$clusters$leiden$groups

      con$clusters$leiden$groups <- as.factor(subclusts)
      con[["embedding"]] <- orig_embed[names(subclusts),]

      # get subtype DE results heamap
      myde <- con$getDifferentialGenes(groups=as.factor(subclusts),append.auc=TRUE,z.threshold=0,upregulated.only=TRUE)
      subc_de_hmap <- plotDEheatmap_conos(con, groups=as.factor(subclusts), de=myde, container,
                                          row.label.font.size=8)

      # make heatmap into a grob
      subc_hmap_grob <- grid::grid.grabExpr(draw(subc_de_hmap,annotation_legend_side = "bottom"))

      # store the plot
      container$plots$subtype_de[[ct_res]] <- subc_hmap_grob
      all_plots[[j]] <- subc_hmap_grob

      # restore embedding
      con$clusters$leiden$groups <- orig_clusts
      con[["embedding"]] <- orig_embed

    } else {
      all_plots[[j]] <- container$plots$subtype_de[[ct_res]]
    }
  }

  return(all_plots)

}

#' Get a figure to display subclusterings at multiple resolutions
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param all_ctypes character A vector of the cell types to include
#' @param all_res numeric A vector of resolutions matching the all_ctypes parameter
#' @param n_col numeric The number of columns to organize the figure into (default=3)
#'
#' @return the project container with the figure in container$plots$subc_umap_fig and
#' the individual umap plots in container$plots$subc_umaps
#' @export
get_subclust_umap <- function(container,all_ctypes,all_res,n_col=3) {

  all_plts <- list()
  plots_store <- list()
  for (i in 1:length(all_ctypes)) {
    ctype <- all_ctypes[i]
    res <- all_res[i]
    con <- container[["embedding"]]
    ct_res <- paste0(ctype,':',as.character(res))
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
    outlier_up_lim_x <- qt_x[2] + 1.5 * iqr_x
    outlier_down_lim_x <- qt_x[1] - 1.5 * iqr_x
    outlier_up_lim_y <- qt_y[2] + 1.5 * iqr_y
    outlier_down_lim_y <- qt_y[1] - 1.5 * iqr_y

    # make sure not too many points will get thrown out
    n_throw_out <- sum(con[["embedding"]][,1] > outlier_up_lim_x)
    while (n_throw_out > 100) {
      xlimits <- outlier_up_lim_x - outlier_down_lim_x
      move_by <- .05 * xlimits
      outlier_up_lim_x <- outlier_up_lim_x + move_by
      n_throw_out <- sum(con[["embedding"]][,1] > outlier_up_lim_x)
    }

    n_throw_out <- sum(con[["embedding"]][,1] < outlier_down_lim_x)
    while (n_throw_out > 100) {
      xlimits <- outlier_up_lim_x - outlier_down_lim_x
      move_by <- .05 * xlimits
      outlier_down_lim_x <- outlier_down_lim_x - move_by
      n_throw_out <- sum(con[["embedding"]][,1] < outlier_down_lim_x)
    }

    n_throw_out <- sum(con[["embedding"]][,2] > outlier_up_lim_y)
    while (n_throw_out > 100) {
      ylimits <- outlier_up_lim_y - outlier_down_lim_y
      move_by <- .05 * ylimits
      outlier_up_lim_y <- outlier_up_lim_y + move_by
      n_throw_out <- sum(con[["embedding"]][,2] > outlier_up_lim_y)
    }

    n_throw_out <- sum(con[["embedding"]][,2] < outlier_down_lim_y)
    while (n_throw_out > 100) {
      ylimits <- outlier_up_lim_y - outlier_down_lim_y
      move_by <- .05 * ylimits
      outlier_down_lim_y <- outlier_down_lim_y - move_by
      n_throw_out <- sum(con[["embedding"]][,2] < outlier_down_lim_y)
    }

    subc_embed_plot <- con$plotGraph()
    subc_embed_plot <- subc_embed_plot +
      ggtitle(paste0(ctype,' res = ',as.character(res))) +
      xlab('UMAP 1') +
      ylab('UMAP 2') +
      xlim(outlier_down_lim_x,outlier_up_lim_x) +
      ylim(outlier_down_lim_y,outlier_up_lim_y) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.title.y = element_text(size = rel(.8)),
            axis.title.x = element_text(size = rel(.8)))

    all_plts[[i]] <- subc_embed_plot
    plots_store[[ct_res]] <- subc_embed_plot

    # reset to original embedding
    con$clusters$leiden$groups <- orig_clusts
    con[["embedding"]] <- orig_embed
  }
  container$plots$subc_umaps <- plots_store
  container$plots$subc_umap_fig <- cowplot::plot_grid(plotlist=all_plts,
                                                ncol=n_col, scale = 0.95)

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
get_indv_subtype_associations <- function(container, donor_props, factor_select) {
  reg_stats_all <- list()
  for (j in 1:ncol(donor_props)) {
    prop_test <- donor_props[,j,drop=FALSE]
    colnames(prop_test) <- 'ilr1'
    rownames(prop_test) <- rownames(donor_props)
    
    # compute regression statistics
    reg_stats <- compute_associations(prop_test,container$tucker_results[[1]],"adj_pval")
    names(reg_stats) <- as.character(c(1:ncol(container$tucker_results[[1]])))
    reg_stats_all[[paste0("K",j,"_")]] <- reg_stats
  }

  reg_stats_all <- unlist(reg_stats_all)

  parsed_name <- sapply(names(reg_stats_all),function(x){
    return(as.numeric(strsplit(x,split="_.")[[1]][2]))
  })
  reg_stats_all <- reg_stats_all[parsed_name==factor_select]

  return(reg_stats_all)
}


#' Plot donor proportions for each factor
#'
#' @param donor_props data.frame Donor proportions as output from compute_donor_props()
#' @param donor_scores data.frame Donor scores from tucker results
#' @param significance numeric F-Statistics as output from compute_associations()
#' @param ctype_mapping character The cell types corresponding with columns of donor_props (default=NULL)
#' @param stat_type character Either "fstat" to get F-Statistics, "adj_rsq" to get adjusted
#' R-squared values, or "adj_pval" to get adjusted pvalues (default='adj_pval')
#' @param n_col numeric The number of columns to organize the plots into (default=2)
#'
#' @return plots of donor proportions for each cell type vs donor factor scores for each factor
plot_donor_props <- function(donor_props, donor_scores, significance,
                             ctype_mapping=NULL, stat_type='adj_pval', n_col=2) {
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

    colnames(tmp2)[2] <- 'cell_types'

    if (stat_type=='fstat') {
      plot_stat_name <- 'F-Statistic'
      round_digits <- 3
    } else if (stat_type=='adj_rsq') {
      plot_stat_name <- 'adj r-sq'
      round_digits <- 3
    } else if (stat_type == 'adj_pval') {
      plot_stat_name <- 'adj p-val'
      round_digits <- 4
    }

    p <- ggplot(tmp2, aes(x=dscore,y=donor_proportion,color=cell_types)) +
      # stat_summary(fun.data=mean_cl_normal) +
      geom_smooth(method='lm', formula= y~x) +
      ggtitle(paste0("Factor ",as.character(f))) +
      labs(color = "Cell Type") +
      xlab("Donor Factor Score") +
      ylab("Cell Type Proportion") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),legend.position="bottom") +
      # annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
      #          label=paste0(plot_stat_name,': ',round(significance[f],digits=round_digits)))
      annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
               label=paste0(plot_stat_name,': ',format(significance[f], scientific = TRUE, digits=2)))

    legend <- cowplot::get_legend(
      p + theme(legend.box.margin = margin(0, 0, 30, 0))
    )

    p <- p + theme(legend.position="none")

    all_plots[[f]] <- p
  }

  fig <- cowplot::plot_grid(plotlist=all_plots, ncol=n_col)

  fig <- cowplot::plot_grid(fig, legend, ncol = 1, rel_heights = c(1, .1))

  return(fig)
}

#' Get plot for associations between subcluster proportions for each major cell
#' type and each factor
#'
#' @param res data.frame Regression statistics for each subcluster analysis
#' @param n_col numeric The number of columns to organize the plots into (default=2)
#'
#' @return plots of regression statistics for each subtypes at varying clustering
#' resolutions and for each factor
#' @export
plot_subclust_associations <- function(res,n_col=2) {

  stat_type <- colnames(res)[1]

  if (stat_type == 'adj_pval') {
    res[,stat_type] <- -log10(res[,stat_type])
  }

  if (stat_type=='fstat') {
    y_axis_name <- 'F-Statistic'
  } else if (stat_type=='adj_rsq') {
    y_axis_name <- 'adj r-sq'
  } else if (stat_type == 'adj_pval') {
    y_axis_name <- '-log10(adj p-val)'
  }

  num_factors <- length(unique(res$factor))
  ctypes <- unique(res$ctype)
  plot_list <- list()

  for (f in 1:num_factors) {
    factor_name <- paste0("Factor ",as.character(f))
    res_factor <- res[res$factor==factor_name,]

    p <- ggplot(res_factor,aes_string(x='resolution',y=stat_type,color='ctype')) +
      geom_line() +
      xlab("Leiden Resolution") +
      ylab(y_axis_name) +
      labs(color = "Cell Type") +
      ggtitle(factor_name) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position="bottom")

    # if plotting r-squared change y-limits to 0-1
    if (stat_type == 'adj_rsq') {
      p <- p + ylim(c(-.1,1))
    }

    # if plotting -log10 pvals draw significance line
    if (stat_type == 'adj_pval') {
      p <- p + geom_hline(yintercept=-log10(.01), linetype="dashed", color = "red")
    }

    legend <- cowplot::get_legend(
      p + theme(legend.box.margin = margin(0, 0, 30, 0))
    )

    p <- p + theme(legend.position="none")

    plot_list[[factor_name]] <- p

  }

  fig <- cowplot::plot_grid(plotlist=plot_list, ncol=n_col)

  fig <- cowplot::plot_grid(fig, legend, ncol = 1, rel_heights = c(1, .1))

  return(fig)
}




#' Plot a heatmap of differential genes. Code is adapted from Conos package.
#' https://github.com/kharchenkolab/conos/blob/master/R/plot.R
#' 
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
    ha <- ComplexHeatmap::Heatmap(x, name='expression', row_title=" ", row_title_gp = gpar(fontsize = 50), col=pal, row_labels=convert_gn(container,rownames(x)), cluster_rows=FALSE, cluster_columns=FALSE, show_row_names=is.null(labeled.gene.subset), show_column_names=FALSE, top_annotation=ha , left_annotation=ra, border=border,  show_heatmap_legend=show_heatmap_legend, row_names_gp = grid::gpar(fontsize = row.label.font.size), column_split=groups[colnames(x)], row_split=rannot[,1], row_gap = unit(split.gap, "mm"), column_gap = unit(split.gap, "mm"), ...);
  } else {
    ha <- ComplexHeatmap::Heatmap(x, name='expression', col=pal,
                                  row_labels=convert_gn(container,rownames(x)),
                                  row_title=" ", row_title_gp = gpar(fontsize = 50),
                                  cluster_rows=FALSE, cluster_columns=FALSE,
                                  show_row_names=is.null(labeled.gene.subset),
                                  show_column_names=FALSE, top_annotation=ha,
                                  left_annotation=ra, border=border,
                                  show_heatmap_legend=show_heatmap_legend,
                                  row_names_gp = grid::gpar(fontsize = row.label.font.size), ...);
    # ha <- ComplexHeatmap::Heatmap(x, name='expression', col=pal,
    #                               row_labels=convert_gn(container,rownames(x)),
    #                               row_title=" ", row_title_gp = gpar(fontsize = 50),
    #                               cluster_rows=FALSE, cluster_columns=FALSE,
    #                               show_row_names=is.null(labeled.gene.subset),
    #                               show_column_names=FALSE, top_annotation=ha,
    #                               left_annotation=ra, border=border,
    #                               show_heatmap_legend=show_heatmap_legend,
    #                               width = unit(15, "cm"),
    #                               height = unit(15, "cm"),
    #                               row_names_gp = grid::gpar(fontsize = row.label.font.size), ...);
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
















