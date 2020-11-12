
utils::globalVariables(c("dscore", "donor_proportion", "ctypes"))

#' Compute associations between donor factor scores and donor proportions of cell subtypes
#' @import conos
#' @importFrom pagoda2 basicP2proc
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
#' @return
#' @export
get_subclusters <- function(container,ctype,resolution,min_cells_group=50,small_clust_action='merge') {
  con <- container$embedding
  
  # using leiden community detection
  clusts <- findSubcommunities(con,method=leiden.community, resolution=resolution, target.clusters=ctype)
  
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
#' @return
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
      # x_mean <- mean(x_y[,1])
      # y_mean <- mean(x_y[,2])
      x_mean <- median(x_y[,1])
      y_mean <- median(x_y[,2])
      return(c(x_mean,y_mean))
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
  panel.preprocessed <- lapply(panel, basicP2proc, n.cores=ncores,
                               min.cells.per.gene=0, n.odgenes=2e3,
                               get.largevis=FALSE, make.geneknn=FALSE)
  
  con <- Conos$new(panel.preprocessed, n.cores=ncores)
  
  # build graph
  con$buildGraph()
  
  # make umap embedding
  con$embedGraph(method="UMAP", min.dist=0.01, spread=15, n.cores=ncores, min.prob.lower=1e-3)
  
  # assign ctype names to the cells
  con$findCommunities(method=leiden.community, resolution=1)
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

#' Run differential expression with cell subtypes included
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ctype character The cell type for which subtypes are to be investigated
#' @param resolution numeric The clustering resolution that was used to generate
#' the clustering
#' @param plot_subclusters logical TRUE to generate an embedding plot colored by
#' clusters and subclusters for the specified ctype (default=TRUE)
#'
#' @return the project container with the DE results in 
#' container$subcluster_de$<ctype>$<resolution> and plot results located in
#' container$plots$subclusters$<ctype>$<resolution>
#' @export
run_subcluster_de <- function(container,ctype,resolution,plot_subclusters=TRUE) {
  ncores <- container$experiment_params$ncores
  
  con <- container[["embedding"]]
  resolution_name <- paste0('res:',as.character(resolution))
  subclusts <- container$subclusters[[ctype]][[resolution_name]]
  
  orig_groups <- con$clusters$leiden$groups
  new_groups <- replace_groups_with_subclusts(con,subclusts,ctype)
  
  con$clusters$leiden$groups <- new_groups
  de.info <- con$getDifferentialGenes(n.cores=ncores, append.auc=FALSE)
  
  # convert gene names in de results dataframe
  for (de_ctype in names(de.info)) {
    de.info[[de_ctype]]$Gene <- convert_gn(container,de.info[[de_ctype]]$Gene)
  }
  
  if (plot_subclusters) {
    subc_plot <- con$plotGraph()
    container$plots$subclusters[[ctype]][[resolution_name]] <- subc_plot
  }
  
  con$clusters$leiden$groups <- orig_groups
  
  container$subcluster_de[[ctype]][[resolution_name]] <- de.info
  return(container)
}

#' Title
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ctype character The cell type for which subtypes are to be investigated
#' @param resolution numeric The clustering resolution that was used to generate
#' the clustering
#' @param use_scores logical Set to TRUE to plot donor scores of a factor from 
#' tucker results as colors on the umap (default=FALSE)
#' @param factor_use numeric The factor to plot scores of if use_scores is TRUE
#' (default=NULL)
#'
#' @return The embedding plot for the cell type
#' @export
plot_subclusts_only <- function(container,ctype,resolution,use_scores=FALSE,
                                factor_use=NULL) {
  
  con <- container[["embedding"]]
  resolution_name <- paste0('res:',as.character(resolution))
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
  
  # REMOVE OUTLIERS HERE!
  
  subc_plot <- con$plotGraph()
  
  # make subtype association plot
  subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
  scMinimal <- container$scMinimal_ctype[[ctype]]
  sub_meta_tmp <- scMinimal$metadata[names(subclusts),]

  # get donor proportions of subclusters
  donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)

  subtype_associations <- get_indv_subtype_associations(container,donor_props)
  
  
  # if (use_scores) {
  #   if (is.null(factor_use)) {
  #     stop('need to specify factor_use if plotting donor scores')
  #   }
  #   
  #   donor_scores <- container$tucker_results[[1]]
  #   
  #   # first limit cells in subclusts to those that we actually have scores for
  #   donor_vec <- container$scMinimal_full$metadata[names(subclusts),'donors']
  #   subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]
  #   con$clusters$leiden$groups <- as.factor(subclusts)
  #   con[["embedding"]] <- con[["embedding"]][names(subclusts),]
  #   
  #   # now get the donor_scores for the donors of these cells
  #   donor_vec <- donor_vec[donor_vec %in% rownames(donor_scores)]
  #   d_scores <- donor_scores[as.character(donor_vec),factor_use]
  #   
  #   ### testing multiplying scores time proportions
  #   subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
  #   scMinimal <- container$scMinimal_ctype[[ctype]]
  #   sub_meta_tmp <- scMinimal$metadata[names(subclusts),]
  #   
  #   # get donor proportions of subclusters
  #   donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)
  #   
  #   get_indv_subtype_associations(container,donor_props)
  #   
  #   # donor_props_transform <- donor_props
  #   
  #   # props_means <- colMeans(donor_props)
  #   # donor_props_transform <- sweep(donor_props,MARGIN=2,props_means,'/')
  #   # donor_props_transform <- donor_props
  #   # min_max_trans <- function(mycol) {
  #   #   mymin <- min(mycol)
  #   #   mymax <- max(mycol)
  #   #   tfed <- c()
  #   #   for (i in mycol) {
  #   #     tmp <- (i - mymin) / (mymax - mymin)
  #   #     tfed <- c(tfed,tmp)
  #   #   }
  #   #   return(tfed)
  #   # }
  #   
  #   # for (j in 1:ncol(donor_props)) {
  #   #   cmin <- min(donor_props[,j])
  #   #   cmax <- max(donor_props[,j])
  #   #   for (i in 1:nrow(donor_props)) {
  #   #     donor_props[i,j] <- (donor_props[i,j] - cmin) / (cmax - cmin)
  #   #   }
  #   # }
  #   # donor_props_transform <- donor_props
  #   
  #   # one below calculates zsc
  #   for (j in 1:ncol(donor_props)) {
  #     cmean <- mean(donor_props[,j])
  #     csd <- sd(donor_props[,j])
  #     for (i in 1:nrow(donor_props)) {
  #       donor_props[i,j] <- (donor_props[i,j] - cmean) / csd
  #     }
  #   }
  #   donor_props_transform <- donor_props
  #   
  #   donor_props_vec <- sapply(1:length(d_scores), function(x) {
  #     d <- names(d_scores)[x]
  #     sub_num <- subclusts_num[x]
  #     return(donor_props_transform[d,sub_num])
  #   })
  #   d_scores <- d_scores * donor_props_vec
  #   ###
  #   
  #   names(d_scores) <- names(subclusts)
  #   
  #   subc_plot <- con$plotGraph(colors=d_scores)
  # } else {
  #   subc_plot <- con$plotGraph()
  # }
  
  # reset the embedding and clusters
  con$clusters$leiden$groups <- orig_clusts
  con[["embedding"]] <- orig_embed
  
  return(subc_plot)
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
    reg_stats_all[[paste0("ct_",j,"_")]] <- reg_stats
  }
  
  reg_stats_all <- unlist(reg_stats_all)
  reg_stats_all <- stats:::p.adjust(reg_stats_all, method = 'fdr')
  
  parsed_name <- sapply(names(reg_stats_all),function(x){
    return(as.numeric(strsplit(x,split="_")[[1]][3]))
  })
  print(parsed_name)
  reg_stats_all <- reg_stats_all[parsed_name==factor_select]
  
  return(reg_stats_all)
}


#' Replace leiden groups in con object with subclusters identified for one cell type
#'
#' @param con conos Object for the dataset with umap projection and groups as cell types
#' @param subclusts numeric Named vector of subcluster assignments for one major cell type
#' only. Names should be the cell barcodes.
#' @param sub_ctype character The name of the major cell type that subclusts belong to
#'
#' @return a factor of the new groups with subclusters included
replace_groups_with_subclusts <- function(con,subclusts,sub_ctype) {
  
  # add back the large group lable to subtype numbers
  subclusts <- sapply(subclusts,function(x) {
    return(paste0(sub_ctype,'_',as.character(x)))
  })
  
  # convert groups to char so can replace by name
  groups <- as.character(con$clusters$leiden$groups)
  names(groups) <- names(con$clusters$leiden$groups)
  groups[names(subclusts)] <- subclusts
  
  # convert back to factor for con
  groups <- as.factor(groups)
  
  return(groups)
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


plot_subclust_de_hmap <- function(container,ctypes,resolutions) {
  ncores <- container$experiment_params$ncores
  
  con <- container[["embedding"]]
  orig_groups <- con$clusters$leiden$groups
  
  for (i in 1:length(ctypes)) {
    ctype <- ctypes[i]
    res <- resolutions[i]
    resolution_name <- paste0('res:',as.character(res))
    subclusts <- container$subclusters[[ctype]][[resolution_name]]
    new_groups <- replace_groups_with_subclusts(con,subclusts,ctype)
    con$clusters$leiden$groups <- new_groups
  }
  
  de.info <- con$getDifferentialGenes(n.cores=ncores, append.auc=TRUE,
                                      z.threshold=0, upregulated.only=TRUE)

  myhmap <- plotDEheatmap_conos(con, groups=con$clusters$leiden$groups, de=de.info, n.genes.per.cluster = 5,
                row.label.font.size = 7)
  

  con$clusters$leiden$groups <- orig_groups

  container$subcluster_de_heatmap <- myhmap
  return(container)
}



plotDEheatmap_conos <- function(con, groups, de=NULL, n.genes.per.cluster=5,
                                expression.quantile=0.99,pal=colorRampPalette(c('dodgerblue1','grey95','indianred1'))(1024),
                                ordering='-AUC',column.metadata=NULL, remove.duplicates=TRUE, show.cluster.legend=TRUE,
                                show_heatmap_legend=FALSE, row.label.font.size=10, ...) {
  
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
  
  de <- lapply(de,function(x) x%>%dplyr::arrange(!!rlang::parse_expr(ordering))%>%head(n.genes.per.cluster))
  de <- de[unlist(lapply(de, nrow))>0]
  
  gns <- lapply(de,function(x) as.character(x$Gene)) %>% unlist
  sn <- function(x) setNames(x,x)
  expl <- lapply(de,function(d) do.call(rbind,lapply(sn(as.character(d$Gene)),function(gene) conos:::getGeneExpression(con,gene))))
  
  
  exp <- do.call(rbind,expl)
  # limit to cells that were participating in the de
  exp <- na.omit(exp[,colnames(exp) %in% names(na.omit(groups))])
  
  # transform expression values
  x <- t(apply(as.matrix(exp), 1, function(xp) {
    if(expression.quantile<1) {
      qs <- quantile(xp,c(1-expression.quantile,expression.quantile))
      if(diff(qs)==0) { # too much, set to adjacent values
        xps <- unique(xp)
        if(length(xps)<3) { qs <- range(xp) } # only two values, just take the extremes
        xpm <- median(xp)
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
  
  o <- order(groups[colnames(x)])
  x=x[,o]
  
  rownames(x) <- convert_gn(container,rownames(x))
  
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
  
  
  column.metadata.colors <- list();
  
  # make sure cluster colors are defined
  if(is.null(column.metadata.colors[['clusters']])) {
    uc <- unique(annot$clusters);
    column.metadata.colors$clusters <- setNames(rainbow(length(uc)),uc)
  }

  if(remove.duplicates) { x <- x[!duplicated(rownames(x)),] }
  
  # draw heatmap
  ha <- ComplexHeatmap::HeatmapAnnotation(df=annot,col=column.metadata.colors,show_legend=show.cluster.legend)
  
  ha <- ComplexHeatmap::Heatmap(x, name='expression', col=pal, cluster_rows=FALSE, cluster_columns=FALSE, show_column_names=FALSE, top_annotation=ha, row_names_gp = grid::gpar(fontsize = row.label.font.size));

  return(ha)
}

















