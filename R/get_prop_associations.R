
# library(umap)
# library(gmodels)

# needs to have tucker results slot as not null so check that
get_prop_associations <- function(container) {
  umap_all <- container$scMinimal_full$umap

  # loop through cell types
  for (ct in container$experiment_params) {
    scMinimal <- container[["scMinimal_ctype"]][[ct]]

    # get umap coordinates
    umap_ctype <- umap_all[colnames(scMinimal$data_sparse),]

    # loop through increasing k values
    for (k in 1:10) {
      # compute associations for unshuffled cells
      # run kmeans
      cell_clusters <- stats::kmeans(umap_ctype,num_clust,nstart=20,iter.max=100,
                                     algorithm = "MacQueen")$cluster
      # get donor proportions of subclusters
      donor_props <- compute_donor_props(cell_clusters,scMinimal$metadata)
      fstats_real <- compute_associations(donor_props,container$tucker_results[[1]])
      # fstats_real <- compute_associations(container$tucker_results[[1]], scMinimal$metadata,
      #                                     umap_ctype, k)

      # compute associations for shuffled cells
      umap_shuffled <- umap_ctype
      rownames(umap_shuffled) <- sample(rownames(umap_shuffled))

      # run kmeans
      cell_clusters <- stats::kmeans(umap_shuffled,num_clust,nstart=20,iter.max=100,
                                     algorithm = "MacQueen")$cluster
      # get donor proportions of subclusters
      donor_props <- compute_donor_props(cell_clusters,scMinimal$metadata)
      fstats_shuffled <- compute_associations(donor_props,container$tucker_results[[1]])


      # store association results


    }
  }
}

get_all_associations <- function(container) {
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

  # compute associations
  sig_res <- compute_associations(donor_props,donor_scores)

  # plot results
  prop_figure <- plot_donor_props(donor_props,donor_scores,ctypes,sig_res)

  # save results
  container$plots$ctype_prop_factor_associations <- prop_figure

  return(container)
}

reduce_dimensions <- function(scMinimal) {
  data_use <- t(scMinimal$data_sparse)
  data_use <- scale(data_use, center = TRUE)
  data_reduced <- gmodels::fast.prcomp(data_use, scale. = FALSE, center = FALSE)
  data_reduced <- umap::umap(data_reduced$x[,1:20])$layout
  return(data_reduced)
}


compute_associations <- function(donor_props, donor_scores) {

  all_fstats <- c()
  # loop through factors
  for (f in 1:ncol(donor_scores)) {
    # compute association between donor proportions and donor scores
    tmp <- cbind(donor_scores[,f],donor_props[rownames(donor_scores),])
    colnames(tmp)[1] <- "dscore"

    # construct the model
    prop_model <- paste0("dscore ~ ", paste(colnames(donor_props),collapse=" + "))

    # run lm (need to figure out how to specify multiple explanatory vars)
    lmres <- stats::lm(prop_model, data=tmp)

    # extract f-statistic
    fstat <- summary(lmres)$fstatistic[[1]]
    all_fstats <- c(all_fstats,fstat)
  }
  return(all_fstats)
}


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
    props <- table(tmp_clusts)/length(tmp_clusts)
    names(props) <- sapply(names(props),function(x) {
      paste0('K',as.character(x))
    })
    for (j in 1:length(props)) {
      donor_props[d,names(props)[j]] <- props[j]
    }
  }
  return(donor_props)
}



plot_donor_props <- function(donor_props,donor_scores,ctype_mapping,significance) {
  all_plots <- list()
  # loop through factors
  for (f in 1:ncol(donor_scores)) {
    # compute association between donor proportions and donor scores
    tmp <- cbind(donor_scores[,f],donor_props[rownames(donor_scores),])
    colnames(tmp)[1] <- "dscore"

    # need to reshape the matrix
    tmp2 <- reshape2::melt(data=tmp, id.vars = 'dscore',
                           variable.name = 'ctypes', value.name = 'donor_proportion')

    if (!is.null(ctype_mapping)) {
      tmp2$ctypes <- sapply(tmp2$ctypes,function(x){
        return(ctype_mapping[x])
      })
    }

    # NEED TO CHANGE BACK CTYPE NAMES HERE... USE MAPPING (PROBABLY NEED PASS IT IN)
    p <- ggplot2::ggplot(tmp2, aes(x=dscore,y=donor_proportion,color=ctypes)) +
      geom_line() +
      ggtitle(paste0("Factor ",as.character(f))) +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlab("Donor Factor Score") +
      ylab("Cell Type Proportion") +
      annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
               label=paste0('F-Statistic: ',round(significance[f],digits=2)))

    all_plots[[f]] <- p
  }

  p_total <- ggpubr::ggarrange(plotlist = all_plots,
                         ncol = 1)
}

plot_subclust_associations <- function() {

}















