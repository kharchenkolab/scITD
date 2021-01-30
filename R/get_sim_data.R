

#' Add batch effect to real dataset
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses (default=NULL)
#' @param nbatches numeric The number of batches to create. Donors are randomly
#' assigned to batches. (default=2)
#' @param batch_override list A list of vectors specifying donors in each batch
#' (default=NULL)
#'
#' @return either project container with batch added to normalized data or the counts
#' matrix with batch effect added to it
#' @export
get_hybrid_sim <- function(container,nbatches=2,batch_override=NULL) {

  # get function from splatter
  getLNormFactors <- utils::getFromNamespace("getLNormFactors", "splatter")

  ## can also later try only making some genes affected by batch
  nGenes <- nrow(container$scMinimal_full$count_data)
  batch_factors <- c()
  for (i in 1:nbatches) {
    # batch.facs <- getLNormFactors(nGenes, 1, 0.5, .1, .1)
    batch.facs <- getLNormFactors(nGenes, 1, 0.5, .5, .1)
    batch_factors <- cbind(batch_factors,batch.facs)
  }

  colnames(batch_factors) <- sapply(1:nbatches,function(x){
    paste0('Batch',as.character(x))
  })

  rownames(batch_factors) <- rownames(container$scMinimal_full$count_data)

  # assign cells to batches by random donor groups
  container$scMinimal_full$metadata$batch <- NULL
  if (is.null(batch_override)) {
    donors_rem <- as.character(unique(container$scMinimal_full$metadata$donors))
    num_d <- length(donors_rem)
    for (i in 1:(nbatches-1)) {
      d_assigns <- sample(donors_rem,round(num_d/nbatches))
      container$scMinimal_full$metadata[container$scMinimal_full$metadata$donors %in% d_assigns,'batch'] <- paste0('Batch',as.character(i))
      donors_rem <- donors_rem[!(donors_rem %in% d_assigns)]
    }
    d_assigns <- donors_rem
    container$scMinimal_full$metadata[container$scMinimal_full$metadata$donors %in% d_assigns,'batch'] <- paste0('Batch',as.character(i+1))
  } else {
    for (i in 1:nbatches) {
      d_assigns <- batch_override[[i]]
      container$scMinimal_full$metadata[container$scMinimal_full$metadata$donors %in% d_assigns,'batch'] <- paste0('Batch',as.character(i))
    }
  }


  # multiply expression by batch factor for assigned cells
  new_dat <- c()
  for (i in 1:nbatches) {
    batch_name <- paste0('Batch',as.character(i))
    cells_alter <- rownames(container$scMinimal_full$metadata)[container$scMinimal_full$metadata$batch==batch_name]
    data_alter <- container$scMinimal_full$count_data[,cells_alter]
    data_alter <- data_alter * batch_factors[,batch_name]
    new_dat <- cbind(new_dat,data_alter)
  }
  # reorder columns to match original dataframe
  new_dat <- new_dat[,colnames(container$scMinimal_full$count_data)]
  container$scMinimal_full$count_data <- new_dat

  return(container)

}


#' Add batch effect to real dataset
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses (default=NULL)
#' @param nbatches numeric The number of batches to create. Donors are randomly
#' assigned to batches. (default=2)
#'
#' @return either project container with batch added to normalized data or the counts
#' matrix with batch effect added to it
#' @export
get_hybrid_sim_nonlin <- function(container,nbatches=2) {

  # get function from splatter
  getLNormFactors <- utils::getFromNamespace("getLNormFactors", "splatter")

  ## can also later try only making some genes affected by batch
  nGenes <- ncol(container$scMinimal_ctype[[1]]$pseudobulk)
  batch_factors <- c()
  for (i in 1:nbatches) {
    # batch.facs <- getLNormFactors(nGenes, 1, 0.5, .1, .1)
    batch.facs <- getLNormFactors(nGenes, 1, 0.5, .35, .1)
    batch_factors <- cbind(batch_factors,batch.facs)
  }

  colnames(batch_factors) <- sapply(1:nbatches,function(x){
    paste0('Batch',as.character(x))
  })

  rownames(batch_factors) <- container$tensor_data[[2]]

  scMinimal <- container$scMinimal_ctype[[1]]
  meta <- scMinimal$metadata

  meta$batch <- NULL
  donors_rem <- as.character(unique(meta$donors))
  num_d <- length(donors_rem)
  for (i in 1:(nbatches-1)) {
    d_assigns <- sample(donors_rem,round(num_d/nbatches))
    meta[meta$donors %in% d_assigns,'batch'] <- paste0('Batch',as.character(i))
    donors_rem <- donors_rem[!(donors_rem %in% d_assigns)]
  }
  d_assigns <- donors_rem
  meta[meta$donors %in% d_assigns,'batch'] <- paste0('Batch',as.character(i+1))
  scMinimal$metadata <- meta
  meta_assigns <- meta


  for (i in 2:length(container$experiment_params$ctypes_use)) {
    ct <- container$experiment_params$ctypes_use[i]
    scMinimal <- container$scMinimal_ctype[[ct]]
    meta <- scMinimal$metadata
    meta$batch <- NULL
    meta$batch <- sapply(meta$donors,function(x) {
      meta_assigns[meta_assigns$donors==x,'batch'][1]
    })
    scMinimal$metadata <- meta
  }

  # multiply expression by batch factor for assigned cells
  for (ct in container$experiment_params$ctypes_use) {
    scMinimal <- container$scMinimal_ctype[[ct]]
    mydat <- scMinimal$pseudobulk
    new_dat <- c()
    for (i in 1:nbatches) {
      batch_name <- paste0('Batch',as.character(i))
      donors_alter <- unique(scMinimal$metadata$donors[scMinimal$metadata$batch==batch_name])
      data_alter <- scMinimal$pseudobulk[donors_alter,]
      data_alter <- t(data_alter) * batch_factors[,batch_name]
      new_dat <- cbind(new_dat,data_alter)
    }
    new_dat <- t(new_dat)
    # reorder columns to match original dataframe
    new_dat <- new_dat[rownames(mydat),colnames(mydat)]
    scMinimal$pseudobulk <- new_dat
  }

  # need to also add batch metadata to overall metadata matrix
  # need to make sure the full data is limited to the cells used in analysis
  all_cells <- c()
  for (ct in container$experiment_params$ctypes_use) {
    cells_in_ctype <- rownames(container$scMinimal_ctype[[ct]]$metadata)
    all_cells <- c(all_cells,cells_in_ctype)
  }

  container$scMinimal_full$metadata <- container$scMinimal_full$metadata[all_cells,]
  container$scMinimal_full$count_data <- container$scMinimal_full$count_data[,all_cells]
  container$scMinimal_full$metadata$donors <- factor(container$scMinimal_full$metadata$donors,
                                                     levels=unique(container$scMinimal_full$metadata$donors))
  container$scMinimal_full$metadata$batch <- sapply(container$scMinimal_full$metadata$donors,function(x) {
    meta_assigns[meta_assigns$donors==x,'batch'][1]
  })


  return(container)

}




