
#' Create a container to store all data and results for the project
#'
#' @param scMinimal environment A sub-container for the project typically
#' consisting of gene expression data in its raw and processed forms as well
#' as metadata
#' @param ctypes_use character Names of the cell types to use for the analysis
#' (default=NULL)
#' @param gn_convert data.frame Gene identifier -> gene name conversions table.
#' Gene identifiers used in counts matrices should appear in the first column and
#' the corresponding gene symbols should appear in the second column. Can remain
#' NULL if the identifiers are already gene symbols. (default=NULL)
#' @param scale_var logical TRUE to scale the gene expression variance across donors
#' for each cell type. If FALSE then all genes are scaled to unit variance across
#' donors for each cell type. (default=TRUE)
#' @param var_scale_power numeric Exponent of normalized variance that is
#' used for variance scaling. Variance for each gene
#' is initially set to unit variance across donors (for a given cell type).
#' Variance for each gene is then scaled by multiplying the unit scaled values
#' by each gene's normalized variance (where the effect of the mean-variance
#' dependence is taken into account) to the exponent specified here. (default=NULL)
#' @param rotate_modes character The names of the tensor modes to rotate with
#' ICA during Tucker decomposition. Can include 'donors', 'genes', and/or 'ctypes'
#' (default='donors')
#' @param tucker_type character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
#' with sparsity constraints (default='regular')
#' @param rotation_type character Set to 'ica' to perform ICA rotation on resulting donor factor
#' matrix and loadings. Otherwise set to 'varimax' to perform varimax rotation. (default='ica')
#' @param ranks numeric The number of donor, gene, and cell type ranks, respectively,
#' to decompose to using Tucker decomposition (default=NULL)
#' @param ncores numeric Number of cores to use (default=4)
#' @param rand_seed numeric Random seed to use for all stochastic analyses to make
#' results reproducible (default=10)
#'
#' @return project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @export
make_new_container <- function(scMinimal, ctypes_use=NULL, gn_convert=NULL, scale_var=TRUE,
                               var_scale_power=NULL, tucker_type='regular', rotation_type='ica',
                               rotate_modes='donors', ranks=NULL, ncores = 4, rand_seed=10) {
  
  container <- new.env()
  container$scMinimal_full <- scMinimal
  container$donors <- scMinimal$donors
  container$ctypes <- scMinimal$ctypes
  if (!is.null(gn_convert)) {
    rownames(gn_convert) <- gn_convert[,1]
  }
  container$gn_convert <- gn_convert
  container$experiment_params <- list(ctypes_use=ctypes_use,
                                      scale_var=scale_var,
                                      var_scale_power=var_scale_power,
                                      rotate_modes=rotate_modes,
                                      tucker_type=tucker_type,
                                      rotation_type=rotation_type,
                                      ranks=ranks, ncores=ncores,
                                      rand_seed=rand_seed)
  container$experiment_params$run_check <- FALSE
  container$all_vargenes <- c()
  container$scMinimal_ctype <- list()
  container$tensor_data <- NULL
  container$tucker_results <- NULL
  container$stability_results <- NULL
  container$gene_score_associations <- NULL
  container$plots <- NULL
  container$rank_determination_results <- NULL
  container$var_scale_results <- NULL

  return(container)
}


#' Get union of significantly variable genes from all cell types
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param vargenes character Significantly variable genes from one cell type
#'
#' @return the project container with updated slot container$all_vargenes
#' @export
set_vargenes_container <- function(container,vargenes) {
  vargenes <- unique(vargenes)

  # need to make sure vargenes union genes in all cell types
  for (ct in container$experiment_params$ctypes_use) {
    scMinimal <- container$scMinimal_ctype[[ct]]
    genes_in_ctype <- rownames(scMinimal$data_sparse)
    vargenes <- intersect(vargenes,genes_in_ctype)
  }

  container$all_vargenes <- vargenes

  return(container)
}

#' Add an environment for a cell type to the container
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param scMinimal environment A sub-container for the project with one
#' cell type's gene expression data in its raw and processed forms as well
#' as metadata
#'
#' @return the container with added scMinimal for a cell type in a slot under
#' container$scMinimal_ctypes
#' @export
add_ctype_data_to_container <- function(container,scMinimal) {
  if (length(scMinimal$ctypes) > 1) {
    # throw error because should only be one cell type
    print("object contains more than one cell type")
  } else {
    container$scMinimal_ctype[[scMinimal$ctypes]] <- scMinimal
  }
  return(container)
}

#' Set experiment parameters
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param ctypes_use character Names of the cell types to use for the analysis
#' (default=NULL)
#' @param scale_var logical TRUE to scale the gene expression variance across donors
#' for each cell type. If FALSE then all genes are scaled to unit variance across
#' donors for each cell type. (default=NULL)
#' @param ranks numeric The number of donor, gene, and cell type ranks, respectively,
#' to decompose to using Tucker decomposition (default=NULL)
#' @param rotate_modes character The names of the tensor modes to rotate with
#' ICA during Tucker decomposition. Can include 'donors', 'genes', and/or 'ctypes'
#' (default=NULL)
#' @param tucker_type character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
#' with sparsity constraints (default=NULL)
#' @param rotation_type character Set to 'ica' to perform ICA rotation on resulting donor factor
#' matrix and loadings. Otherwise set to 'varimax' to perform varimax rotation. (default=NULL)
#' @param var_scale_power numeric Exponent of normalized variance that is
#' used for variance scaling. Variance for each gene
#' is initially set to unit variance across donors (for a given cell type).
#' Variance for each gene is then scaled by multiplying the unit scaled values
#' by each gene's normalized variance (where the effect of the mean-variance
#' dependence is taken into account) to the exponent specified here. (default=NULL)
#' @param ncores numeric Number of cores to use (default=NULL)
#'
#' @return the project container with updated experiment parameters in
#' container$experiment_params
#' @export
set_experiment_params <- function(container, ctypes_use=NULL, scale_var=NULL,
                                  ranks=NULL, rotate_modes=NULL, tucker_type=NULL,
                                  rotation_type=NULL, var_scale_power=NULL, 
                                  ncores=NULL) {
  # if user/code enters a value for a param then reset its value
  if (!is.null(ctypes_use)) {
    container$experiment_params$ctypes_use <- ctypes_use
  }
  if (!is.null(scale_var)) {
    container$experiment_params$scale_var <- scale_var
  }
  if (!is.null(ranks)) {
    container$experiment_params$ranks <- ranks
  }
  if (!is.null(rotate_modes)) {
    container$experiment_params$rotate_modes <- rotate_modes
  }
  if (!is.null(tucker_type)) {
    container$experiment_params$tucker_type <- tucker_type
  }
  if (!is.null(rotation_type)) {
    container$experiment_params$rotation_type <- rotation_type
  }
  if (!is.null(var_scale_power)) {
    container$experiment_params$var_scale_power <- var_scale_power
  }
  if (!is.null(ncores)) {
    container$experiment_params$ncores <- ncores
  }
  return(container)
}


#' Extract metadata for sex information if not provided already
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#'
#' @return the scMinimal object with sex metadata added to the metadata
#' @export
identify_sex_metadata <- function(container) {
  scMinimal <- container$scMinimal_full

  dge_sparse <- t(scMinimal$data_sparse)

  # get donor means for each gene in the dataset
  donor_meta <- as.factor(scMinimal$metadata$donors)
  means_all <- get_means(dge_sparse,donor_meta,table(donor_meta))
  means_all <- means_all[2:nrow(means_all),]

  # convert rownames to gene symbols using provided mapping
  gn_names <- convert_gn(container, colnames(means_all))

  y_ndx <- which(gn_names == 'RPS4Y1')
  x_ndx <- which(gn_names == 'XIST')
  y_mean <- mean(means_all[,y_ndx])
  x_mean <- mean(means_all[,x_ndx])

  make_note <- FALSE
  scMinimal$metadata$sex <- NA
  for (i in 1:nrow(scMinimal$metadata)) {
    d <- scMinimal$metadata$donors[i]
    if (means_all[d,y_ndx] > y_mean && means_all[d,x_ndx] < x_mean) {
      scMinimal$metadata$sex[i] <- 'M'
    } else if (means_all[d,y_ndx] < y_mean && means_all[d,x_ndx] > x_mean) {
      scMinimal$metadata$sex[i] <- 'F'
    } else {
      scMinimal$metadata$sex[i] <- 'A'
      make_note <- TRUE
    }
  }

  if (make_note) {
    print('Some assignments are ambiguous and are labeled A in the metadata.
          We recommend correcting these manually or providing the sex metadata when
          instantiating scMinimal.')
  }

  container$scMinimal_full <- scMinimal

  return(container)

}










