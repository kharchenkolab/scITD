#' @import Matrix
NULL

#' Partition main gene by cell matrix into per cell type matrices with significantly
#' variable genes only
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param method character The method used to select significantly variable
#' genes across donors within a cell type (default="anova")
#' @param thresh numeric Pvalue threshold to use for gene significance
#' (default=0.01)
#'
#' @return the project container with scMinimal environments added for each
#' cell type. The scMinimal environments contain expression count matrices
#' with columns as the cells belonging to the respective cell type and rows
#' as the genes which were identified as significantly variable across donors
#' in any cell type. This is in the "count_data_sparse" slot. The
#' normalized, log-transformed count matrices with the same dimensions are
#' also available in the "data_sparse slot". Metadata is also present
#' in these sub-containers for the corresponding included cells and is in the
#' "metadata" slot.
#' @export
get_ctype_vargenes <- function(container, method="anova", thresh=0.01) {

  # get ctype matrices
  container <- get_ctype_data(container)

  # set random seed to work with mclapply
  RNGkind("L'Ecuyer-CMRG")
  set.seed(container$experiment_params$rand_seed)

  ncores <- container$experiment_params$ncores

  var_res <- data.frame(matrix(ncol=3, nrow=0))

  # check that ctypes_use param has been set
  if (is.null(container$experiment_params$ctypes_use)) {
    stop("ctypes_use parameter from container$experiment_params is NULL. Use set_experiment_params()")
  }

  for (ct in container$experiment_params$ctypes_use) {
    print(ct)
    if (method == "anova") {
      pvals <- vargenes_anova(container$scMinimal_ctype[[ct]], ncores)
    } else if (method == "empir") {
      pvals <- vargenes_shuffle(container$scMinimal_ctype[[ct]], 10000, ncores)
    }
    var_res <- rbind(var_res,cbind(names(pvals), pvals, rep(ct,length(pvals))))
  }
  colnames(var_res) <- c('genes', 'pvalues', 'ctypes')
  var_res <- as.data.frame(var_res)
  var_res$pvalues <- as.numeric(as.character(var_res$pvalues))
  var_res$genes <- as.character(var_res$genes)
  var_res$ctypes <- as.character(var_res$ctypes)
  var_res$padj <- p.adjust(var_res$pvalues, "fdr")
  vargenes_all <- var_res$genes[var_res$padj < thresh]
  container <- set_vargenes_container(container, vargenes_all)

  for (ct in container$experiment_params$ctypes_use) {
    ctype_res <- var_res[var_res$ctypes == ct,]
    ctype_vargenes <- ctype_res[ctype_res$padj < thresh, 'genes']
    container$scMinimal_ctype[[ct]]$vargenes <- ctype_vargenes
  }

  # reduce ctype data to only significantly variable genes
  container <- reduce_to_vargenes(container)

  return(container)
}

#' Create an environment for each cell type
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param make_clean logical TRUE to apply minimum thresholds for number of cells
#' expressing a gene and number of cells per donor (default=TRUE)
#'
#' @return the project container with the scMinimal environments added into
#' the container$scMinimal_ctypes slot
#' @export
get_ctype_data <- function(container, make_clean=TRUE) {
  for (ct in container$experiment_params$ctypes_use) {
    ctype_sub <- subset_scMinimal(container$scMinimal_full, ctypes_use=ct)
    if (make_clean) {
      ctype_sub <- clean_data(ctype_sub)
    }
    container <- add_ctype_data_to_container(container, ctype_sub)
  }

  ## need to ensure all ctype matrices have same genes...

  # first get intersection of genes in all ctypes not removed by cleaning
  ct1 <- container$experiment_params$ctypes_use[1]
  g1 <- rownames(container$scMinimal_ctype[[ct]]$data_sparse)
  genes_in_all <- g1
  for (ct in container$experiment_params$ctypes_use) {
    ctype_genes <- rownames(container$scMinimal_ctype[[ct]]$data_sparse)
    genes_in_all <- intersect(genes_in_all,ctype_genes)
  }

  # now only keep the intersection of genes
  for (ct in container$experiment_params$ctypes_use) {
    scMinimal <- container$scMinimal_ctype[[ct]]
    scMinimal <- subset_scMinimal(scMinimal, make_copy=FALSE, genes_use=genes_in_all)
  }

  return(container)
}

#' Reduce each cell type's expression matrix to just the significantly variable genes
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#'
#' @return the project container with gene by cell matrices in the scMinimal_ctype
#' slots reduced to contain only the genes that were found to be significantly
#' variable across donors in at least one cell type
#' @export
reduce_to_vargenes <- function(container) {
  vargenes <- container$all_vargenes
  for (ct in container$experiment_params$ctypes_use) {
    scMinimal <- container$scMinimal_ctype[[ct]]
    scMinimal <- subset_scMinimal(scMinimal, make_copy=FALSE, genes_use=vargenes)
  }
  return(container)
}
