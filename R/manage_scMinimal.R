
#' Create an scMinimal object
#'
#' @param data_sparse sparseMatrix Matrix of normalized, log-transformed
#' counts with genes as rows and cells as columns
#' @param count_data sparseMatrix Matrix of raw counts with genes as rows
#' and cells as columns
#' @param meta_data data.frame Metadata with cells as rows and variables
#' as columns. Number of rows in metadata should equal number of columns
#' in data_sparse
#'
#' @return scMinimal object
#' @export
instantiate_scMinimal <- function(data_sparse,count_data,meta_data) {
  scMinimal <- new.env()
  scMinimal$data_sparse <- data_sparse
  scMinimal$count_data_sparse <- count_data
  scMinimal$data_residuals <- NULL
  scMinimal$data_means <- NULL
  scMinimal$metadata <- meta_data
  scMinimal$vargenes <- c()
  scMinimal$associated_genes <- c()
  scMinimal$donors <- levels(scMinimal$metadata$donors)
  scMinimal$ctypes <- levels(scMinimal$metadata$ctypes)
  return(scMinimal)
}

#' Convert Seurat object to scMinimal object
#'
#' @param s_obj Seurat object that has been cleaned and includes the counts
#' as well as normalized, log-transformed counts. The meta.data should include
#' a column with the header 'sex' and values of 'M' or 'F'. The metadata should
#' also have a column with the header 'ctypes' with the corresponding names of
#' the cell types as well as a column with header 'donors' that contains
#' identifiers for each donor.
#'
#' @return scMinimal object
#' @export
seurat_to_scMinimal <- function(s_obj) {
  # check if metadata has minimum required columns
  if (sum(c("sex", "donors", "ctypes") %in% colnames(s_obj@meta.data)) != 3) {
    stop("Metadata does not have columns labeled sex, donors, and ctypes")
  }
  scMinimal <- new.env()
  scMinimal$data_sparse <- methods::as(as.matrix(Seurat::GetAssayData(s_obj)),'sparseMatrix')
  scMinimal$count_data_sparse <- methods::as(as.matrix(s_obj@assays$RNA@counts),'sparseMatrix')
  scMinimal$data_residuals <- NULL
  scMinimal$data_means <- NULL
  scMinimal$metadata <- s_obj@meta.data
  scMinimal$donors <- levels(s_obj@meta.data$donors)
  scMinimal$ctypes <- levels(s_obj@meta.data$ctypes)
  return(scMinimal)
}

#' Subset an scMinimal object by specified genes, donors, cells, or cell types
#'
#' @param scMinimal environment A sub-container for the project typically
#' consisting of gene expression data in its raw and processed forms as well
#' as metadata
#' @param make_copy logical If TRUE, puts the subsetted scMinimal in a new object.
#' If FALSE, overwrites the input scMinimal object with the subsetted one. (default=TRUE)
#' @param ctypes_use character The cell types to keep (default=NULL)
#' @param cells_use character Cell barcodes for the cells to keep (default=NULL)
#' @param donors_use character The donors to keep (default=NULL)
#' @param genes_use character The genes to keep (default=NULL)
#'
#' @return a subsetted scMinimal object
#' @export
subset_scMinimal <- function(scMinimal, make_copy=T, ctypes_use=NULL,
                             cells_use=NULL, donors_use=NULL, genes_use=NULL) {
  data_sparse <- scMinimal$data_sparse
  count_data_sparse <- scMinimal$count_data_sparse
  metadata <- scMinimal$metadata
  if (!is.null(ctypes_use)) {
    data_sparse <- data_sparse[,which(metadata$ctypes %in% ctypes_use)]
    count_data_sparse <- count_data_sparse[,which(metadata$ctypes %in% ctypes_use)]
    metadata <- metadata[which(metadata$ctypes %in% ctypes_use),]
  } else if (!is.null(cells_use)) {
    data_sparse <- data_sparse[,cells_use]
    count_data_sparse <- count_data_sparse[,cells_use]
    metadata <- metadata[cells_use,]
  } else if (!is.null(donors_use)) {
    data_sparse <- data_sparse[,which(metadata$donors %in% donors_use)]
    count_data_sparse <- count_data_sparse[,which(metadata$donors %in% donors_use)]
    metadata <- metadata[which(metadata$donors %in% donors_use),]
  }

  if (!is.null(genes_use)) {
    # only use specified genes that are in the df
    genes_use <- genes_use[genes_use %in% rownames(data_sparse)]
    data_sparse <- data_sparse[genes_use,]
    count_data_sparse <- count_data_sparse[genes_use,]
  }

  metadata$donors <- factor(metadata$donors,levels=as.character(unique(metadata$donors)))
  metadata$ctypes <- factor(metadata$ctypes,levels=as.character(unique(metadata$ctypes)))

  if (make_copy) {
    scMinimal_sub <- instantiate_scMinimal(data_sparse,count_data_sparse,metadata)
  } else {
    scMinimal$metadata <- metadata
    scMinimal$data_sparse <- data_sparse
    scMinimal$count_data_sparse <- count_data_sparse
    scMinimal_sub <- scMinimal
  }
  return(scMinimal_sub)
}







































