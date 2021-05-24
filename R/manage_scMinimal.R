
#' Create an scMinimal object
#'
#' @param count_data sparseMatrix Matrix of raw counts with genes as rows
#' and cells as columns
#' @param meta_data data.frame Metadata with cells as rows and variables
#' as columns. Number of rows in metadata should equal number of columns
#' in count matrix.
#' @param metadata_cols character The names of the metadata columns to use
#' (default=NULL)
#' @param metadata_col_nm character New names for the selected metadata columns
#' if wish to change their names. If NULL, then the preexisting column names are
#' used. (default=NULL)
#'
#' @return scMinimal object
#' @export
instantiate_scMinimal <- function(count_data, meta_data, metadata_cols=NULL, metadata_col_nm=NULL) {
  scMinimal <- new.env()
  scMinimal$count_data <- count_data

  # subset and rename metadata columns if specified
  if (!is.null(metadata_cols)) {
    meta_data <- meta_data[,metadata_cols]
    if (!is.null(metadata_col_nm)) {
      colnames(meta_data) <- metadata_col_nm
    }
  }
  scMinimal$metadata <- meta_data

  # check if metadata has minimum required columns
  if (sum(c("donors", "ctypes") %in% colnames(meta_data)) != 2) {
    stop("Metadata does not have columns labeled 'donors' and 'ctypes', which is required.")
  }

  return(scMinimal)
}

#' Convert Seurat object to scMinimal object
#'
#' @param seurat_obj Seurat object that has been cleaned and includes the normalized,
#' log-transformed counts. The meta.data should include a column with the header
#' 'sex' and values of 'M' or 'F' if available. The metadata should
#' also have a column with the header 'ctypes' with the corresponding names of
#' the cell types as well as a column with header 'donors' that contains
#' identifiers for each donor.
#' @param metadata_cols character The names of the metadata columns to use
#' (default=NULL)
#' @param metadata_col_nm character New names for the selected metadata columns
#' if wish to change their names. If NULL, then the preexisting column names are
#' used. (default=NULL)
#'
#' @return scMinimal object
#' @export
seurat_to_scMinimal <- function(seurat_obj, metadata_cols=NULL, metadata_col_nm=NULL) {

  scMinimal <- new.env()
  scMinimal$count_data <- seurat_obj@assays$RNA@counts
  meta_data <- seurat_obj@meta.data

  # subset and rename metadata columns if specified
  if (!is.null(metadata_cols)) {
    meta_data <- meta_data[,metadata_cols]
    if (!is.null(metadata_col_nm)) {
      colnames(meta_data) <- metadata_col_nm
    }
  }
  scMinimal$metadata <- meta_data

  # check if metadata has minimum required columns
  if (sum(c("donors", "ctypes") %in% colnames(meta_data)) != 2) {
    stop("Metadata does not have columns labeled donors and ctypes, which is required.")
  }

  return(scMinimal)
}

#' Subset an scMinimal object by specified genes, donors, cells, or cell types
#'
#' @param scMinimal environment A sub-container for the project typically
#' consisting of gene expression data in its raw and processed forms as well
#' as metadata
#' @param ctypes_use character The cell types to keep (default=NULL)
#' @param cells_use character Cell barcodes for the cells to keep (default=NULL)
#' @param donors_use character The donors to keep (default=NULL)
#' @param genes_use character The genes to keep (default=NULL)
#' @param in_place logical If set to TRUE then replaces the input object with the
#' new subsetted object (default=TRUE)
#'
#' @return a subsetted scMinimal object
#' @export
subset_scMinimal <- function(scMinimal, ctypes_use=NULL,cells_use=NULL,
                             donors_use=NULL, genes_use=NULL, in_place=TRUE) {
  count_data <- scMinimal$count_data
  metadata <- scMinimal$metadata
  if (!is.null(ctypes_use)) {
    count_data <- count_data[,which(metadata$ctypes %in% ctypes_use)]
    metadata <- metadata[which(metadata$ctypes %in% ctypes_use),]
  } else if (!is.null(cells_use)) {
    count_data <- count_data[,cells_use]
    metadata <- metadata[cells_use,]
  } else if (!is.null(donors_use)) {
    count_data <- count_data[,which(metadata$donors %in% donors_use)]
    metadata <- metadata[which(metadata$donors %in% donors_use),]
  }

  if (!is.null(genes_use)) {
    # only use specified genes that are in the df
    genes_use <- genes_use[genes_use %in% rownames(count_data)]
    count_data <- count_data[genes_use,]
  }

  metadata$donors <- factor(metadata$donors,levels=as.character(unique(metadata$donors)))
  metadata$ctypes <- factor(metadata$ctypes,levels=as.character(unique(metadata$ctypes)))

  if (in_place) {
    scMinimal$metadata <- metadata
    scMinimal$count_data <- count_data
    scMinimal_sub <- scMinimal
  } else {
    scMinimal_sub <- instantiate_scMinimal(count_data,metadata)
  }

  return(scMinimal_sub)
}


































