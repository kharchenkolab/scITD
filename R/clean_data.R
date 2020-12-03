
#' Clean data to remove genes only expressed in a few cells and donors
#' with very few cells
#'
#' @param scMinimal environment A sub-container for the project typically
#' consisting of gene expression data in its raw and processed forms
#' @param donor_min_cells numeric Minimum threshold for number of cells per
#' donor (default=5)
#' @param gene_min_cells numeric Minimum threshold for number of cells
#' with nonzero expression of a gene (default=5)
#'
#' @return an scMinimal environment with cleaned counts matrix, transformed
#' counts matrix, and metadata fields
#' @export
clean_data <- function(scMinimal, donor_min_cells=5, gene_min_cells=5) {
  donor_counts <- table(scMinimal$metadata$donors)
  donors_keep <- names(donor_counts)[donor_counts > donor_min_cells]
  gene_counts <- rowSums(scMinimal$data_sparse > 0)
  genes_keep <- names(gene_counts)[gene_counts > gene_min_cells]
  clean_scMinimal <- subset_scMinimal(scMinimal, donors_use = donors_keep, genes_use = genes_keep)
  return(clean_scMinimal)
}
































