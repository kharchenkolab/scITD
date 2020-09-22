
#' Convert gene names
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param genes character Vector of the gene identifiers to be converted to
#' gene symbols
#'
#' @return character vector of corresponding gene symbols
#' @export
convert_gn <- function(container, genes) {
  if (!is.null(container$gn_convert)) {
    genes <- container$gn_convert[genes,2]
  }
  return(genes)
}


