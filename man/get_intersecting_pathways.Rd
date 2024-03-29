% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_gsea.R
\name{get_intersecting_pathways}
\alias{get_intersecting_pathways}
\title{Extract the intersection of gene sets which are enriched in two or more cell types for a factor}
\usage{
get_intersecting_pathways(
  container,
  factor_select,
  these_ctypes_only,
  up_down,
  thresh = 0.05
)
}
\arguments{
\item{container}{environment Project container that stores sub-containers
for each cell type as well as results and plots from all analyses}

\item{factor_select}{numeric The factor to investigate}

\item{these_ctypes_only}{character A vector of cell types for which to get gene sets that are
enriched in all of these and not in any other cell types}

\item{up_down}{character Set to "up" to get the gene sets for the positive loading genes. Set
to "down" to get the gene sets for the negative loadings genes.}

\item{thresh}{numeric Pvalue significance threshold for selecting enriched sets (default=0.05)}
}
\value{
A vector of the intersection of pathways that are significantly enriched in two
or more cell types for a factor.
}
\description{
Extract the intersection of gene sets which are enriched in two or more cell types for a factor
}
