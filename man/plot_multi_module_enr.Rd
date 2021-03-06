% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_LR_interact.R
\name{plot_multi_module_enr}
\alias{plot_multi_module_enr}
\title{Generate gene set x ct_module heatmap showing significantly enriched sets}
\usage{
plot_multi_module_enr(
  container,
  ctypes,
  modules,
  sig_thresh = 0.05,
  db_use = "TF"
)
}
\arguments{
\item{container}{environment Project container that stores sub-containers
for each cell type as well as results and plots from all analyses}

\item{ctypes}{character A vector of cell type names corresponding to the module
numbers in mod_select, specifying the modules to compute enrichment for}

\item{modules}{numeric A vector of module numbers corresponding to the cell
types in ctype, specifying the modules to compute enrichment for}

\item{sig_thresh}{numeric P-value threshold for results to include. Only shows
a given gene set if at least one module has a result lower than the threshold.
(default=0.05)}

\item{db_use}{character The database of gene sets to use. Database
options include "GO", "Reactome", "KEGG", "BioCarta", "Hallmark", "TF", and
"immuno". More than one database can be used. (default="GO")}
}
\value{
the heatmap plot of enrichment results
}
\description{
Generate gene set x ct_module heatmap showing significantly enriched sets
}
