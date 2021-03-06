% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_factors.R
\name{get_compare_go_enrich}
\alias{get_compare_go_enrich}
\title{Compute enriched gene sets for the previously run factor comparison}
\usage{
get_compare_go_enrich(container, ctype, direc, db_use = "GO")
}
\arguments{
\item{container}{environment Project container that stores sub-containers
for each cell type as well as results and plots from all analyses}

\item{ctype}{character The cell type to compute enrichment for}

\item{direc}{numeric Can either be set to 1 or -1. If compare_type was "same",
use 1. If compare_type was "different", then use 1 or -1 to choose which of the
two compared factors to compute enrichment for.}

\item{db_use}{character The database of gene sets to use. Database
options include "GO", "Reactome", "KEGG", and "BioCarta". More than
one database can be used. (default="GO")}
}
\value{
a named vector of adjusted p-values for the gene sets tested
}
\description{
Compute enriched gene sets for the previously run factor comparison
}
