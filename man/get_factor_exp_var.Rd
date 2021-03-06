% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_tucker_ica.R
\name{get_factor_exp_var}
\alias{get_factor_exp_var}
\title{Get explained variance for each cell type for one factor}
\usage{
get_factor_exp_var(container, factor_use)
}
\arguments{
\item{container}{environment Project container that stores sub-containers
for each cell type as well as results and plots from all analyses}

\item{factor_use}{numeric The factor to investigate}
}
\value{
explained variance for each cell type in a list
}
\description{
Get explained variance for each cell type for one factor
}
