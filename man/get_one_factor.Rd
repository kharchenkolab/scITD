% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_tucker_ica.R
\name{get_one_factor}
\alias{get_one_factor}
\title{Get the donor scores and loadings matrix for a single-factor}
\usage{
get_one_factor(container, factor_select)
}
\arguments{
\item{container}{environment Project container that stores sub-containers
for each cell type as well as results and plots from all analyses}

\item{factor_select}{numeric The number corresponding to the factor to extract}
}
\value{
A list with the first element as the donor scores and the second element
as the corresponding loadings matrix for one factor.
}
\description{
Get the donor scores and loadings matrix for a single-factor
}
\examples{
f1_res <- get_one_factor(test_container, factor_select=1)
}
