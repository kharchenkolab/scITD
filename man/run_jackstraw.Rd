% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_jackstraw.R
\name{run_jackstraw}
\alias{run_jackstraw}
\title{Run tensor-based jackstraw to get gene_cell type combinations that are significantly
associated with donor scores for factors extracted by Tucker decomposition}
\usage{
run_jackstraw(
  container,
  ranks,
  n_fibers = 100,
  n_iter = 500,
  tucker_type = "regular",
  rotation_type = "ica"
)
}
\arguments{
\item{container}{environment Project container that stores sub-containers
for each cell type as well as results and plots from all analyses}

\item{ranks}{numeric The number of donor, gene, and cell type ranks, respectively,
to decompose to using Tucker decomposition.}

\item{n_fibers}{numeric The number of fibers the randomly shuffle in each iteration
(default=100)}

\item{n_iter}{numeric The number of shuffling iterations to complete (default=500)}

\item{tucker_type}{character Set to 'regular' to run regular tucker or to 'sparse' to run tucker
with sparsity constraints (default='regular')}

\item{rotation_type}{character Set to 'ica' to perform ICA rotation on resulting donor factor
matrix and loadings. Otherwise set to 'varimax' to perform varimax rotation. (default='ica')}
}
\value{
the project container with adjusted pvalues in container$gene_score_associations
}
\description{
Run tensor-based jackstraw to get gene_cell type combinations that are significantly
associated with donor scores for factors extracted by Tucker decomposition
}
\examples{
test_container <- run_jackstraw(test_container, ranks=c(2,4,2), n_fibers=6,
n_iter=50, tucker_type='regular', rotation_type='ica')
}
