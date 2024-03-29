% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manage_scMinimal.R
\name{seurat_to_scMinimal}
\alias{seurat_to_scMinimal}
\title{Convert Seurat object to scMinimal object. Generally, this should be done through calling the
make_new_container() wrapper function.}
\usage{
seurat_to_scMinimal(seurat_obj, metadata_cols = NULL, metadata_col_nm = NULL)
}
\arguments{
\item{seurat_obj}{Seurat object that has been cleaned and includes the normalized,
log-transformed counts. The meta.data should include a column with the header
'sex' and values of 'M' or 'F' if available. The metadata should
also have a column with the header 'ctypes' with the corresponding names of
the cell types as well as a column with header 'donors' that contains
identifiers for each donor.}

\item{metadata_cols}{character The names of the metadata columns to use
(default=NULL)}

\item{metadata_col_nm}{character New names for the selected metadata columns
if wish to change their names. If NULL, then the preexisting column names are
used. (default=NULL)}
}
\value{
An scMinimal object holding counts and metadata for a project.
}
\description{
Convert Seurat object to scMinimal object. Generally, this should be done through calling the
make_new_container() wrapper function.
}
