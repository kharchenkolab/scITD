% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_gsea.R
\name{ht_clusters}
\alias{ht_clusters}
\title{Visualize the similarity matrix and the clustering. Adapted from simplifyEnrichment package.
https://github.com/jokergoo/simplifyEnrichment/blob/master/R/ht_clusters.R}
\usage{
ht_clusters(
  mat,
  cl,
  dend = NULL,
  col = c("white", "red"),
  draw_word_cloud = is_GO_id(rownames(mat)[1]) || !is.null(term),
  term = NULL,
  min_term = 5,
  order_by_size = FALSE,
  exclude_words = character(0),
  max_words = 10,
  word_cloud_grob_param = list(),
  fontsize_range = c(4, 16),
  column_title = NULL,
  ht_list = NULL,
  use_raster = TRUE,
  ...
)
}
\arguments{
\item{mat}{A similarity matrix.}

\item{cl}{Cluster labels inferred from the similarity matrix, e.g. from `cluster_terms` or `binary_cut`.}

\item{dend}{Used internally.}

\item{col}{A vector of colors that map from 0 to the 95^th percentile of the similarity values.}

\item{draw_word_cloud}{Whether to draw the word clouds.}

\item{term}{The full name or the description of the corresponding GO IDs.}

\item{min_term}{Minimal number of functional terms in a cluster. All the clusters
with size less than ``min_term`` are all merged into one separated cluster in the heatmap.}

\item{order_by_size}{Whether to reorder clusters by their sizes. The cluster
that is merged from small clusters (size < ``min_term``) is always put to the bottom of the heatmap.}

\item{exclude_words}{Words that are excluded in the word cloud.}

\item{max_words}{Maximal number of words visualized in the word cloud.}

\item{word_cloud_grob_param}{A list of graphic parameters passed to `word_cloud_grob`.}

\item{fontsize_range}{The range of the font size. The value should be a numeric vector with length two.
The minimal font size is mapped to word frequency value of 1 and the maximal font size is mapped
to the maximal word frequency. The font size interlopation is linear.}

\item{column_title}{Column title for the heatmap.}

\item{ht_list}{A list of additional heatmaps added to the left of the similarity heatmap.}

\item{use_raster}{Whether to write the heatmap as a raster image.}

\item{...}{other parameters}
}
\value{
A list containing a `ComplexHeatmap::HeatmapList-class` object and GO term ordering.
}
\description{
Visualize the similarity matrix and the clustering. Adapted from simplifyEnrichment package.
https://github.com/jokergoo/simplifyEnrichment/blob/master/R/ht_clusters.R
}
