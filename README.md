
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->
[![Build Status](https://travis-ci.com/kharchenkolab/scITD.svg?branch=master)](https://travis-ci.com/github/kharchenkolab/scITD)
<!-- badges: end -->

# scITD

Single-Cell Interpretable Tensor Decomposition (scITD) employs the
Tucker tensor decomposition to extract multi-cell type expression
programs from single-cell RNA-sequencing data. This tool is most useful
for scRNA-seq datasets derived from many source donors. By rearranging
the factor matrices and core tensor traditionally extracted with Tucker
we create a biologically meaningful result that can be interpreted with
standard methods such as gene set enrichment analysis (GSEA). To further
improve interpretability of the extracted factors we apply independent
component analysis (ICA) to rotate factors toward independence. We also
implement several methods to determine the appropriate ranks to
decompose the tensor to, as this can be one of the main challenges in
using such a method. Further, a jackstraw-like method has been
implemented to identify genes that are significant in each of the
extracted factors. Overall, this package provides the basic tools
necessary to extracting multi-cell type processes from scRNA-seq
datasets that will be useful in understanding complex and heterogeneous
diseases.

## Installation

You can install scITD with the following function call:

``` r
devtools::install_github("kharchenkolab/scITD")
```

## Walkthrough

Follow the [guided
tutorial](https://htmlpreview.github.io/?https://raw.githubusercontent.com/kharchenkolab/scITD//master/doc/walkthrough.html)
to learn how to use scITD. The tutorial walks through the standard
pipeline and applies it to a dataset of PBMC’s from 45 healthy donors.
The dataset is originally from a paper by [van der Wijst et
al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5905669/).
