
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![<kharchenkolab>](https://circleci.com/gh/kharchenkolab/scITD.svg?style=svg)](https://circleci.com/gh/circleci/circleci-docs)
<!-- badges: end -->

# scITD

Single-Cell Interpretable Tensor Decomposition (scITD) is capable of
extracting multi-cell type expression programs from single-cell
RNA-sequencing data. This tool is useful for scRNA-seq datasets derived
from many donors. The multicellular processes can be interpreted with
standard methods such as gene set enrichment analysis (GSEA) and
ligand-receptor analysis. Overall, this package provides the basic tools
necessary to extracting multicellular expression processes from
scRNA-seq datasets that will be useful in understanding complex and
heterogeneous diseases.

## Installation

You can install scITD with the following:

``` r
devtools::install_github("kharchenkolab/scITD")
```

## Walkthrough

Follow the [walkthrough](http://pklab.med.harvard.edu/jonathan/) to
learn how to use scITD. The tutorial walks through the standard
processing pipeline and applies it to a dataset of PBMC’s from 45
healthy donors.

## Citation

If you find `scITD` useful for your publication, please cite:

    Jonathan Mitchel, Evan Biederstedt, and Peter Kharchenko (2020). scITD: Single-Cell
    Interpretable Tucker Decomposition. R package version 0.1.0.
