
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![<kharchenkolab>](https://circleci.com/gh/kharchenkolab/scITD.svg?style=svg)](https://app.circleci.com/pipelines/github/kharchenkolab/scITD)
[![CRAN status](https://www.r-pkg.org/badges/version/scITD)](https://cran.r-project.org/package=scITD)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/scITD)](https://cran.r-project.org/package=scITD)
<!-- badges: end -->

<img src="https://github.com/kharchenkolab/scITD/blob/develop/inst/scITD_logo.png" align="right" height="200">

# scITD

- [Introduction](#introduction)
- [Installation](#installation)
- [Walkthrough](#walkthrough)
- [Citation](#citation)

## Introduction

Single-Cell Interpretable Tensor Decomposition (scITD) is computational
method capable of extracting multicellular gene expression programs that
vary across donors or samples. The approach is premised on the idea that
higher-level biological processes often involve the coordinated actions
and interactions of multiple cell types. Given single-cell expression
data from multiple heterogenous samples, scITD aims to detect these
joint patterns of dysregulation impacting multiple cell types. This
method has a wide range of potential applications, including the study
of inter-individual variation at the population-level, patient
sub-grouping/stratification, and the analysis of sample-level batch
effects. The multicellular information provided by our method allows one
to gain a deeper understanding of the ways that cells might be
interacting or responding to certain stimuli. To enable such insights,
we also provide an integrated suite of downstream data processing tools
to transform the scITD output into succinct, yet informative summaries
of the data.

<img src="https://github.com/kharchenkolab/scITD/blob/develop/inst/scITD_overview_v2.jpg" align="center" height="275">


## Installation

The package has several dependencies from Bioconductor. To install these:
``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "edgeR", "sva", "Biobase"))
```

Then, to install scITD from CRAN:

``` r
install.packages('scITD')
```

To install the latest version of scITD from GitHub:

``` r
devtools::install_github("kharchenkolab/scITD")
```

## Walkthrough

Follow the [walkthrough](https://pklab.med.harvard.edu/jonathan/) to
learn how to use scITD. The tutorial introduces the standard processing
pipeline and applies it to a dataset of PBMC’s from 45 healthy donors.

We also created a tutorial for running [ligand-receptor analysis](https://pklab.med.harvard.edu/jonathan/LR_analysis.html). This uses the same dataset as the main walkthrough.

## Citation

If you find scITD useful for your publication, please cite:

[Jonathan Mitchel, M. Grace Gordon, Richard K. Perez, Evan Biederstedt, Raymund Bueno, Chun Jimmie Ye, Peter V. Kharchenko (2022). Tensor        decomposition reveals coordinated multicellular patterns of transcriptional variation that distinguish and stratify disease individuals.        bioRxiv 2022.](https://www.biorxiv.org/content/10.1101/2022.02.16.480703v1)
