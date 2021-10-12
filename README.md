
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![<kharchenkolab>](https://circleci.com/gh/kharchenkolab/scITD.svg?style=svg)](https://app.circleci.com/pipelines/github/kharchenkolab/scITD)
<!-- badges: end -->

# scITD

Single-Cell Interpretable Tensor Decomposition (scITD) is computational method capable of
extracting multicellular gene expression programs from single-cell
RNA-sequencing data. These programs may be represented by patterns of genes which are relevant in one or more unique cell types. Specifically, our tool enables one to find such patterns that are most variable across donors in multi-donor single-cell dataset collections. Therefore, scITD has a wide range of potential applications, including the study of population-level inter-individual variation, patient sub-grouping/stratification, and interrogating sample-level batch effects. The multicellular information provided by our method allows one to gain a deeper understanding of the ways that cells might be interacting or responding cetain stimuli. To enable such insights, we also provide an integrated suite of downstream data processing tools to transform the scITD output into succinct, yet informative summaries of the data.

## Installation

Many of our visualizations require the development version of ComplexHeatmap package, so this should be installed first:
    
``` r
devtools::install_github("kharchenkolab/scITD")
```  

Then, install scITD with the following:

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
