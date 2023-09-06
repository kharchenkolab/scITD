## Upcoming

## [1.0.4] - 2023-09-06

### Changed
- Updated readme to include instructions for installing other dependencies
- Fixed error in LR analysis caused by genes with 0 expression

## [1.0.3] - 2023-01-23

### Added
- New function for users to easily get gene significance p-values for a factor
- Additional url for new LR analysis tutorial in the README

### Changed
- Enabled `compute_LR_interact` to run when data contains alternate gene symbols
- Included warning/catch in `plot_loadings_annot` if no genes are significant
- Fixed bug in `determine_ranks_tucker` with variance scaling in wrong order
- Updated README figure

## [1.0.2] - 2022-03-23

### Added
- Added `count_word()` function from simplifyEnrichment package as the remoevd it
- Included necessary suggests packages for `count_word()` 

### Changed
- Updated README to include our preprint citation and new figure

## [1.0.1] - 2022-01-28

### Added
- Optional parameters `min_gs_size` and `max_gs_size` to `run_gsea_one_factor()`
- New function project_new_data() to project a multicellular pattern onto new data
- Warning if cell-level metadata is included as donor-level metadata
- Check that cell type names in `ctypes_use` are spelled correctly

### Changed
- Checks for empty batches before running ComBat batch correction
- Fixed label issue with LR output when only have associations with 1 factor

## [1.0.0] - 2021-10-08

### Added
- Option to use custom list of genes in the tensor
- Check for matching cell names between count matrix and metadata

### Removed
- Removed `optimize_var_scale_power()` and helpers as they're no longer used
- Removed `compare_factors()` and helpers as they're no longer used

### Changed
- Fixed `anova` vargenes method for genes with zero expression across donors
- Updated `run_jackstraw()` and `get_min_sig_genes()` to work with new rotations
- Fixed bug in `plot_subclust_associations()`
- Enabled storage of results from `get_subtype_prop_associations()`
- Updated README to no longer require dev version of ComplexHeatmap
- Updated DESCRIPTION authors and description
