## Upcoming

## [1.0.0] - 2021-10-08

### Added
- Option to use custom list of genes in the tensor
- Check for matching cell names between count matrix and metadata

### Removed
- Removed `optimize_var_scale_power()` and helpers as they're no longer used

### Changed
- Updated `run_jackstraw()` and `get_min_sig_genes()` to work with new rotations
- Fixed bug in `plot_subclust_associations()`
- Enabled storage of results from `get_subtype_prop_associations()`
- Updated README to no longer require dev version of ComplexHeatmap
- Updated DESCRIPTION authors and description
