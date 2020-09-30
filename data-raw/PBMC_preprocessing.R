library(Seurat)
library(Matrix)
library(ggplot2)

# add rows/columns and concatenate all matrices
full_mat <- NULL # to store full dataframe
for (i in 1:8) {
  matrix_dir <- paste0(file.path(find.package('scITD'),'extdata','count_matrices_per_lane'),
                       "/lane_", as.character(i))

  barcode.path <- paste0(matrix_dir, "/barcodes.tsv")
  features.path <- paste0(matrix_dir, "/genes.tsv")
  matrix.path <- paste0(matrix_dir, "/matrix.mtx")

  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)

  # append lane name to barcodes
  barcode.names <- lapply(barcode.names[,1],function(x){
    tmp <- strsplit(x,"-")[[1]][1]
    tmp <- paste0(tmp,"_lane",as.character(i))
    return(tmp)
  })
  colnames(mat) = barcode.names
  rownames(mat) = feature.names$V1

  #append to total df
  if (is.null(full_mat)) {
    full_mat <- mat
  } else {
    full_mat <- cbind(full_mat,mat)
  }
}

# load the cell type annotations for cells that pass QC
ctypes <- read.table(file = file.path(find.package('scITD'),'extdata','barcodes_to_cell_types','barcodes_to_cell_types.tsv'),
           sep = '\t', header = TRUE, stringsAsFactors = F)

# remove cells that were not annotated
mask <- colnames(full_mat) %in% ctypes$cell_barcode
qc_mat <- full_mat[,mask]

# create seurat object
pbmc <- CreateSeuratObject(counts = qc_mat)

# load donor meta data and add it to seurat object
donors <- read.table(file = file.path(find.package('scITD'),'extdata','cell_barcodes','pilot3_persons.tsv'),
                     sep = '\t', header = FALSE, stringsAsFactors = F)
pbmc@meta.data$donors <- donors$V2

# add lane meta data to seurat object
lanes <- lapply(names(Idents(pbmc)),function(x){
  tmp <- strsplit(x,"_")[[1]][2]
  return(tmp)
})
pbmc@meta.data$lanes <- as.factor(unlist(lanes))

# add cell type annotation to metadata
pbmc@meta.data$ctypes <- as.factor(ctypes$cell_type)

# normalize data
pbmc <- NormalizeData(pbmc)

# remove two donors with only one cell
Idents(pbmc) <- pbmc@meta.data$donors
pbmc <- subset(pbmc, idents = c("1_LLDeep_0203","1_LLDeep_0198"), invert=TRUE)

# get median number cells per donor at fine resolution
Idents(pbmc) <- pbmc@meta.data$ctypes
ctypes <- levels(Idents(pbmc))
for (ctype in ctypes) {
  pbmc_sub <- subset(pbmc, idents = ctype)
  Idents(pbmc_sub) <- pbmc_sub@meta.data$donors
  print(ctype)
  print(table(Idents(pbmc_sub)))
  print(median(table(Idents(pbmc_sub))))
}

# remove "CD56(bright) NK", "mDC", "Megakaryocyte", "ncMonocyte", "pDC", and "Plasma"
# as they all have less than median of 10 cells per individual
pbmc <- subset(pbmc,idents = c("CD56(bright) NK", "mDC",
                               "Megakaryocyte", "ncMonocyte",
                               "pDC", "Plasma"), invert = TRUE)

# make donor metadata a factor
pbmc@meta.data$donors <- factor(pbmc@meta.data$donors,levels=as.character(unique(pbmc@meta.data$donors)))

# reduce dataset size for vignette
gene_counts <- rowSums(as.matrix(GetAssayData(pbmc)))
pbmc_sub <- subset(pbmc,features = rownames(pbmc)[gene_counts>10])
Idents(pbmc_sub) <- pbmc_sub@meta.data$donors
pbmc_sub <- subset(pbmc_sub,idents = c('s1','s2','s5','s6','s7','s9','s10','s11','s13','s14','s15',
                                       's18','s20','s23','s25','s26','s29','s30','s33','s35','s38',
                                       's40','s43','s44','s45')) #selecting ~half the donors (most variable ones)
Idents(pbmc_sub) <- pbmc_sub@meta.data$ctypes
pbmc_sub <- subset(pbmc_sub,idents = c('CD4+ T','cMonocyte','CD56(dim) NK')) #selecting just a few interesting cell types with most cells
pbmc_sub@meta.data$orig.ident <- NULL
pbmc_sub@meta.data$nCount_RNA <- NULL
pbmc_sub@meta.data$nFeature_RNA <- NULL
pbmc_sub_transformed <- methods::as(as.matrix(Seurat::GetAssayData(pbmc_sub)),'sparseMatrix')
pbmc_sub_meta <- pbmc_sub@meta.data
pbmc_sub_counts <- methods::as(as.matrix(pbmc_sub@assays$RNA@counts),'sparseMatrix')

# # save pbmc_sub_transformed, pbmc_sub_meta, and pbmc_sub_counts in data as .RData files for vignette (compress = TRUE)
# save(pbmc_sub_transformed,file='/home/jmitchel/scITD/data/pbmc_sub_transformed.RData',compress = "xz")
# save(pbmc_sub_meta,file='/home/jmitchel/scITD/data/pbmc_sub_meta.RData',compress = "xz")
# save(pbmc_sub_counts,file='/home/jmitchel/scITD/data/pbmc_sub_counts.RData',compress = "xz")















