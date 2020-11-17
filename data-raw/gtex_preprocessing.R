library(cmapR)
library(data.table)
library(stringr)

# To reproduce these preprocessing steps download the following files from the GTEX site:
# - GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct
# - GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
# - GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

# Place the downloaded files in a single directory and change the base_dir string below
# to the location where the files are stored

base_dir <- '/home/jmitchel/data/gtex/'

gtex_tpm <- parse_gctx(paste0(base_dir,'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'))
gtex_tpm <- gtex_tpm@mat

donor_meta <- as.data.frame(fread(paste0(base_dir,"GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")))
sample_meta <- as.data.frame(fread(paste0(base_dir,"GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")))

# add column of donor meta data
sample_meta$donors <- sapply(sample_meta$SAMPID, function(x) {
  strsplit(x,split='-')[[1]][2]
})

# # convert sample ids to have dots instead of dashes
# sample_meta$SAMPID <- sapply(sample_meta$SAMPID, function(x) {
#   str_replace_all(x,'-','.')
# })

# subset meta data to only have samples present in the expression matrix
sample_meta <- sample_meta[sample_meta$SAMPID %in% colnames(gtex_tpm),]

# get number of unique donors per tissue
ttypes <- unique(sample_meta$SMTSD)

tissue_counts <- list()
donors_each_tissue <- list()
for (tis in ttypes) {
  donors_for_tis <- unique(sample_meta$donors[sample_meta$SMTSD==tis])
  donors_each_tissue[[tis]] <- donors_for_tis
  num_donors <- length(donors_for_tis)
  tissue_counts[[tis]] <- num_donors
}

# get top n tissues with most donors
n <- 7
tissue_counts <- unlist(tissue_counts)
tissue_counts <- tissue_counts[order(tissue_counts,decreasing=TRUE)]
top_n_tissues <- tissue_counts[1:n]

print(top_n_tissues)

# see how many donors are in the intersection of all those tissues...
donors_in_all <- donors_each_tissue[[names(top_n_tissues)[1]]]
for (i in 2:length(top_n_tissues)) {
  donors_in_all <- intersect(donors_in_all,donors_each_tissue[[names(top_n_tissues)[i]]])
}

print(length(donors_in_all))


# limit the meta data matrix to just the donors at the intersection as well as the tissues in all...
sample_meta_sub <- sample_meta[(sample_meta$donors %in% donors_in_all),]
sample_meta_sub <- sample_meta_sub[(sample_meta_sub$SMTSD %in% names(top_n_tissues)),]

donor_meta$SUBJID <- sapply(donor_meta$SUBJID, function(x) {
  strsplit(x,split='-')[[1]][2]
})

donor_meta_sub <- donor_meta[(donor_meta$SUBJID %in% donors_in_all),]

# add donor meta data to sample meta data
sample_meta_sub$sex <- sapply(sample_meta_sub$donors, function(x) {
  donor_sex <- donor_meta_sub[donor_meta_sub$SUBJID == x, 'SEX']
  return(donor_sex)
})
sample_meta_sub$age <- sapply(sample_meta_sub$donors, function(x) {
  donor_age <- donor_meta_sub[donor_meta_sub$SUBJID == x, 'AGE']
  return(donor_age)
})
sample_meta_sub$dthhrdy <- sapply(sample_meta_sub$donors, function(x) {
  donor_dth <- donor_meta_sub[donor_meta_sub$SUBJID == x, 'DTHHRDY']
  return(donor_dth)
})

# remove unneeded variables for final meta data df to be used in the analysis
gtex_meta <- sample_meta_sub[c('SAMPID','SMTSD','donors','sex','age','dthhrdy')]
rownames(gtex_meta) <- gtex_meta$SAMPID
gtex_meta$SAMPID <- NULL
colnames(gtex_meta)[1] <- 'ctypes'
print(head(gtex_meta))

# subset the expression matrix to just these donors with intersecting samples
gtex_tpm_sub <- gtex_tpm[,rownames(gtex_meta)]

# seeing if biomaRt does better...
library(biomaRt)
feature.names <- rownames(gtex_tpm_sub)
feature.names.clean <- sapply(strsplit(as.character(feature.names),"\\."),"[[",1)
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
features.converted <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
      filters = 'ensembl_gene_id', values=feature.names.clean, mart=ensembl)

features.converted <- features.converted[features.converted$hgnc_symbol!="",]

features.converted <- features.converted[!duplicated(features.converted$ensembl_gene_id),]

converted_ordered <- sapply(1:length(feature.names.clean), function(i) {
  hgnc <- NA
  ensg <- feature.names.clean[i]
  if (ensg %in% features.converted$ensembl_gene_id) {
    hgnc <- features.converted[features.converted$ensembl_gene_id==ensg,2]
  }
  return(hgnc)
})

converted_ordered <- unlist(converted_ordered)
feature.names <- feature.names[!is.na(converted_ordered)]
converted_ordered <- converted_ordered[!is.na(converted_ordered)]

feature.names.final <- cbind(feature.names,unlist(converted_ordered))
# saveRDS(feature.names.final,file='/home/jmitchel/data/gtex/genes.rds')
# feature.names.final <- readRDS(file='/home/jmitchel/data/gtex/genes.rds')

# reduce number of genes in gtex counts to ones with gene symbols
gtex_tpm_sub <- gtex_tpm_sub[feature.names.final[,1],]

# log transform the tpm values
gtex_tpm_sub_transform <- log2(gtex_tpm_sub + 1)
# saveRDS(gtex_tpm_sub_transform,file='/home/jmitchel/data/gtex/7_tissues_counts.rds')
# gtex_tpm_sub_transform <- readRDS(file='/home/jmitchel/data/gtex/7_tissues_counts.rds')

# make dthhrdy score a factor
gtex_meta$dthhrdy <- as.factor(gtex_meta$dthhrdy)
# saveRDS(gtex_meta,file='/home/jmitchel/data/gtex/7_tissues_meta.rds')
# gtex_meta <- readRDS(file='/home/jmitchel/data/gtex/7_tissues_meta.rds')























