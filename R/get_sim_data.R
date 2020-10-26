
#' Generate simulated data with subgroups of donors expressing multi-cell type processes
#' @import splatter
#' @importFrom SingleCellExperiment colData rowData counts
#'
#' @param data_for_p_est matrix Gene by cell matrix of counts. It is recommended to input a
#' small subsampled matrix for faster runtime. Need an input for either this parameter or the
#' params parameter. (default=NULL)
#' @param prev_params SplatParams Precomputed parameters from Splatter (default=NULL)
#' @param donors_total numeric Total number of donors to include in the simulation
#' @param cells_per_donor numeric Average number of cells per donor to generate. Total
#' dataset size will be cells_per_donor * donors_total number of cells.
#' @param n_processes numeric The number of multi-cell type processes/factors to generate.
#' Each process will be overexpressed in a subgroup of donors.
#' @param donors_per_process numeric The size of each subgroup of donors to go along with
#' each multi-cell type process. donors_per_process * n_processes should not exceed
#' donors_total.
#' @param n_ctypes numeric The number of different cell types to be generated (default=2)
#' @param n_genes numeric Number of genes to include in dataset. Should be considered alongside
#' de_prob, as too few DE genes will make it hard to detect a signal. (default=2000)
#' @param de_prob numeric The probability that each gene gets upregulated in each cell
#' type in each process (default=.05)
#' @param de_strength numeric This is roughly the average effect size for all differentially
#' expressed genes used to generate the distinc processes (default=2)
#' @param factor_overlap logical Set to TRUE to have ~20% of the DE genes be found across
#' multiple factors (default=TRUE)
#' @param rseed numeric The random seed to use (default=10)
#'
#' @return Normalized gene by cell count matrix and metadata data.frame. The DE genes can be
#' identified by the rownames of the normalized count matrix. GroupX_GeneY indicates that the
#' Gene Y was made to be upregulated in GroupX. The group numbers are ordered such that they
#' describe cell types 1-N for factor 1, then cell types 1-N for factor 2, and so on.
#' @export
get_sim_data <- function(data_for_p_est=NULL, prev_params=NULL, donors_total, cells_per_donor, n_processes, donors_per_process,
                         n_ctypes=2, n_genes=2000, de_prob=0.05, de_strength=2, factor_overlap=TRUE, rseed=10) {

  if (!is.null(data_for_p_est)) {
    params <- splatEstimate(data_for_p_est)
  } else if (!is.null(prev_params)) {
    params <- prev_params
  } else {
    params <- newSplatParams()
  }


  n_cells <- donors_total * cells_per_donor

  params <- setParam(params, "batchCells", n_cells)
  params <- setParam(params, "nGenes", n_genes)
  params <- setParam(params, "seed", rseed)

  # get general non-DE genes
  scsim <- splatSimulate(params)
  meta <- colData(scsim)
  sim_counts <- counts(scsim)
  lib_sizes <- colSums(sim_counts)

  # normalize sim counts
  sim_counts <- sweep(sim_counts,MARGIN=2,lib_sizes,'/')

  n_de_genes_per_group <- round(de_prob * n_genes)
  if (factor_overlap) {
    n_overlap_de_genes_per_group <- round(n_de_genes_per_group*.4)
    n_de_genes_per_group <- round(n_de_genes_per_group*.8)
  }

  n_groups <- (n_ctypes*n_processes) + 1

  # randomly assign cells by name to groups
  g_names <- sapply(1:n_groups,function(x){
    paste0('Group',x)
  })
  group_fracs <- rep((donors_per_process/n_ctypes)/donors_total,n_groups-1)
  group_fracs <- c(group_fracs,1-sum(group_fracs))
  meta$Group <- sample(g_names,nrow(meta),prob = group_fracs,replace = TRUE)

  # add data for each group
  de_groups <- data.frame(matrix(ncol=0,nrow=0))
  params <- setParam(params, "batchCells", n_cells*1.25)
  for (g in 1:(n_groups-1)) {
    # get the fractional counts for de genes specific to the current group
    de_df <- sim_helper(g,group_fracs[g],n_de_genes_per_group,params,meta,de_strength)

    # store the de data
    de_groups <- rbind(de_groups,de_df)
  }
  for (ct in 1:(n_ctypes)) {
    ct_groups <- seq(ct,n_groups-1,n_ctypes)

    # get the fractional counts for de genes found in same cell type but multiple processes
    de_df <- sim_helper(ct_groups,n_processes*group_fracs[1],
                        n_overlap_de_genes_per_group,params,meta,de_strength)

    # store the de data
    de_groups <- rbind(de_groups,de_df)
  }

  # replace part of old sim data with de data
  n_replace <- nrow(de_groups)
  sim_counts <- rbind(sim_counts[1:(nrow(sim_counts)-n_replace),],de_groups)
  sim_counts <- log1p(sim_counts*10000)

  # assign cells to cell types
  meta$ctypes <- NA
  ctype_names <- sapply(1:n_ctypes,function(x){paste0('ct',as.character(x))})
  group_ctypes <- rep(ctype_names,n_processes)
  for (i in 1:(n_groups-1)) {
    group_name <- paste0('Group',i)
    group_mask <- meta$Group==group_name
    meta$ctypes[group_mask] <- group_ctypes[i]
  }
  # randomly assign last group cells to cell types
  ncells_base <- sum(is.na(meta$ctypes))
  assigns <- sample(ctype_names,ncells_base,replace=T)
  meta$ctypes[is.na(meta$ctypes)] <- assigns

  # assign donors to cells
  meta$donors <- NA
  group_donors <- list()
  track <- 1
  for (i in 1:n_processes) {
    d_range <- c(track,track + donors_per_process - 1)
    for (j in 1:n_ctypes) {
      group_donors[[length(group_donors)+1]] <- d_range
    }
    track <- d_range[2] + 1
  }
  group_donors[[length(group_donors)+1]] <- c(track,donors_total)

  for (i in 1:(n_groups)) {
    group_name <- paste0('Group',i)
    group_mask = meta$Group==group_name
    ncells_group = sum(group_mask)
    d_range <- sapply(group_donors[[i]][1]:group_donors[[i]][2], function(x){paste0('s',as.character(x))})
    assigns <- sample(d_range,ncells_group,replace=T)
    meta$donors[group_mask] <- assigns
  }

  sim_meta <- meta[,c('Group','donors','ctypes')]

  return(list(sim_counts,sim_meta))
}

#' Get data for DE genes that will correspond with the main simulation
#'
#' @param groups numeric Splatter group numbers corresponding with the groups to give the DE genes to
#' @param g_frac numeric Decimal indicating the fractional size of the DE groups
#' @param n_de_genes_per_group numeric The desired number of DE genes to get
#' @param params SplatParams Parameters used for the simulation
#' @param meta DataFrame The colData of the simulation
#' @param de_strength numeric An average foldchange for the DE genes
#'
#' @return normalized counts matrix of DE genes for the specified groups
sim_helper <- function(groups,g_frac,n_de_genes_per_group,params,meta,de_strength) {
  scsim_de <- splatSimulateGroups(params,
                                  group.prob = c(g_frac,1-g_frac),
                                  de.prob = c(.1,0),
                                  de.downProb = c(0,0),
                                  de.facLoc = de_strength,
                                  verbose = FALSE)


  de_counts <- counts(scsim_de)

  # divide by total counts
  cell_sizes <- colSums(de_counts)
  cell_sizes <- sapply(cell_sizes,function(x) {
    if (x==0) {
      return(1)
    } else {
      return(x)
    }
  })
  de_counts <- sweep(de_counts,MARGIN=2,cell_sizes,'/')

  # limit to desired number of DE genes
  de_genes <- as.data.frame(rowData(scsim_de)@listData)
  de_genes <- de_genes[de_genes$DEFacGroup1!=1,]
  de_genes <- as.character(de_genes$Gene)
  if (length(de_genes) > n_de_genes_per_group) {
    de_genes <- de_genes[1:n_de_genes_per_group]
  }

  de_counts <- de_counts[de_genes,]

  # get group g cells
  de_meta <- colData(scsim_de)
  de_cells_only <- as.character(de_meta$Cell)[de_meta$Group=='Group1']

  # get cell names to assign group g cells to
  cur_group <- sapply(groups,function(x){paste0('Group',x)})
  group_cell_names <- as.character(meta$Cell)[meta$Group%in%cur_group]
  n_cells_pick <- length(group_cell_names)

  # sample de cells to the number of cells equal to that which are in the current group
  picked_de_cells <- sample(de_cells_only,n_cells_pick)

  # select just the picked cells
  de_counts_de_cells <- de_counts[,picked_de_cells]

  # rename picked de cells
  colnames(de_counts_de_cells) <- group_cell_names

  # select rest of non de cells and rename
  not_de_cells <- as.character(de_meta$Cell)[de_meta$Group=='Group2']
  group_cell_names <- as.character(meta$Cell)[!(meta$Group%in%cur_group)]
  n_cells_pick <- length(group_cell_names)
  picked_not_de_cells <- sample(not_de_cells,n_cells_pick)
  de_counts_not_de_cells <- de_counts[,picked_not_de_cells]
  colnames(de_counts_not_de_cells) <- group_cell_names

  # combine de and not de cells matrixes
  de_df <- cbind(de_counts_de_cells,de_counts_not_de_cells)

  # reorder columns to match that of orig data
  de_df <- de_df[,as.character(meta$Cell)]

  # rename genes so that they are unique
  rownames(de_df) <- sapply(rownames(de_df),function(x){
    group_names_concat <- paste(cur_group, collapse = '_')
    paste0(group_names_concat,"_",x)
  })

  return(de_df)
}







