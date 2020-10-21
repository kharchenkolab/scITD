
#' Generate simulated data with subgroups of donors expressing multi-cell type processes
#' @import splatter
#' @importFrom SingleCellExperiment colData rowData counts
#'
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
#' @param rseed numeric The random seed to use (default=10)
#'
#' @return normalized gene by cell count matrix, raw gene by cell count matrix,
#' metadata data.frame, and matrix of differentially expressed genes
#' @export
#'
#' @examples
#' get_sim_data(donors_total=50, cells_per_donor=200, n_processes=2, donors_per_process=10,
#' n_ctypes=2, n_genes=2000, de_prob=0.05, de_strength=2, rseed=10)
get_sim_data <- function(donors_total, cells_per_donor, n_processes, donors_per_process,
                         n_ctypes=2, n_genes=2000, de_prob=0.05, de_strength=2, rseed=10){
  params <- newSplatParams()
  params <- setParam(params, "batchCells", donors_total * cells_per_donor)
  params <- setParam(params, "nGenes", n_genes)
  params <- setParam(params, "seed", rseed)

  n_groups <- n_ctypes*n_processes
  group_fracs <- rep((donors_per_process/n_ctypes)/donors_total,n_groups)
  group_fracs <- c(group_fracs,1-sum(group_fracs))

  scsim <- splatSimulateGroups(params, group.prob = group_fracs,
                               de.prob = c(rep(de_prob,n_groups),0),
                               de.downProb = c(rep(0,n_groups), 0),
                               de.facLoc = de_strength,
                               verbose = FALSE)

  # assign cells to cell types
  cell_meta <- colData(scsim)
  cell_meta$ctypes <- NA
  ctype_names <- sapply(1:n_ctypes,function(x){paste0('ct',as.character(x))})
  group_ctypes <- rep(ctype_names,n_processes)
  for (i in 1:n_groups) {
    group_name <- paste0('Group',i)
    group_mask = cell_meta$Group==group_name
    cell_meta$ctypes[group_mask] <- group_ctypes[i]
  }
  # randomly assign last group cells to cell types
  ncells_base <- sum(is.na(cell_meta$ctypes))
  assigns <- sample(ctype_names,ncells_base,replace=T)
  cell_meta$ctypes[is.na(cell_meta$ctypes)] <- assigns


  # assign donors to cells
  cell_meta$donors <- NA
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

  for (i in 1:(n_groups+1)) {
    group_name <- paste0('Group',i)
    group_mask = cell_meta$Group==group_name
    ncells_group = sum(group_mask)
    d_range <- sapply(group_donors[[i]][1]:group_donors[[i]][2], function(x){paste0('s',as.character(x))})
    assigns <- sample(d_range,ncells_group,replace=T)
    cell_meta$donors[group_mask] <- assigns
  }

  # extract the DE genes
  sim_var_act <- as.data.frame(rowData(scsim)@listData)
  sim_var_act <- sim_var_act[,c(1,5:ncol(sim_var_act))]
  rownames(sim_var_act) <- sim_var_act$Gene
  sim_var_act$Gene <- NULL
  sim_var_act <- sim_var_act != 1

  sim_counts <- counts(scsim)
  sim_norm_counts <- normalize_counts(sim_counts)
  sim_meta <- cell_meta[,c('Group','donors','ctypes')]

  return(list(sim_norm_counts,sim_counts,sim_meta,sim_var_act))
}


