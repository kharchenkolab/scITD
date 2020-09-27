
#' Form the tensor and scale the data
#'
#' @param container environment Project container that stores sub-containers
#' for each cell type as well as results and plots from all analyses
#' @param var_scale_power numeric Exponent of normalized variance that is
#' used for variance scaling. Variance for each gene
#' is initially set to unit variance across donors (for a given cell type).
#' Variance for each gene is then scaled by multiplying the unit scaled values
#' by each gene's normalized variance (where the effect of the mean-variance
#' dependence is taken into account) to the exponent specified here.
#' If NULL, uses var_scale_power from container$experiment_params. (default=NULL)
#'
#' @return the project container with tensor data added in the
#' container$tensor_data slot
#' @export
form_tensor <- function(container, var_scale_power=NULL) {
  ctypes_use <- container$experiment_params$ctypes_use
  scale_var <- container$experiment_params$scale_var

  if (scale_var) {
    if (is.null(var_scale_power)) {
      var_scale_power <- container$experiment_params$var_scale_power
      if (is.null(var_scale_power)) {
        print('need to supply var_scale_power - cannot run')
        # throw error here
      }
    }
  }

  # first get donors present in all ctypes
  donors_in_all <- container$scMinimal_ctype[[ctypes_use[1]]]$donors
  for (i in 2:length(ctypes_use)) {
    ct <- ctypes_use[i]
    donors_in_all <- intersect(donors_in_all,container$scMinimal_ctype[[ct]]$donors)
  }

  tnsr <- array(NA, dim = c(length(donors_in_all),
                            ncol(container$scMinimal_ctype[[1]]$data_means),
                            length(ctypes_use)))

  for (i in 1:length(ctypes_use)) {
    ct <- ctypes_use[i]
    donor_means <- container$scMinimal_ctype[[ct]]$data_means
    donor_means <- donor_means[donors_in_all,]

    # center with unit variance
    donor_means <- scale(donor_means, center=TRUE)

    # to change gene variance across donors by normalized variability
    if (scale_var) {
      scale_factor <- get_normalized_variance(container$scMinimal_ctype[[ct]])
      donor_means <- apply(donor_means,MARGIN=1,function(x) {
        x * (scale_factor ** var_scale_power)
      })
      donor_means <- t(donor_means)
    }

    tnsr[, ,i] <- as.matrix(donor_means)
  }

  container$tensor_data <- list(donors_in_all, colnames(donor_means), ctypes_use, tnsr)
  return(container)
}














