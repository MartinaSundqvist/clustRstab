#' Generate perturbed data sets by subsampling varibles
#'
#' @export
subVar_data <- function(i, data, param){
  # param = proportion
  p <- ncol(data)
  perturbed_data <- data[, sample.int(p, round(p*param))]
  perturbed_data
}

randProj_data <- function(i, data, param){
  # param = dimension
}

noised_data <- function(i, data, param) {
  # param = noise
}


# subsampled_obs <- function(data, proportion){
# More complicated cuz needs to be taken into account when doing the classification comparaisons.
# Thus, a vector of the observations IDs needs to be given by this function
# so that the intersection of observations can be extracted for the clust comp.
#}
