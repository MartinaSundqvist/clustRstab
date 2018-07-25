


randProj_data <- function(data, param, nsim){
  # param = dimension
}

subVar_data <- function(data, param, nsim){
  # param = proportion
  p <- ncol(data)
  perturbed_data <- data[, sample.int(p, round(p*param))]
  perturbed_data
}

subVar_data <- function(data, param, nsim){
  # param = proportion
  p <- ncol(data)
  perturbed_data <- data[, sample.int(p, round(p*param))]
  perturbed_data
}


# subsampled_obs <- function(data, proportion){
# More complicated cuz needs to be taken into account when doing the classification comparaisons.
# Thus, a vector of the observations IDs needs to be given by this function
# so that the intersection of observations can be extracted for the clust comp.
#}
# mclapply(1:nsim, getOneClassif, props = props, clustering.method = clustering.method, mc.cores=mc.cores)

noised_data <- function(data, param, nsim) {
  # param = noise
}