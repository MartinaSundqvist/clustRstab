
#' call to hclust formatted for the clusterStability package
#'
#' @importFrom stats hclust dist
#' @export
clusterStab_hclust <- function(data, vec_ngroup, options_clustering) {

  options <- options_default_hclust
  options[names(options_clustering)] <- options_clustering

  hc_out <- hclust(dist(data, method = options$distance), method = options$aggregation)
  cl_out <- cutree(hc_out, vec_ngroup)
  cl_out
}

#' call to kmeans formatted for the clusterStability package
#'
#' @importFrom stats kmeans
#' @export

clusterStab_kmeans <- function(data, vec_ngroup, options_clustering) {

  options <- options_default_kmeans
  options[names(options_clustering)] <- options_clustering

  cl_out <- sapply(vec_ngroup, function(k) kmeans(data, k, iter.max = options$iter.max,
                                                  nstart = options$nstart,
                                                  algorithm = options$algorithm)$cluster)
  cl_out
}

#
# clusterStab_kmeans <- function(data, vec_ngroup, ...) {
#   options <- list(...)
#   cl_out <- sapply(vec_ngroup, function(k) kmeans(data, k)$cluster)
#   cl_out
# }


#' call to mclust formatted for the clusterStability package
#'
#' @import mclust
#' @export
clusterStab_gmm <- function(data, vec_ngroup, options_clustering) {

  options <- options_default_gmm
  options[names(options_clustering)] <- options_clustering

  # mclust()
  cl_out <- sapply(vec_ngroup, function(k) Mclust(data, G = k, modelNames = options$gmmModel)$classification)
  cl_out
}




