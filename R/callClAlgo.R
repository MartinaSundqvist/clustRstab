# Any function that:
#  - takes as imput a n x p dataset
#  - gives as output a n x 1 classirifaction in $k$ classes (A named $n$ vector with cluster belongings)

#' @export
clAlgoHCWard <- function(data, k){
  d <- dist(x = data)
  clTree <- hclust(d = d, method = "ward.D2")
  cl <- cutree(tree = clTree, k = k)
  cl
}

#' @export
clAlgoKmeans <- function(data, k){
  cl <- kmeans(x = data, centers = k, nstart = 10)$cluster
  cl
}

# clAlgoMclust <- function(data, k, ...){
#   # optModelName <- Mclust(data)$modelName
#   # cl <- Mclust(data, G = k, modelNames = optModelName, ...)$classification
#   cl <- Mclust(data, G = k, modelNames = "EEI", ...)$classification
#   cl
# }
# La convergence de Mclust semble dependre le jeux de donnée sousjacent et ne converge pas à chaque fois
# Portant il faut fixer le model name une fois pour tous, car sinon incomparable !








# Old version
#' #' call to hclust formatted for the clusterStability package
#' #'
#' #' @importFrom stats hclust dist
#' #' @export
#' clusterStab_hclust <- function(data, k, options_clustering) {
#'
#'   options <- options_default_hclust
#'   options[names(options_clustering)] <- options_clustering
#'
#'   hc_out <- hclust(dist(data, method = options$distance), method = options$method)
#'   cl_out <- cutree(hc_out, vec_ngroup)
#'   cl_out
#' }
#'
#' #' call to kmeans formatted for the clusterStability package
#' #'
#' #' @importFrom stats kmeans
#' #' @export
#'
#' clusterStab_kmeans <- function(data, vec_ngroup, options_clustering = list()) {
#'
#'   options <- options_default_kmeans
#'   options[names(options_clustering)] <- options_clustering
#'
#'   cl_out <- sapply(vec_ngroup, function(k) kmeans(data, k, iter.max = options$iter.max,
#'                                                   nstart = options$nstart,
#'                                                   algorithm = options$algorithm)$cluster)
#'   cl_out
#' }
#'
#' #' call to mclust formatted for the clusterStability package
#' #'
#' #' @import mclust
#' #' @export
#' clusterStab_gmm <- function(data, vec_ngroup, options_clustering) {
#'
#'   options <- options_default_gmm
#'   options[names(options_clustering)] <- options_clustering
#'
#'   # mclust()
#'   cl_out <- sapply(vec_ngroup, function(k) Mclust(data, G = k, modelNames = options$modelNames)$classification)
#'   cl_out
#' }





