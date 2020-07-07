#' A function to compute cluster stability
#'
#' clustRstab() computes cluster stability as the arithmetic mean of the cluster comparision
#' scores. It also computes sd, min, et max du score (by K). It makes a simple plot
#'
#' @param data a data frame dsfdsf dqsdqs
#' @param kVec a vector for the successive number of group considered
#' @param perturbedDataFun a function for perturbing the original data
#' @param clAlgo a function: the clustering method used
#' @param clCompScore dqsdsq
#' @param param a list, dqsdqs
#' @param nsim a scalar, the number of simulation performed to estimate the clustergin stability
#'
#' @import aricode
#' @export
clustRstab <-  function(
  data,
  kVec,
  perturbedDataFun = clustRstab::subsample,
  clAlgo = clustRstab::clAlgoKmeans,
  clCompScore = aricode::MARI,
  nsim = 100,
  param
  ){

  perturbedDataList <- lapply(1:nsim, function(i) getPerturbedData(data, perturbedDataFun, param, i))

  clsList  <- getCl(perturbedDataList, data = data, kVec = kVec, clAlgo = clAlgo)

  clsByK   <- getClsByK(clsList = clsList, kVec = kVec)

  MARIStab <- getScore(clsByKList = clsByK, clCompScore = clCompScore, kVec = kVec, nsim = nsim)[1,]

  plot(y = MARIStab, x = kVec, type = "l", ylim = c(0,1))
}

