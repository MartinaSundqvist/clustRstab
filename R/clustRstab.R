#' A function to compute cluster stability
#'
#' clustRstab() computes cluster stability as the arithmetic mean of the cluster comparision
#' scores. It also computes sd, min, et max du score (by K). It makes a simple plot
#'
#' @param data a data frame dsfdsf dqsdqs
#' @param kVec a vector for the successive number of group considered
#' @param perturbedDataFun a function for perturbing the original data
#' @param clAlgo a function: the clustering method used
#' @param clCompScore dqsds
#' @param nsim a scalar, the number of simulation performed to estimate the clustergin stability
#' @param baseLineCorrection dada
#' @param plot dazadz
#' @param nProp azdzd
#' @param pProp fsdfsd
#'
#' @import aricode
#' @import mclust
#' @import RPEnsemble
#' @import ggplot2
#' @import parallel
#' @importFrom stats cutree dist hclust kmeans rnorm sd
#'
#' @export
clustRstab <-  function(
  data,
  kVec,
  typeOfComp = "all",
  perturbedDataFun = clustRstab::subSample,
  clAlgo = clustRstab::clAlgoKmeans,
  clCompScore = aricode::MARI,
  nsim = 100,
  baseLineCorrection = FALSE,
  plot = TRUE,
  nProp = 0.7,
  pProp = 0.9,
  noiseGaussianMean = 0,
  noiseGaussianSD = 1,
  randProjMethod = "Haar",
  randProjDim = 20,
  mc.cores = 2
  ){


  if (isFALSE(typeOfComp %in% c("toInitial", "all", "random")))
    stop("invalid type of comparison:", paste("", typeOfComp), ". Valid comparisons:", paste("", c("toInitial", "all", "random") ))

  if (isTRUE(all.equal(aricode::ARI, clCompScore)) & isTRUE(baseLineCorrection)){
    warning("ARI do not need to be corrected by chance")}

  if (isTRUE(all.equal(aricode::MARI, clCompScore)) & isTRUE(baseLineCorrection)){
    warning("MARI do not need to be corrected by chance")}

  if (isTRUE(all.equal(aricode::NID, clCompScore)) & isFALSE(baseLineCorrection)){
    warning("NID need to be corrected by chance")}


  perturbedDataList <- getNsimPerturbedDataSets(data = data,
                                                perturbedDataFun = perturbedDataFun,
                                                nProp = nProp,
                                                pProp = pProp,
                                                noiseGaussianMean = noiseGaussianMean,
                                                noiseGaussianSD = noiseGaussianSD,
                                                randProjMethod = randProjMethod,
                                                randProjDim = randProjDim,
                                                nsim = nsim)

  clsList  <- getCl(perturbedDataList = perturbedDataList,
                    data = data,
                    kVec = kVec,
                    clAlgo = clAlgo,
                    mc.cores = mc.cores)


  clsByK   <- getClsByK(clsList = clsList,
                        kVec = kVec)

  clustRstab <- getScore(typeOfComp = typeOfComp,
                       clCompScore = clCompScore,
                       clsByKList = clsByK,
                       kVec = kVec,
                       nsim = nsim,
                       data = data,
                       clAlgo = clAlgo,
                       baseLineCorrection = baseLineCorrection,
                       mc.cores = mc.cores)


  if (plot == TRUE){
    clustRstabPlot <- clustRstabPlot(clustRstabObj = clustRstab, kVec = kVec, nsim = nsim)
    print(clustRstabPlot)
    #plot(y = clustRstab[1,], x = kVec, type = "l", ylim = c(0,1))
  }
  clustRstab
}
