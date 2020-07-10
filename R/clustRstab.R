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
#' @importFrom stats cutree dist hclust kmeans rnorm sd
#'
#' @export
clustRstab <-  function(
  data,
  kVec,
  typeOfComp = "all",
  perturbedDataFun = clustRstab::subsample,
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
  randProjDim = 10
  ){

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
                    clAlgo = clAlgo)


  clsByK   <- getClsByK(clsList = clsList,
                        kVec = kVec)

  if (baseLineCorrection == FALSE){
    clustRstab <- getScore(typeOfComp = typeOfComp,
                       clCompScore = clCompScore,
                       clsByKList = clsByK,
                       kVec = kVec,
                       nsim = nsim,
                       data = data,
                       clAlgo = clAlgo)
  }

  if (baseLineCorrection == TRUE){
    clustRstabInit <- getScore(typeOfComp = typeOfComp,
                           clCompScore = clCompScore,
                           clsByKList = clsByK,
                           kVec = kVec,
                           nsim = nsim,
                           data = data,
                           clAlgo = clAlgo)

    clustRstabBaseLine <- getBaseLineScore(typeOfComp = typeOfComp,
                                     clCompScore = clCompScore,
                                     clsByKList = clsByK,
                                     kVec = kVec,
                                     nsim = nsim,
                                     data = data,
                                     clAlgo = clAlgo)

    clustRstab <- as.data.frame(data.matrix(clustRstabInit) - data.matrix(clustRstabBaseLine))
    clustRstab[1,] <- clustRstab[1,] + 1
    colnames(clustRstab) <- colnames(clustRstabInit)
    rownames(clustRstab) <- rownames(clustRstabInit)
  }
  if (plot == TRUE){
    clustRstabPlot <- clustRstabPlot(clustRstabObj = clustRstab, kVec = kVec, nsim = nsim)
    print(clustRstabPlot)
    #plot(y = clustRstab[1,], x = kVec, type = "l", ylim = c(0,1))
  }
  clustRstab
}
