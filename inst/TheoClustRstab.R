#' A function to compute cluster stability
#'
#' clustRstab() computes cluster stability as the arithmetic mean of the cluster comparision
#' scores. It also computes sd, min, et max du score (by K). It makes a simple plot
#'
#
#'
#'
#'


source("R/clAlgoFuns.R")
source("R/clustRstabPlot.R")
# source("R/getCl.R")
source("R/getClsByK.R")
source("R/getScore.R")

# GET CL
getCl2 <- function(perturbedDataList, data, kVec, clAlgo, mc.cores){
  n <- nrow(data)
  lapply(1:length(perturbedDataList), function(i) {
    cl <- sapply(kVec, function(k){
      clAlgo(data = perturbedDataList[[i]], k = k)})
    rownames(cl) <- 1:n
    mat <- as.data.frame(matrix(ncol = length(kVec), nrow = n, NA)) # For when all samples have not been sampled!
    mat[rownames(cl), ] <- cl
    mat <- cbind(mat, rep(paste0("df.", i), n))
    colnames(mat) <- c(kVec, "df")
    mat
  })
}


clAlgoKmeans2 <- function(data, k){
  cl <- kmeans(x = data, centers = k)$cluster
  cl
}

# GET THEO CLUSTER
TheoClustRstab <-  function(
  dataList,
  data = dataList[[1]],
  kVec,
  typeOfComp = "all",
  clAlgo = clAlgoKmeans2,
  clCompScore = aricode::NID,
  nsim = 100,
  baseLineCorrection = FALSE,
  plot = TRUE,
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


  # perturbedDataList <- getNsimPerturbedDataSets(data = data,
  #                                               perturbedDataFun = perturbedDataFun,
  #                                               nProp = nProp,
  #                                               pProp = pProp,
  #                                               noiseGaussianMean = noiseGaussianMean,
  #                                               noiseGaussianSD = noiseGaussianSD,
  #                                               randProjMethod = randProjMethod,
  #                                               randProjDim = randProjDim,
  #                                               nsim = nsim)

  clsList  <- getCl2(perturbedDataList = dataList,
                    data = dataList[[1]],
                    kVec = kVec,
                    clAlgo = clustRstab::clAlgoKmeans,
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
