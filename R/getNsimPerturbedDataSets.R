
# ------------------------------------------------------------------------
# 2. A function to generate perturbed versions of the dataset
#-------------------------------------------------------------------------

# getPerturbedData() generates perturbed datasets by subsampling varaibles and/or observations
# - The observations / variables indexes are sorted in increasing order
getNsimPerturbedDataSets <- function(data, perturbedDataFun, nsim,
                                     nProp, pProp,
                                     noiseGaussianMean, noiseGaussianSD,
                                     randProjMethod, randProjDim,
                                     ...){
  perturbedDataList <- lapply(1:nsim, function(i) perturbedDataFun(data,
                                                                   nProp = nProp, pProp = pProp,
                                                                   noiseGaussianMean = noiseGaussianMean,
                                                                   noiseGaussianSD = noiseGaussianMean,
                                                                   randProjMethod = randProjMethod,
                                                                   randProjDim = randProjDim, ...))
  perturbedDataList
}

