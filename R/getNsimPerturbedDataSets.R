
# ------------------------------------------------------------------------
# 2. A function to generate perturbed versions of the dataset
#-------------------------------------------------------------------------

# getPerturbedData() generates perturbed datasets by subsampling varaibles and/or observations
# - The observations / variables indexes are sorted in increasing order
getNsimPerturbedDataSets <- function(data, perturbedDataFun, param, nsim){
  perturbedDataList <- lapply(1:nsim, function(i) perturbedDataFun(data, param))
  perturbedDataList
}

