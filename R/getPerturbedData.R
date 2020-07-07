
# ------------------------------------------------------------------------
# 2. A function to generate perturbed versions of the dataset
#-------------------------------------------------------------------------

# getPerturbedData() generates perturbed datasets by subsampling varaibles and/or observations
# - The observations / variables indexes are sorted in increasing order
getPerturbedData <- function(data, perturbedDataFun, param, i, ...){
  perturbedData <- perturbedDataFun(data, param, ...)
  perturbedData
}
