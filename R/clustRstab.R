# ------------------------------------------------------------------------
# 6. A function to compute cluster stability (by K)
#-------------------------------------------------------------------------

# clustRstab() computes cluster stability as the arithmic mean of the cluster comparaision scores (here MARI)
# It also computes sd, min, et max du score (by K)
# It makes a simple plot

clusteRstab <-  function(data, kVec, perturbedDataFun, clAlgo, clCompScore = MARI, param, nsim){
  perturbedDataList <- lapply(1:nsim, function(i) getPerturbedData(data, perturbedDataFun, param, i))
  clsList <- getCl(perturbedDataList, data = data, kVec = kVec, clAlgo = clAlgo)
  clsByK <- getClsByK(clsList = clsList, kVec = kVec)
  MARIStab <- getScore(clsByKList = clsByK, clCompScore = clCompScore, kVec = kVec, nsim = nsim)[1,]
  plot(y = MARIStab, x = kVec, type = "l", ylim = c(0,1))
}
