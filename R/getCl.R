# ------------------------------------------------------------------------
# 3. A function to compute clusterings for K = 2,...,Kmax for each perturbed dataset
#-------------------------------------------------------------------------

# getCl() compute a clustering from the perturbed dataset into 2:Kmax (=kVec) clusters.
# - This function take into account that not all observations are present, and organizes clustering in a n x length(kVec) df,
# - where NAs are added when the observation index has not been sampled

getCl <- function(perturbedDataList, data, kVec, clAlgo, mc.cores){
  n <- nrow(data)
  lapply(1:length(perturbedDataList), function(i) {
    cl <- sapply(kVec, function(k){
    clAlgo(data = perturbedDataList[[i]], k = k)})
    mat <- as.data.frame(matrix(ncol = length(kVec), nrow = n, NA)) # For when all samples have not been sampled!
    mat[rownames(cl), ] <- cl
    mat <- cbind(mat, rep(paste0("df.", i), n))
    colnames(mat) <- c(kVec, "df")
    mat
  })
}
