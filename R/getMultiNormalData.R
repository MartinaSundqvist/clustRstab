# ------------------------------------------------------------------------
# 1. A function to create multinormal data
#-------------------------------------------------------------------------
getMultiNormalData <- function(nb.grs, obsPerGr, vec.mu, vec.sd,
                               n.col.signal, n.col.noise = 0){

  g <- list()
  for (gr in 1:nb.grs){
    g[[gr]] <- replicate(n = n.col.signal, rnorm(n=obsPerGr, mean = vec.mu[gr], sd = vec.sd[gr]))
  }

  mnData <- do.call(rbind, g)
  mnData <- cbind(mnData, replicate(n.col.noise, rnorm(nb.grs*obsPerGr, 0,1)))
  return(as.data.frame(mnData))
}
