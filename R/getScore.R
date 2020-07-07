# ------------------------------------------------------------------------
# 5. A function to compute cluster comparaisons between all clusterings (by K)
#-------------------------------------------------------------------------

# getScore() computes cluster comparaisons between all clusterings (by K) using MARI

getScore <- function(clsByKList, clCompScore, kVec, nsim){
  listOfComp <- combn(nsim,2)
  scoresRes <-do.call(cbind, lapply(clsByKList, function(k){
    score <- apply(listOfComp, MARGIN = 2, function(c) {
      intSec <- intersect(which(!is.na(k[,c[1]])), which(!is.na(k[,c[2]])))
      clCompScore(k[intSec,c[1]], k[intSec,c[2]])
    } )
    c(mean(score), sd(score), min(score), max(score))
  }))
  rownames(scoresRes) <- c("mean", "sd", "min", "max")
  colnames(scoresRes) <- kVec
  scoresRes
}
