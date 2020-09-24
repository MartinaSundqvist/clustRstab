# ------------------------------------------------------------------------
# 5. A function to compute cluster comparaisons between all clusterings (by K)
#-------------------------------------------------------------------------

# getScore() computes cluster comparaisons between all clusterings (by K) using MARI

#' @importFrom utils combn
getScore <- function(typeOfComp,
                     clCompScore,
                     clsByKList,
                     kVec,
                     nsim,
                     data,
                     clAlgo,
                     baseLineCorrection,
                     mc.cores){

  TYPEOFCOMP <- c("toInitial", "all", "random")
  iMeth <- pmatch(typeOfComp, TYPEOFCOMP)
  # if (is.na(iMeth))
  #   stop("invalid type of comparison:", paste("", typeOfComp), ". Valid comparisons:", paste("", TYPEOFCOMP ))
  # if (iMeth == -1)
  #   stop("ambiguous type of comparison", paste("", typeOfComp))

   if (iMeth == 1){
    scoresRes <- do.call(cbind, lapply(clsByKList, function(k){
      clInitial <- clAlgo(data, max(k, na.rm = T))
      score <- sapply(1:nsim, function(c) {
        intSec <- which(!is.na(k[,c]))
        if (baseLineCorrection == TRUE){

         if (isTRUE(all.equal(aricode::NID, clCompScore))) {
            # clCompScore(clInitial[intSec], k[intSec, c]) -
            #  mean(replicate(10, clCompScore(clInitial[intSec], sample(k[intSec, c]))))  +
            #  1
           clCompScore(clInitial[intSec], k[intSec, c]) /
             mean(replicate(10, clCompScore(clInitial[intSec], sample(k[intSec, c]))))

          }

        else {
              clCompScore(clInitial[intSec], k[intSec, c]) -
            mean(replicate(10, clCompScore(clInitial[intSec], sample(k[intSec, c]))))
          }
        }

        else if (baseLineCorrection == FALSE){
          clCompScore(clInitial[intSec], k[intSec, c])
        }
      } )
      c(mean(score), sd(score), min(score), max(score))
    }))
  }


 if (iMeth == 2){
  listOfComp <- combn(nsim,2)
  scoresRes <-do.call(cbind, lapply(clsByKList, function(k){
    score <- apply(listOfComp, MARGIN = 2, function(c) {
      intSec <- intersect(which(!is.na(k[,c[1]])), which(!is.na(k[,c[2]])))
      if (baseLineCorrection == TRUE){
        if (isTRUE(all.equal(aricode::NID, clCompScore))) {
          # clCompScore(k[intSec,c[1]], k[intSec,c[2]]) -
          #   mean(replicate(10, clCompScore(k[intSec,c[1]], sample(k[intSec,c[2]])))) +
          #   1
          clCompScore(k[intSec,c[1]], k[intSec,c[2]]) /
            mean(replicate(10, clCompScore(k[intSec,c[1]], sample(k[intSec,c[2]]))))
        }
        else {
          clCompScore(k[intSec,c[1]], k[intSec,c[2]]) -
            mean(replicate(10, clCompScore(k[intSec,c[1]], sample(k[intSec,c[2]]))))
        }
      }
      else if (baseLineCorrection == FALSE){
      clCompScore(k[intSec,c[1]], k[intSec,c[2]])
      }
    } )
    c(mean(score), sd(score), min(score), max(score))
  }))
 }


  if (iMeth == 3){
    listOfComp <- combn(nsim,2)[,sample(1:ncol(combn(nsim,2)), size = ncol(combn(nsim,2))/2)]
    scoresRes <-do.call(cbind, lapply(clsByKList, function(k){
      score <- apply(listOfComp, MARGIN = 2, function(c) {
        intSec <- intersect(which(!is.na(k[,c[1]])), which(!is.na(k[,c[2]])))
        if (baseLineCorrection == TRUE){
          if (isTRUE(all.equal(aricode::NID, clCompScore))) {
            # clCompScore(k[intSec,c[1]], k[intSec,c[2]]) -
            #   mean(replicate(10, clCompScore(k[intSec,c[1]], sample(k[intSec,c[2]])))) +
            #   1
            clCompScore(k[intSec,c[1]], k[intSec,c[2]]) /
              mean(replicate(10, clCompScore(k[intSec,c[1]], sample(k[intSec,c[2]]))))


          }
          else {
            clCompScore(k[intSec,c[1]], k[intSec,c[2]]) -
              mean(replicate(10, clCompScore(k[intSec,c[1]], sample(k[intSec,c[2]]))))
          }
        }
        else if (baseLineCorrection == FALSE){
          clCompScore(k[intSec,c[1]], k[intSec,c[2]])
        }
      } )
      c(mean(score), sd(score), min(score), max(score))
    }))
  }

  rownames(scoresRes) <- c("mean", "sd", "min", "max")
  colnames(scoresRes) <- kVec
  scoresRes
}


