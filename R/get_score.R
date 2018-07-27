#' Generate perturbed data sets by subsampling varibles
#'
#' @import aricode parallel
#' @export

score.pairs <- function(cl1, cl2, mc.cores) {
  do.call(rbind, mclapply(pairs, function(pair) {
  classif1 <- Classifs[[pair[1]]]
  classif2 <- Classifs[[pair[2]]]
  if (length(nb.grs) == 1) {classif1 <- matrix(classif1); classif2 <- matrix(classif2)}
  return(sapply(1:length(nb.grs), function(i) f.score(classif1[,i], classif2[,i]))) #compute NID score
}, mc.cores=3))
  }