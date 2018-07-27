#' Get score
#'
#' @import aricode parallel
#' @export


get.score <- function(clust.comp.score, Classifs, nsim, mc.cores) {
  pairs <- combn(x = nsim, m = 2, simplify=FALSE)
  do.call(rbind, mclapply(pairs, function(pair) {
    classif1 <- Classifs[[pair[1]]]
    classif2 <- Classifs[[pair[2]]]
    if (length(Classifs[[1]]) == 1) {classif1 <- matrix(classif1); classif2 <- matrix(classif2)}
    score.pairs = sapply(1:ncol(Classifs[[1]]), function(i) clust.comp.score(classif1[,i], classif2[,i]))
    score.pairs[is.nan(score.pairs)] <- 0 # Since group 1 rende des NaNs, but I'm not sure that this should be put here !!!! cuz differ from NID, ARI and maybe others?
    return(score.pairs)
    }, mc.cores = mc.cores))
}

#' Get score
#' @import aricode parallel
#' @export

get.score.boot <- function(clust.comp.score, Classifs, nsim, mc.cores){
  pairs <- combn(x = nsim, m = 2, simplify=FALSE)
  score.pairs.boot <- do.call(rbind, mclapply(pairs, function(pair) {
    classif1 <- Classifs[[pair[1]]]
    classif2 <- Classifs[[pair[2]]]
    if (length(Classifs[[1]]) == 1) {classif1 <- matrix(classif1); classif2 <- matrix(classif2)}
    score.pairs.boot <- sapply(1:ncol(Classifs[[1]]), function(i) mean(
      replicate(10, clust.comp.score(sample(classif1[,i]), classif2[,i]))))
    score.pairs.boot[is.nan(score.pairs.boot)] <- 0
    return(score.pairs.boot)
    }, mc.cores=mc.cores))
}






#
# score.pairs <- do.call(rbind, mclapply(pairs, function(pair) {
#   classif1 <- Classifs[[pair[1]]]
#   classif2 <- Classifs[[pair[2]]]
#   if (length(nb.grs) == 1) {classif1 <- matrix(classif1); classif2 <- matrix(classif2)}
#   return(sapply(1:length(nb.grs), function(i) f.score(classif1[,i], classif2[,i]))) #compute NID score
# }, mc.cores=mc.cores))
#
#
# get.score <- function(pair) {
#   classif1 <- Classifs[[pair[1]]]
#   classif2 <- Classifs[[pair[2]]]
#   if (length(nb.grs) == 1) {classif1 <- matrix(classif1); classif2 <- matrix(classif2)}
#   return(sapply(1:length(nb.grs), function(i) f.score(classif1[,i], classif2[,i]))) #compute NID score
# }
#
#
# score.pairs <- function(cl1, cl2, mc.cores) {
#   do.call(rbind, mclapply(combn(x = 5, m = 2, simplify=FALSE), get.score, pair, mc.cores=3))
# }

#
# score.pairs <- do.call(rbind, mclapply(pairs, function(pair) {
#   classif1 <- Classifs[[pair[1]]]
#   classif2 <- Classifs[[pair[2]]]
#   if (length(nb.grs) == 1) {classif1 <- matrix(classif1); classif2 <- matrix(classif2)}
#   return(sapply(1:length(nb.grs), function(i) f.score(classif1[,i], classif2[,i]))) #compute NID score
# }, mc.cores=mc.cores))
#
#
# for (pair in pairs){
#   classif1 <- Classifs[[pair[1]]]
#   classif2 <- Classifs[[pair[2]]]
#   if (length(nb.grs) == 1) {classif1 <- matrix(classif1); classif2 <- matrix(classif2)}
#   print(sapply(1:length(nb.grs), function(i) NID(classif1[,i], classif2[,i])))
#
# }







