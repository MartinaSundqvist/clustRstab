# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Perform clustering stability
#'
#' @param data a data frame
#'
#' @import aricode mclust parallel
#' @export
#'
source("R/get_perturbed_data.R")
source("R/calls_to_clustering.R")
clust.stab.index <- function(data, nb.grs = 1:10, data.generator, data.genenerator.param, clustering.method, nsim =100, mc.cores = 3, score = "NID", ...){
  if (score == "NID") f.score = NID
  if (score == "ARI") f.score = ARI


  # 1. Import data
  X <- data

  # 2. Generate perturbed datasets
  Perturbed_data <- mclapply(1:nsim, data.generator, data = X, param = data.genenerator.param, mc.cores=mc.cores)

  # 3. Compute clusterings for each perturbed dataset for K = vec_ngroup
  Classifs <- mclapply(Perturbed_data,
                       clustering.method,
                       vec_ngroup = nb.grs,
                       mc.cores=mc.cores)


  # 4. Compute pairwise comparaisons

  #   4.a Define all the possible pairs of comparaison
  pairs <- combn(x = nsim, m = 2, simplify=FALSE) # pair des jeux de données, pour m = 2 : (1,2), (1,3), (1,4),...,(n-1, n), returned in a list
  nb.pairs <-  length(pairs) # x*(x-1)/m

  #   4.b Compute pairwise comparisons
  score.pairs <- do.call(rbind, mclapply(pairs, function(pair) {
      classif1 <- Classifs[[pair[1]]]
      classif2 <- Classifs[[pair[2]]]
    if (length(nb.grs) == 1) {classif1 <- matrix(classif1); classif2 <- matrix(classif2)}
    return(sapply(1:length(nb.grs), function(i) f.score(classif1[,i], classif2[,i]))) #compute NID score
  }, mc.cores=mc.cores)) # use a nb of mc.cores

  score.pairs[is.nan(score.pairs)] <- 0 # Since group 1 rende des NaNs

  # 5. Permute labels and compare clusters in order to have a measure of randomly shared information

  score.pairs.boot <- do.call(rbind, mclapply(pairs, function(pair) {
      classif1 <- Classifs[[pair[1]]]
      classif2 <- Classifs[[pair[2]]]
    if (length(nb.grs) == 1) {classif1 <- matrix(classif1); classif2 <- matrix(classif2)}
    return(sapply(1:length(nb.grs), function(i) mean(
      replicate(10, f.score(sample(classif1[,i]), classif2[,i])))
    ))
  }, mc.cores=mc.cores))
  score.pairs.boot[is.nan(score.pairs.boot)] <- 0

  # 6. Normalize cluster comparaisons by the permuted comparaisons

  # Instab index and sd
  mean <- colMeans(score.pairs)
  sd <- apply(score.pairs, 2, sd)

  # Normalized Instab index and sd
  if (score == "NID"){ mean.normalized <- colMeans(score.pairs - score.pairs.boot) + 1}
  else if (score == "ARI"){ mean.normalized <- colMeans(score.pairs - score.pairs.boot)}
  sd.normalized <- apply(t(score.pairs) - colMeans(score.pairs.boot), 1, sd)

  # mean.normalized <- colMeans(score.pairs - score.pairs.boot) + 1

  # 7. Return all needed info
  return(list(instab = mean,
              instab.sd = sd,
              instab.norm = mean.normalized,
              instab.norm.sd = sd.normalized,
              nsim = nsim,
              data.genenerator.param = data.genenerator.param,
              method = clustering.method,
              #instab.plot = instab.plot,
              nb.grs = nb.grs))

}



# TEST
# dat <- iris[,1:4]
# clust.stab.index(dat, data.generator = subVar_data, nsim = 10, data.genenerator.param =0.5, clustering.method = clusterStab_kmeans, score = "NID")









##############################################################################################################################
# Initial function

# instab.index2 <- function(data, nb.grs = 1:10, data.gen, clustering.method, nsim =100, mc.cores = 3, score = "NID", ...){
#
#
#   # 1. Import data
#   X <- data
#   p <- ncol(X) # variables
#   n <- nrow(X) # observation
#   if (score == "NID") f.score = NID
#   if (score == "ARI") f.score = ARI
#
#   getOneClassif <- function(i, G=nb.grs, props, clustering.method){
#     # 2. Produce new datasets
#     #     2.a Subsample obs/var
#     Xsub <- X[sample.int(n, round(n*props[1])), sample.int(p, round(p*props[2]))]
#     Xsub <- X[, sample.int(p, round(p*props[2]))]
#     #     2.b Introduce noise
#     #     2.c Project in lower space
#
#     # 3. Cluster each new dataset by using a specific clustering algo
#     #   if (method == 'GMM'){classif <- sapply(G, function(g) Mclust(Xsub, G = g, ...)$classification)} # modelNames="EII", If model name not choosen the running will be very long
#     #   else if (method == 'HAC') {return(hclust(dist(Xsub, ...), ...))} # method="ward.D2", otherwise par default "complete", the distance can also be choosen if eucledian is not wanted
#     #   else if (method == 'Kmeans') {classif <- sapply(G, function(g) kmeans(Xsub, centers = g, ...)$cluster)}
#     #   # else if (method == 'user_function') {classif <- }
#     # }
#     #
#
#     if (clustering.method == 'GMM'){classif <- sapply(G, function(g) Mclust(Xsub, G = g, modelNames="EII")$classification)} # , If model name not choosen the running will be very long
#     else if (clustering.method == 'HAC') {return(hclust(dist(Xsub), method="ward.D2"))} # , otherwise par default "complete", the distance can also be choosen if eucledian is not wanted
#     else if (clustering.method == 'Kmeans') {classif <- sapply(G, function(g) kmeans(Xsub, centers = g)$cluster)}
#     # else if (method == 'user_function') {classif <- }
#   }
#
#   perturbed_data <- mclapply(1:nsim, subVar_data, data = X)
#
#   # Get One Classif for each nb of groups and each simulated data set
#   Classifs <- mclapply(1:nsim, getOneClassif, props = props, clustering.method = clustering.method, mc.cores=mc.cores)
#
#   # 4. Compute pairwise comparaisons
#
#   #   4.a Define all the possible pairs of comparaison
#   pairs <- combn(x = nsim, m = 2, simplify=FALSE) # pair des jeux de données, pour m = 2 : (1,2), (1,3), (1,4),...,(n-1, n), returned in a list
#   nb.pairs <-  length(pairs) # x*(x-1)/m
#
#   #   4.b Compute pairwise comparisons
#   score.pairs <- do.call(rbind, mclapply(pairs, function(pair) {
#     if (clustering.method =='HAC') {
#       classif1 <- cutree(Classifs[[pair[1]]], nb.grs) #for each couple of trees (pair1 et pair 2), cut the tree at the levels 1 to n-1, je change pour 1:10
#       classif2 <- cutree(Classifs[[pair[2]]], nb.grs)
#     }
#     else {
#       classif1 <- Classifs[[pair[1]]]
#       classif2 <- Classifs[[pair[2]]]
#     }
#     if (length(nb.grs) == 1) {classif1 <- matrix(classif1); classif2 <- matrix(classif2)}
#     return(sapply(1:length(nb.grs), function(i) f.score(classif1[,i], classif2[,i]))) #compute NID score
#   }, mc.cores=mc.cores)) # use a nb of mc.cores
#
#   score.pairs[is.nan(score.pairs)] <- 0 # Since group 1 rende des NaNs
#
#   # 5. Permute labels and compare clusters in order to have a measure of randomly shared information
#
#   score.pairs.boot <- do.call(rbind, mclapply(pairs, function(pair) {
#     if (clustering.method =='HAC'){
#       classif1 <- cutree(Classifs[[pair[1]]], nb.grs) #for each couple of trees (pair1 et pair 2), cut the tree at the levels 1 to n-1, je change pour 1:10
#       classif2 <- cutree(Classifs[[pair[2]]], nb.grs)
#     }
#     else {
#       classif1 <- Classifs[[pair[1]]]
#       classif2 <- Classifs[[pair[2]]]
#     }
#     if (length(nb.grs) == 1) {classif1 <- matrix(classif1); classif2 <- matrix(classif2)}
#
#     return(sapply(1:length(nb.grs), function(i) mean(
#       replicate(10, f.score(sample(classif1[,i]), classif2[,i])))
#     ))
#   }, mc.cores=mc.cores))
#   score.pairs.boot[is.nan(score.pairs.boot)] <- 0
#
#   # 6. Normalize cluster comparaisons by the permuted comparaisons
#
#   # Instab index and sd
#   mean <- colMeans(score.pairs)
#   sd <- apply(score.pairs, 2, sd)
#
#   # Normalized Instab index and sd
#   if (score == "NID"){ mean.normalized <- colMeans(score.pairs - score.pairs.boot) + 1}
#   else if (score == "ARI"){ mean.normalized <- colMeans(score.pairs - score.pairs.boot)}
#   sd.normalized <- apply(t(score.pairs) - colMeans(score.pairs.boot), 1, sd)
#
#   # mean.normalized <- colMeans(score.pairs - score.pairs.boot) + 1
#
#   # 7. Return all needed info
#   return(list(instab = mean,
#               instab.sd = sd,
#               instab.norm = mean.normalized,
#               instab.norm.sd = sd.normalized,
#               nsim = nsim,
#               props = props,
#               method = clustering.method,
#               #instab.plot = instab.plot,
#               nb.grs = nb.grs))
# }

