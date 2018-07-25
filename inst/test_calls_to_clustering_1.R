rm(list = ls())
library(clusterStability)


classif <- clusterStab_kmeans(iris[,1:3], 1:10, iterations = 100,
                              initialisations = 1 ,
                              algorithm = "Hartigan-Wong")

classif <- clusterStab_gmm(iris[,1:3], 1:10, gmmModel = "EEI")

var.subsample <- function(data, prop){
  p <- ncols(data)
  data <- data[, sample.int(p, round(p*prop))]
}


my.function <- function(data, prop, data.generator, nsim, vec_ngroup, clustering.method){
  mclapply(1:nsim, function(i)
    data <- data.generator(data, prop = prop),
    clustering.method(data, vec_ngroup, distance = "euclidean", aggregation = "ward.D2" ),
    mc.cores = 10
  )
}



classif <- my.function( data = iris[,1:4], prop = 0.5, data.generator = var.subsample, nsim = 20, vec_ngroup = 1:10, clustering.method = clusterStab_hclust)


pairs <- combn(x = nsim, m = 2, simplify=FALSE) # pair des jeux de données, pour m = 2 : (1,2), (1,3), (1,4),...,(n-1, n), returned in a list
nb.pairs <-  length(pairs) # x*(x-1)/m


