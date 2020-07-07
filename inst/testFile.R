# This code tests the architecture of the clustRstab package
# It uses:
# (1) subsampling of obs. et vars. in order to pertube the data
# (2) k-means to compute the clustering
# (3) MARI to compute cluster comparisons
# (4) Cluster comparisions are computed between all obtained clusterings (for K = 2,..,Kmax)

# Libraries
library(aricode)

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



# ------------------------------------------------------------------------
# 2. A function to generate perturbed versions of the dataset
#-------------------------------------------------------------------------

# getPerturbedData() generates perturbed datasets by subsampling varaibles and/or observations
# - The observations / variables indexes are sorted in increasing order
getPerturbedData <- function(data, nProp, pProp, i){
  n <- dim(data)[1]
  p <- dim(data)[2]
  nIndex <- sort(sample(n, nProp*n))
  pIndex <- sort(sample(p, pProp*p))
  data[nIndex, pIndex]
  }


# ------------------------------------------------------------------------
# 3. A function to compute clusterings for K = 2,...,Kmax for each perturbed dataset
#-------------------------------------------------------------------------

# getCl() compute a clustering from the perturbed dataset into 2:Kmax (=kVec) clusters.
# - This function take into account that not all observations are present, and organizes clustering in a n x length(kVec) df,
# - where NAs are added when the observation index has not been sampled

getCls <- function(pertubedDataList, kVec, n){
  lapply(pertubedDataList, function(i) {
    cl <- sapply(kVec, function(k)
      kmeans(i, centers = k)$cluster
     )
    mat <- as.data.frame(matrix(ncol = length(kVec), nrow = n, NA)) # For when all samples have not been sampled!
    mat[rownames(cl), ] <- cl
    mat <- cbind(mat, rep(paste0("df.", i), n))
    colnames(mat) <- c(kVec, "df")
    mat
  })
}


# ------------------------------------------------------------------------
# 4. A function to execute the perturbation and clustering of the data, and organize it by K
#-------------------------------------------------------------------------

# getClsByK() perturbes the data and compute clustering, organize the results after corresponding K

# - Get a list where each item corresponds to the classification of perturbed dataset df.x. cols = kVec, rowns = 1:n
# - From this create a list of length(kVec) items, where each item is n*nsim x 2 df
# - These df.s are then split each n (corresponding to each perturbed df), and cocanated in a n x length(kVec) df
# - These tables are then concatenated by rbind to a n*nsim x length(kVec) + 1 data.frame

getClsByK <- function(clsList, kVec){
  lapply(lapply(kVec, function(k)
    do.call(what = rbind,
          args = clsList
          )[,c(as.character(k), "df")]), function(j)
            sapply(split(x = j, f = j$df), function(i)
              i[,-ncol(i)]
              )
    )
}

# ------------------------------------------------------------------------
# 5. A function to compute cluster comparaisons between all clusterings (by K)
#-------------------------------------------------------------------------

# getScore() computes cluster comparaisons between all clusterings (by K) using MARI

getScore <- function(clsByK, kVec, nsim){
  listOfComp <- combn(nsim,2)
  scoresRes <-do.call(cbind, lapply(clsByK, function(k){
    score <- apply(listOfComp, MARGIN = 2, function(c) {
    intSec <- intersect(which(!is.na(k[,c[1]])), which(!is.na(k[,c[2]])))
    MARI(k[intSec,c[1]], k[intSec,c[2]])
    } )
    c(mean(score), sd(score), min(score), max(score))
  }))
  rownames(scoresRes) <- c("mean", "sd", "min", "max")
  colnames(scoresRes) <- kVec
  scoresRes
}

# ------------------------------------------------------------------------
# 6. A function to compute cluster stability (by K)
#-------------------------------------------------------------------------

# clustRstab() computes cluster stability as the arithmic mean of the cluster comparaision scores (here MARI)
# It also computes sd, min, et max du score (by K)
# It makes a simple plot

clusteRstab <-  function(data, kVec, nProp, pProp, nsim){
  n                 <- nrow(data)
  pertubedDataList  <- lapply(1:nsim, function(i) getPerturbedData(dat, nProp, pProp,i))
  clsList           <- getCls(pertubedDataList = pertubedDataList, kVec = kVec, n = n)
  clsByK            <- getClsByK(clsList = clsList, kVec = kVec)
  MARIStab          <- getScore(clsByK = clsByK, kVec = kVec, nsim = nsim)
  plot(y = MARIStab[1,], x = kVec, type = "l")
  MARIStab
}


# ------------------------------------------------------------------------
# 7. Execution
#-------------------------------------------------------------------------

# 7a. Simulation Parameters For Intitial Data
trueK <- 5
nPerGr <- 20
grMeans <- c(2,4,8,-4,20)
sim.sd <- 1
grSD <- rep(sim.sd, trueK)
vars <- 10

# 7.b Cluster Stability Parammeters
Kmax1 <- 8
kVec2 <- 2:Kmax1
nsim2 = 50
nProp2 <- 0.9
pProp2 <- 0.5

# 7.c Generate the Inital DAta
dat <- getMultiNormalData(nb.grs = trueK,
                           obsPerGr = nPerGr,
                           vec.mu = grMeans,
                           vec.sd = grSD,
                           n.col.signal = vars)

# 7.d Compute Cluster Stability
start_time <- Sys.time()
print(clusteRstab(data= dat,
                  kVec = kVec2,
                  nProp = nProp2,
                  pProp =pProp2,
                  nsim = nsim2))
end_time <- Sys.time()
print(end_time - start_time)

# Mettre ça dans une function paralelle (?)
# system.time({ sleep_for_a_minute() })


# ------------------------------------------------------------------------
# OLD !!!! !!!! !!!! !!!!
# -------------------------------------------------------------------------
#
#
#
# # This code tests the architecture of the clustRstab package
# # It uses:
# # (1) subsampling of obs. et vars. in order to pertube the data
# # (2) k-means to compute the clustering
# # (3) MARI to compute cluster comparisons
# # (4) Cluster comparisions are computed between all obtained clusterings (for K = 2,..,Kmax)
#
# # Libraries
# library(aricode)
#
#
# # ------------------------------------------------------------------------
# # 1. A function to create multinormal data
# #-------------------------------------------------------------------------
# getMultiNormalData <- function(nb.grs, obsPerGr, vec.mu, vec.sd,
#                                n.col.signal, n.col.noise = 0){
#
#   g <- list()
#   for (gr in 1:nb.grs){
#     g[[gr]] <- replicate(n = n.col.signal, rnorm(n=obsPerGr, mean = vec.mu[gr], sd = vec.sd[gr]))
#   }
#
#   mnData <- do.call(rbind, g)
#   mnData <- cbind(mnData, replicate(n.col.noise, rnorm(nb.grs*obsPerGr, 0,1)))
#   return(as.data.frame(mnData))
# }
#
#
#
# # ------------------------------------------------------------------------
# # 2. A function to generate perturbed versions of the dataset
# #-------------------------------------------------------------------------
#
# # getPerturbedData() generates perturbed datasets by subsampling varaibles and/or observations
# # - The observations / variables indexes are sorted in increasing order
# getPerturbedData <- function(data, nProp, pProp, i){
#   n <- dim(data)[1]
#   p <- dim(data)[2]
#   nIndex <- sort(sample(n, nProp*n))
#   pIndex <- sort(sample(p, pProp*p))
#   data[nIndex, pIndex]
# }
#
#
# # ------------------------------------------------------------------------
# # 3. A function to compute clusterings for K = 2,...,Kmax for each perturbed dataset
# #-------------------------------------------------------------------------
#
# # getCl() compute a clustering from the perturbed dataset into 2:Kmax (=kVec) clusters.
# # - This function take into account that not all observations are present, and organizes clustering in a n x length(kVec) df,
# # - where NAs are added when the observation index has not been sampled
#
# getCl <- function(data, kVec, nProp, pProp,nsim){
#   n <- nrow(data)
#   lapply(1:nsim, function(i) {
#     cl <- sapply(kVec, function(k)
#       kmeans(getPerturbedData(data, nProp, pProp,i),
#              centers = k)$cluster
#     )
#     mat <- as.data.frame(matrix(ncol = length(kVec), nrow = nrow(data), NA)) # For when all samples have not been sampled!
#     mat[rownames(cl), ] <- cl
#     # mat <- cbind(mat, 1:n, rep(paste0("df.",i), n))
#     # colnames(mat) <- c(kVec, "id", "df")
#     mat <- cbind(mat, rep(paste0("df.", i), n))
#     colnames(mat) <- c(kVec, "df")
#     mat
#   })
# }
#
#
# # ------------------------------------------------------------------------
# # 4. A function to execute the perturbation and clustering of the data, and organize it by K
# #-------------------------------------------------------------------------
#
# # getClsByK() perturbes the data and compute clustering, organize the results after corresponding K
#
# # - Get a list where each item corresponds to the classification of perturbed dataset df.x. cols = kVec, rowns = 1:n
# # - From this create a list of length(kVec) items, where each item is n*nsim x 2 df
# # - These df.s are then split each n (corresponding to each perturbed df), and cocanated in a n x length(kVec) df
# # - These tables are then concatenated by rbind to a n*nsim x length(kVec) + 1 data.frame
#
# getClsByK <- function(data, kVec, nProp, pProp, nsim){
#   lapply(lapply(kVec, function(k)
#     do.call(what = rbind,
#             args = getCl(data = data,
#                          kVec = kVec,
#                          nProp = nProp,
#                          pProp = pProp,
#                          nsim = nsim)
#     )[,c(as.character(k), "df")]), function(j)
#       sapply(split(x = j, f = j$df), function(i)
#         i[,-ncol(i)]
#       )
#   )
# }
#
# # ------------------------------------------------------------------------
# # 5. A function to compute cluster comparaisons between all clusterings (by K)
# #-------------------------------------------------------------------------
#
# # getScore() computes cluster comparaisons between all clusterings (by K) using MARI
#
# getScore <- function(clsByK, kVec, nsim){
#   listOfComp <- combn(nsim,2)
#   scoresRes <-do.call(cbind, lapply(clsByK, function(k){
#     score <- apply(listOfComp, MARGIN = 2, function(c) {
#       intSec <- intersect(which(!is.na(k[,c[1]])), which(!is.na(k[,c[2]])))
#       MARI(k[intSec,c[1]], k[intSec,c[2]])
#     } )
#     c(mean(score), sd(score), min(score), max(score))
#   }))
#   rownames(scoresRes) <- c("mean", "sd", "min", "max")
#   colnames(scoresRes) <- kVec
#   scoresRes
# }
#
# # ------------------------------------------------------------------------
# # 6. A function to compute cluster stability (by K)
# #-------------------------------------------------------------------------
#
# # clustRstab() computes cluster stability as the arithmic mean of the cluster comparaision scores (here MARI)
# # It also computes sd, min, et max du score (by K)
# # It makes a simple plot
#
# clusteRstab <-  function(data, kVec, nProp, pProp, nsim){
#   clsByK <- getClsByK(data = data, kVec = kVec, nProp = nProp, pProp = pProp, nsim = nsim)
#   MARIStab <- getScore(clsByK = clsByK, kVec = kVec, nsim = nsim)
#   plot(y = MARIStab[1,], x = kVec, type = "l")
#   MARIStab
# }
#
#
# # ------------------------------------------------------------------------
# # 7. Execution
# #-------------------------------------------------------------------------
#
# # 7a. Simulation Parameters For Intitial Data
# trueK <- 5
# nPerGr <- 20
# grMeans <- c(2,4,8,-4,20)
# sim.sd <- 1
# grSD <- rep(sim.sd, trueK)
# vars <- 100
#
# # 7.b Cluster Stability Parammeters
# Kmax1 <- 8
# kVec1 <- 2:Kmax1
# nsim1 = 50
# nProp1 <- 0.9
# pProp1 <- 0.5
#
# # 7.c Generate the Inital DAta
# dat <- getMultiNormalData(nb.grs = trueK,
#                           obsPerGr = nPerGr,
#                           vec.mu = grMeans,
#                           vec.sd = grSD,
#                           n.col.signal = vars)
#
# # 7.d Compute Cluster Stability
# start_time <- Sys.time()
# print(clusteRstab(data= dat, kVec = kVec1, nProp = nProp1, pProp =pProp1, nsim = nsim1))
# end_time <- Sys.time()
# print(end_time - start_time)
#
# # Mettre ça dans une function paralelle (?)
# # system.time({ sleep_for_a_minute() })
#
#
#
#


