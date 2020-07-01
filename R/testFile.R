
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

data <- getMultiNormalData(4,10,c(2,4,8,-4), rep(1,4), 20)
n <- nrow(data)

getPerturbedData <- function(data, nProp, pProp, i){
  n <- dim(data)[1]
  p <- dim(data)[2]
  nIndex <- sort(sample(n, nProp*n))
  pIndex <- sort(sample(p, pProp*p))
  data[nIndex, pIndex]
  }

Kmax <- 6
kVec <- 2:Kmax
nsim=2
nProp <- 0.8
pProp <- 0.7

# getCl <- function(data, kVec, nProp, pProp,nsim){
#  lapply(1:nsim, function(i) sapply(kVec, function(k) kmeans(getPerturbedData(data, nProp, pProp,i), centers = k)$cluster))
#     }


cls <- getCl(data, kVec, nProp, pProp, nsim)

getCl <- function(data, kVec, nProp, pProp,nsim){
  lapply(1:nsim, function(i) {
    cl <- sapply(kVec, function(k)
      kmeans(getPerturbedData(data, nProp, pProp,i),
             centers = k)$cluster
     )
    mat <- as.data.frame(matrix(ncol = length(kVec), nrow = nrow(data), NA)) # For when all samples have not been sampled!
    mat[rownames(cl), ] <- cl
    colnames(mat) <- kVec
    mat
  })
}



clsDf <- do.call(what = cbind, args = cls)
clsDf2 <- clsDf[,colnames(clsDf) == 2 ]
ListOfComp <- combn(nsim,2)





# clsDf[,colnames(clsDf) == 2]
#
# t <- list()
# g=0
# for(k in kVec){
# print(k)
#     # g=g+1
#   #  t[g] <- clsDf[,colnames(clsDf) == k]
# }
#
#
#
#
# w <- do.call(what = cbind, args = cls)
#
# typeof(w)
# dim(w)
#
# do.call(what = cbind, args = cls)
#
#
# 1:n*(1:n %in% as.numeric(rownames(cls[[1]])))
#
# rownames(cls[[1]]) <- as.numeric(rownames(cls[[1]]))
#
# 1:n*(1:n %in% rownames(cls[[1]]))
# mat <- as.data.frame(matrix(ncol = length(kVec), nrow = n, NA))
# mat[rownames(cls[[1]]),] <- cls[[1]]
#
# mat2 <- as.data.frame(matrix(ncol = length(kVec), nrow = n, NA))
# mat2[rownames(cls[[2]]),] <- cls[[1]]
#
# NID(cls[[1]][intersect(names(cls[[1]][,1]), names(cls[[2]][,1])),1], cls[[2]][intersect(names(cls[[1]][,1]), names(cls[[2]][,1])),1])

