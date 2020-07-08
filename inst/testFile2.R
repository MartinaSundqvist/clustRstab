# This code tests the architecture of the clustRstab package
# It uses:
# (1) subsampling of obs. et vars. in order to pertube the data
# (2) k-means to compute the clustering
# (3) MARI to compute cluster comparisons
# (4) Cluster comparisions are computed between all obtained clusterings (for K = 2,..,Kmax)

# Libraries
# library(aricode)

# ------------------------------------------------------------------------
# 7. Execution
#-------------------------------------------------------------------------

# 7a. Simulation Parameters For Intitial Data
trueK <- 5
nPerGr <- 10
grMeans <- c(2, 6, 10, 20, 15)
sim.sd <- 1
grSD <- rep(sim.sd, trueK)
vars <- 100



# 7.b Cluster Stability Parammeters
Kmax1 <- 8
kVec1 <- 2:Kmax1
nsim1 <- 20
nProp1 <- 0.7
pProp1 <- 0.5
clAlgo1 <- clAlgoGmmEEI
clCompScore1 = aricode::MARI
perturbedDataFun1 = noiseGaussian

getMultiNormalData <- function(nb.grs, obsPerGr, vec.mu, vec.sd,
                               n.col.signal, n.col.noise = 0){

  g <- list()
  for (gr in 1:nb.grs){
    g[[gr]] <- replicate(n = n.col.signal, rnorm(n=obsPerGr, mean = vec.mu[gr], sd = vec.sd[gr]))
  }

  mnData <- do.call(rbind, g)
  mnData <- cbind(mnData, replicate(n.col.noise, rnorm(nb.grs*obsPerGr, 0,1)))

  as.data.frame(mnData)
}

# 7.c Generate the Inital Data
dat <- getMultiNormalData(nb.grs = trueK,
                           obsPerGr = nPerGr,
                           vec.mu = grMeans,
                           vec.sd = grSD,
                           n.col.signal = vars)

# 7.d Compute Cluster Stability
start_time <- Sys.time()
clustRstab(data= dat, perturbedDataFun = perturbedDataFun1,
            kVec = kVec1,
            clCompScore = clCompScore1,
            clAlgo = clAlgo1,
            param = c(0,1),
            nsim = nsim1)

end_time <- Sys.time()
print(end_time - start_time)

