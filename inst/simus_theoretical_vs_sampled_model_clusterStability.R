library(tidyverse)
library(parallel)

#----------------------------------------------------------------
# Get data function
#----------------------------------------------------------------

getMultiNormalData <- function(i, nbGrs, nPerGr, muVec, sdVec,
                               nColSignal){ #nColNoise = 0
  g <- list()
  for (gr in 1:nbGrs){
    g[[gr]] <- replicate(n = nColSignal, rnorm(n=nPerGr, mean = muVec[gr], sd = sdVec[gr]))
  }

  mnData <- do.call(rbind, g)
#  mnData <- cbind(mnData, replicate(n.col.noise, rnorm(nb.grs*obsPerGr, 0,1)))

  as.data.frame(mnData)
}


#----------------------------------------------------------------
# Get dataList for theoretical clsuterstab
#----------------------------------------------------------------

getNsimDataSets <- function(nbGrs, nPerGr, muVec, sdVec,
                            nColSignal, nsim){
  dataList <- lapply(1:nsim, function(i) getMultiNormalData(i, nbGrs, nPerGr, muVec, sdVec,
                                                                     nColSignal))
  dataList
}


#----------------------------------------------------------------
# First simulation
#----------------------------------------------------------------

trueK <- 7
nPerGr <- 50
grMeans <- c(-6,-4,-2,0,2,4,6)
grMeans2 <- c(-6,-4,-2,0,2,5,6)
sim.sd <- rep(0.7, trueK)
#grSD <- rep(sim.sd, trueK)
vars <- 20
nsim = 500

clAlgoKmeans2 <- function(data, k){
  cl <- kmeans(x = data, centers = k)$cluster
  cl
}

# Cluster Stability Parammeters
Kmax <- 20
kVec <- 2:Kmax
clAlgo <- clAlgoKmeans
clCompScore = aricode::NID
typeOfComp <- "all"


#----------------------------------------------------------------
# Theoretical cluster stability
#----------------------------------------------------------------
source("inst/TheoClustRstab.R")
set.seed(67)
dataList <- getNsimDataSets(nbGrs = trueK,
                            nPerGr = nPerGr,
                            muVec = grMeans,
                            sdVec = sim.sd,
                            nColSignal = vars,
                            nsim = nsim)

simu1IM <- TheoClustRstab(dataList = dataList,
                          kVec = kVec,
                          typeOfComp = typeOfComp,
                          clAlgo = clAlgo,
                          clCompScore = clCompScore,
                          nsim = nsim,
                          plot = FALSE)
# Simu2
dataList2 <- getNsimDataSets(nbGrs = trueK,
                            nPerGr = nPerGr,
                            muVec = grMeans2,
                            sdVec = sim.sd,
                            nColSignal = vars,
                            nsim = nsim)

simu2IM <- TheoClustRstab(dataList = dataList2,
                          kVec = kVec,
                          typeOfComp = typeOfComp,
                          clAlgo = clAlgo,
                          clCompScore = clCompScore,
                          nsim = nsim,
                          plot = FALSE)


# First case


#----------------------------------------------------------------
# Subsampled cluster stability
#----------------------------------------------------------------


# First case


data <- getMultiNormalData(nbGrs = trueK,
                           nPerGr = nPerGr,
                            muVec = grMeans,
                           sdVec =  sim.sd,
                           nColSignal =  vars)

pPropVec = seq(0.25, 0.95, 0.1)

simu1SM <-  mclapply(pPropVec, function(p)
  clustRstab(data,
             kVec,
             typeOfComp = "all",
             clAlgo = clAlgoKmeans,
             nsim = nsim,
             clCompScore = aricode::NID,
             nProp = 1,
             pProp = p,
             perturbedDataFun = subSample),
  mc.cores = 14)


# PLOT !

simu1IMPlot <- as.data.frame(t(rbind(nb.grs = kVec, simu1IM)))


plot.clustRstab1 <- do.call(what = rbind,
                            args = lapply(1:length(pPropVec),
                                          function(i) cbind(t(simu1SM[[i]][1:2,]),
                                                            K = 2:Kmax,
                                                            pProp = rep(pPropVec[i], length(2:Kmax))
                                          ))) %>%
  as.tibble() %>%
  # mutate(pProp = as.factor(pProp)) %>%
  ggplot() +
  geom_line(aes(x=K, y=mean,
                group=pProp, colour = pProp)) +
  geom_ribbon(aes(x = K,
                  ymin=mean-sd,
                  ymax=mean+sd,
                  group=pProp), alpha=0.1) +
  geom_vline(xintercept = trueK, color = "black", linetype="dashed") +
  geom_line(data = simu1IMPlot, aes(x=nb.grs, y=mean), colour = "red", size = 1.2) +
  geom_ribbon(data = simu1IMPlot, aes(x = nb.grs,
                                      ymin=mean-sd,
                                      ymax=mean+sd), alpha=0.01, colour = "red", linetype = "dotted") +
  scale_x_continuous(breaks = 2:Kmax) +
  # ggtitle(paste("Cluster stability")) +
  ylab("Cluster Stability (mean NID)") +
  xlab("Number of clusters K") +
  # scale_color_manual(values=c( "#52854C","#4E84C4", "#293352")) +
  ylim(-0.05, 0.4) +
  theme_bw()

print(plot.clustRstab1)


# Second case


data <- getMultiNormalData(nbGrs = trueK,
                           nPerGr = nPerGr,
                           muVec = grMeans2,
                           sdVec =  sim.sd,
                           nColSignal =  vars)

simu2SM <-  mclapply(pPropVec, function(p)
  clustRstab(data,
             kVec,
             typeOfComp = "all",
             clAlgo = clustRstab::clAlgoKmeans,
             nsim = nsim,
             clCompScore = aricode::NID,
             nProp = 1,
             pProp = p,
             perturbedDataFun = subSample),
  mc.cores = 14)


# PLOT !

simu2IMPlot <- as.data.frame(t(rbind(nb.grs = kVec, simu2IM)))


plot.clustRstab2 <- do.call(what = rbind,
                           args = lapply(1:length(pPropVec),
                                         function(i) cbind(t(simu2SM[[i]][1:2,]),
                                                           K = 2:Kmax,
                                                           pProp = rep(pPropVec[i], length(2:Kmax))
                                         ))) %>%
  as.tibble() %>%
  # mutate(pProp = as.factor(pProp)) %>%
  ggplot() +
  geom_line(aes(x=K, y=mean,
                group=pProp, colour = pProp)) +
  geom_ribbon(aes(x = K,
                  ymin=mean-sd,
                  ymax=mean+sd,
                  group=pProp), alpha=0.1) +
  geom_vline(xintercept = trueK, color = "black", linetype="dashed") +
  geom_line(data = simu2IMPlot, aes(x=nb.grs, y=mean), colour = "red", size = 1.2) +
  geom_ribbon(data = simu2IMPlot, aes(x = nb.grs,
                  ymin=mean-sd,
                  ymax=mean+sd), alpha=0.01, colour = "red", linetype = "dotted") +
  scale_x_continuous(breaks = 2:Kmax) +
  # ggtitle(paste("Cluster stability")) +
  ylab("Cluster Stability (mean NID)") +
  xlab("Number of clusters K") +
  # scale_color_manual(values=c( "#52854C","#4E84C4", "#293352")) +
  ylim(-0.05, 0.4) +
  theme_bw()


print(plot.clustRstab2)



#------------------------------------------
# Data representation
#------------------------------------------

# First case
dataReptoPlot <- dataList[[1]] %>%
  select(V1,V2) %>%
  mutate(V1 = as.numeric(V1), V2 = as.numeric(V2), groups = as.factor(rep(1:trueK, each =nPerGr )))
glimpse(dataReptoPlot)

ggplot(dataReptoPlot, aes(V1, V2, colour = groups)) +
  geom_point(alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
  theme_bw() +
  stat_ellipse()


dataReptoPlot2 <- dataList2[[1]] %>%
  select(V1,V2) %>%
  mutate(V1 = as.numeric(V1), V2 = as.numeric(V2), groups = as.factor(rep(1:trueK, each =nPerGr )))
glimpse(dataReptoPlot)

ggplot(dataReptoPlot2, aes(V1, V2, colour = groups)) +
  geom_point(alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
  theme_bw()+
  stat_ellipse()

