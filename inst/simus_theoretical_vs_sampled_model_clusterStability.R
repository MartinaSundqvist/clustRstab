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
sim.sd <- rep(1, trueK)
#grSD <- rep(sim.sd, trueK)
vars <- 20
nsim = 500


# Cluster Stability Parammeters
Kmax <- 20
kVec <- 2:Kmax
clAlgo <- clAlgoKmeans
clCompScore = aricode::NID
typeOfComp <- "all"


#----------------------------------------------------------------
#----------------------------------------------------------------
# SIMULATIONS !!!
#----------------------------------------------------------------
#----------------------------------------------------------------



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


#------------------------------------------
#------------------------------------------
# PLOTS
#------------------------------------------
#------------------------------------------

#------------------------------------------
# Data representation
#------------------------------------------

# First case
p1 <-  dataList[[2]] %>%
    select(V1,V2) %>%
    mutate(V1 = as.numeric(V1),
           V2 = as.numeric(V2),
           groups = as.factor(rep(1:trueK, each = nPerGr))) %>%
  ggplot(aes(V1, V2, colour = groups)) +
    geom_point(alpha = 0.5) +
    scale_colour_brewer(palette = "Dark2") +
    theme_bw() +
    ylab("V2")+
    scale_y_continuous(breaks = seq(-8, 8, by=2))+
    scale_x_continuous(breaks = seq(-8, 8, by=2))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin =  unit(c(5.5,5.5,5.5,10), "pt")) +
  geom_vline(xintercept=seq(-6,6, by = 2), linetype="dotted")

p2 <-  dataList2[[2]] %>%
  select(V1,V2) %>%
  mutate(V1 = as.numeric(V1),
         V2 = as.numeric(V2),
         groups = as.factor(rep(1:trueK, each = nPerGr))) %>%
  ggplot(aes(V1, V2, colour = groups)) +
  geom_point(alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
  theme_bw() +
  ylab("V'2")+
  scale_y_continuous(breaks = seq(-8, 8, by=2))+
  scale_x_continuous(breaks = seq(-8, 8, by=2))+
  theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
        plot.margin =  unit(c(5.5,5.5,5.5,10), "pt"))+
  geom_vline(xintercept=c(seq(-6,2, by = 2),5,6), linetype="dotted")


# library(RColorBrewer)
# brewer.pal(n = 8, name = "Dark2") # To get colors !

p3 <-  dataList[[2]] %>%
  select(V1) %>%
  mutate(V1 = as.numeric(V1),
         groups = as.factor(rep(1:trueK, each = nPerGr))) %>%
  ggplot(aes(V1, color = groups, fill = groups)) +
  geom_density(alpha = 0.3) +
  scale_colour_brewer(palette = "Dark2") +
  theme_bw() +
  xlab("V1")+
  scale_fill_manual( values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")) +
  scale_x_continuous(breaks = seq(-8, 8, by=2)) +
  geom_vline(xintercept=seq(-6,6, by = 2), linetype="dotted")

p4 <-  dataList2[[2]] %>%
  select(V1) %>%
  mutate(V1 = as.numeric(V1),
         groups = as.factor(rep(1:trueK, each = nPerGr))) %>%
  ggplot(aes(V1, color = groups, fill = groups)) +
  geom_density(alpha = 0.3) +
  scale_colour_brewer(palette = "Dark2") +
  theme_bw() +
  xlab("V'1")+
  scale_fill_manual( values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")) +
  scale_x_continuous(breaks = seq(-8, 8, by=2)) +
  geom_vline(xintercept=c(seq(-6,2, by = 2),5,6), linetype="dotted")


#library(ggpubr)
ggarrange(p1, p2, p3, p4, ncol=2, nrow = 2, common.legend = TRUE, legend="right")

#----------------------------------------------------------------
# Global representation of simu study (with facet but not used for chapter)
#----------------------------------------------------------------
# Theoretical results
simu.clustRstabALLTheo <- rbind(
  cbind(t(rbind(nb.grs = kVec, simu1IM)), simuCase = rep(1, length(2:Kmax))),
  cbind(t(rbind(nb.grs = kVec, simu2IM)), simuCase = rep(2, length(2:Kmax)))) %>%
  as.tibble() %>%
  mutate(simuCase = factor(simuCase, levels = c(1, 2), labels = c("Symmetrical", "Unsymmetrical")))



# The plot with the sampled results
plot.clustRstabALL <- rbind(do.call(what = rbind,
                                    args = lapply(1:length(pPropVec),
                                                  function(i) cbind(t(simu1SM[[i]][1:2,]),
                                                                    K = 2:Kmax,
                                                                    nb.vars = rep(pPropVec[i]*vars, length(2:Kmax)),
                                                                    simuCase = rep(1, length(2:Kmax))
                                                  ))),
                            do.call(what = rbind,
                                    args = lapply(1:length(pPropVec),
                                                  function(i) cbind(t(simu2SM[[i]][1:2,]),
                                                                    K = 2:Kmax,
                                                                    nb.vars = rep(pPropVec[i]*vars, length(2:Kmax)),
                                                                    simuCase = rep(2, length(2:Kmax))
                                                  )))
) %>%
  as.tibble() %>%
  mutate(simuCase = factor(simuCase, levels = c(1, 2), labels = c("Symmetrical", "Unsymmetrical"))) %>%
  ggplot() +
  geom_line(aes(x=K, y=mean,
                group=nb.vars, colour = nb.vars)) +
  geom_ribbon(aes(x = K,
                  ymin=mean-sd,
                  ymax=mean+sd,
                  group=nb.vars), alpha=0.1) +
  geom_vline(xintercept = trueK, color = "black", linetype="dashed") +
  geom_line(data = simu.clustRstabALLTheo, aes(x=nb.grs, y=mean), colour = "red", size = 1.2) +
  geom_ribbon(data = simu.clustRstabALLTheo, aes(x = nb.grs,
                                                 ymin=mean-sd,
                                                 ymax=mean+sd), alpha=0.01, colour = "red", linetype = "dotted") +
  scale_x_continuous(breaks = 2:Kmax) +
  # ggtitle(paste("Cluster stability")) +
  ylab("Cluster Stability (mean NID)") +
  xlab("Number of clusters K") +
  # scale_color_manual(values=c( "#52854C","#4E84C4", "#293352")) +
  ylim(-0.05, 0.5) +
  theme_bw()

plot.clustRstabALL + facet_grid(.~ simuCase)
