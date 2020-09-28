# Libraries for cluseter stability
library(clustRstab)
library(clv)
library(clusterStab)
library(ClusterStability)
library(clValid)
library(ConsensusClusterPlus)



# Other libraries
library(tidyverse)
library(dendextend) # Library for plotting dendogram
library(grDevices) # Library for gray colors
library(mclust) # Library for copmputing gaussian mixture model (GMM)
library(factoextra)
library(cluster)
library(parallel)
library(aricode)

## Prepare data
# Define data
data(NCI60)

data <- t(NCI60$expr) %>% # A 59 x 6830 dataframe with gene expression for the 59 of the cancer cell lines
  as.tibble()  %>%
  mutate(type = NCI60$type) %>%
  filter( !grepl("repro", type), !grepl("UNKNOWN", type)) %>%
  select(-type)

# Select genes
geneSD <- apply(data, 2, sd)
geneSelect <- geneSD %>%
  density %>% plot(., main = 'Desity plot', xlab = "Standard Deviation", ylab = "density" , ylim = c(0,2.7))
treshold <- 0.8
data <- data[, geneSD > treshold]
dim(data)
data

seedNb <- 67432
set.seed(seedNb)

#-------------------------------------------------------------------
# clusterStab
#-------------------------------------------------------------------


bh.time <- system.time(clusterStab::benhur(as.matrix(data), # 6 times slower than my function when clustRstab is not parallellized...
                          freq=0.75, # Prop of subsampled obs
                          upper = 9, # Kmax
                          seednum=seedNb, # seed for
                          linkmeth = "ward",
                          distmeth = "euclidean",
                          iterations = 100))


bh <- clusterStab::benhur(as.matrix(data),
                                           freq=0.75, # Prop of subsampled obs
                                           upper = 9, # Kmax
                                           seednum=seedNb, # seed for
                                           linkmeth = "ward",
                                           distmeth = "euclidean",
                                           iterations = 100)


# Comparison is done to initial clustering
# Do figure at the office
# Parallellize the clustRstab fun !


hist(bh)


# Function to "validate" clusterings, not sure how it works or how to interpret the resutls cuz only at 0 all the time...
clusterCOmpNCI602 <- clusterStab::clusterComp(as.matrix(data),
                                              cl = 2,
                                              seednum = seedNb,
                                              B =10,
                                              sub.frac= 0.75, # prop for subsampling varaibles
                                              method="ward",
                                              distmeth = "euclidean")

clusterCOmpNCI603 <- clusterStab::clusterComp(as.matrix(data),
                                                cl = 3,
                                                seednum = seedNb,
                                                B = 10,
                                                sub.frac= 0.75,
                                                method="ward",
                                                distmeth = "euclidean")

#-------------------------------------------------------------------
# clv
#-------------------------------------------------------------------

cls.stab.sim.indNCI60 <- cls.stab.sim.ind(data,
                                 cl.num = 2:15,
                                 rep.num=99,
                                 subset.ratio = 0.75,
                                 clust.method = c("kmeans", "hclust"),
                                 method.type = "ward",
                                 sim.ind.type = c("rand", "jaccard","sim.ind", "dot.pr"))


mean.kmeans.rand.clv <- apply(cls.stab.sim.indNCI60$kmeans$rand, 2, mean)
sd.kmeans.rand.clv <- apply(cls.stab.sim.indNCI60$kmeans$rand, 2, sd)

mean.kmeans.jaccard.clv <- apply(cls.stab.sim.indNCI60$kmeans$jaccard, 2, mean)
sd.kmeans.jaccard.clv <- apply(cls.stab.sim.indNCI60$kmeans$jaccard, 2, sd)


mean.hclust.rand.clv <- apply(cls.stab.sim.indNCI60$hclust.ward$rand, 2, mean)
sd.hclust.rand.clv <- apply(cls.stab.sim.indNCI60$hclust.ward$rand, 2, sd)

mean.hclust.jaccard.clv <- apply(cls.stab.sim.indNCI60$hclust.ward$jaccard, 2, mean)
sd.hclust.jaccard.clv <- apply(cls.stab.sim.indNCI60$hclust.ward$jaccard, 2, sd)


clvRes <- cbind(mean= c(mean.kmeans.rand.clv, mean.kmeans.jaccard.clv, mean.hclust.rand.clv,mean.hclust.jaccard.clv),
                sd = c(sd.kmeans.rand.clv, sd.kmeans.jaccard.clv, sd.hclust.rand.clv, sd.hclust.jaccard.clv),
                clAlgo= rep(1:2, each = 14*2), score = rep(rep(1:2, each = 14),2),
                K = rep(2:15, 4)) %>%
  as.tibble() %>%
  mutate(clAlgo = factor(clAlgo, levels = 1:2, labels = c("kmeans","hclust")),
         score = factor(score, levels = 1:2, labels = c("rand","jaccard"))) %>%
  ggplot() +
  ggplot2::geom_line(aes(x = K, y = mean, lty = clAlgo, color = score))+
  geom_ribbon(aes(x = K,
                  ymin=mean-sd,
                  ymax=mean+sd,
                  group=interaction(score, clAlgo)), alpha=0.2) +
  scale_x_continuous(breaks = 2:15) +
  #ggtitle("Different comparisons") +
  ylab("clv - mean Rand or Jaccard")+
  xlab("Number of clusters K") +
  scale_color_manual(values=c("#52854C","#4E84C4")) +
  ylim(-0.05, 1.05) +
  theme_minimal()

clvRes



#-------------------------------------------------------------------
# ClusterStability
#-------------------------------------------------------------------
devtools::install_version("latticeExtra",repos = "https://www.stats.bris.ac.uk/R/ ", version="0.6-28")


#-------------------------------------------------------------------
# Consensus Clustering
#-------------------------------------------------------------------

ConsensusClusterPlus::ConsensusClusterPlus(d = as.matrix(data),
                                           maxK = 15,
                                           reps = 100,
                                           pItem = 0.75,
                                           pFeature = 0.75,
                                           clusterAlg = "km",
                                           title = "NCI60ConsensusClustering",
                                           seed = seedNb)
