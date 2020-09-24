library(clustRstab)

clustRstab(data,
           perturbedDataFun = subSample,
           nsim = 2,
           kVec = 2:7,
           clAlgo = clAlgoHCWard,
           clCompScore = MARI,
           typeOfComp = "random",
           baseLineCorrection = FALSE,
           plot = TRUE)


############ NEW TEST ################

# Libraries
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

# xfit<-seq(0, max(geneSD),length=length(geneSD))
# yfit<-dnorm(xfit,mean=0.48,sd=0.15)
# lines(xfit, yfit, col="blue", lwd=1, lty = 2 )

# yfit<-dnorm(xfit,mean=0.8,sd=0.45)
# lines(xfit, yfit, col="red", lwd=1, lty = 2)
# abline(v = treshold, lty = 2, col = "blue")
# graphics::legend("topright", legend = c(" ", " "),
#                text.width = strwidth("1,000,000"),
#                lty = 1:2, xjust = 1, yjust = 1, inset = 1/10,
#                title = "Line Types", trace=TRUE)




# Define cancer type
type <- NCI60$type %>% as.tibble() %>% # A vector containing the cancer type labels
  filter(!grepl("repro", value) & !grepl("UNKNOWN", value)) %>%
  pull()


# Compute hclust dendogram
dend <- data %>%
  dist %>% hclust(., method = "ward.D2") %>%
  as.dendrogram

# Set cancer type as dendrogram labels
dendoLabelOrder <- labels(dend) # dendogram label order
labels(dend) <- type[dendoLabelOrder] # cancer type as dendogram labels
labelCols <- as.factor(labels(dend)) # color labels
levels(labelCols) <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
labels_colors(dend) <- as.character(labelCols)


## Get classifications
nbGrs <- 9 # Set the number of groups

HclustCl <- cutree(dend, k = nbGrs)[dendoLabelOrder] # HAC classification
KmeansCl <- as.factor(kmeans(data, nbGrs)$cluster[dendoLabelOrder]) # k-means classification
GmmCl <- as.factor(Mclust(data,
                          G = nbGrs,
                          modelNames = "EII")$classification[dendoLabelOrder]) # GMM classification

# Color classifications
grCols <- gray.colors(nbGrs)
levels(KmeansCl) <- grCols
levels(GmmCl) <- grCols

# Classification bars for plotting
the_bars <- cbind(as.character(KmeansCl), as.character(GmmCl))


# Plot the dendrogram with classification bars
par(mar = c(10,2,1,1))
dend %>% # set("branches_k_color", value = grCols, k=nbGrs) %>%
  set("leaves_pch", HclustCl) %>%
  plot
colored_bars(colors = the_bars, dend = dend, y_shift = -100, rowLabels = c("k-means", "GGM-EII"))

######
# Select the number of groups


# ----------------------------------------------------------------------------
# Clustering and selection of nb of grs
# ---------------------------------------------------------------------------


set.seed(31)

kMax <- 15

# Plot 1: Elbow, within cluster sum of squares, k-means
plot.elbow <- fviz_nbclust(data, kmeans,
                           nstart = 10,
                           method = "wss",
                           k.max = kMax) + theme_minimal() +
  ggtitle("the Elbow Method") + xlab("Number of clusters K")

# Plot 2: Gap-statistic, k-means
gap_stat <- clusGap(data, kmeans,
                    nstart = 10,
                    K.max = kMax,
                    B = 50)

plot.gap <- fviz_gap_stat(gap_stat) +
  theme_minimal() +
  ggtitle("Gap Statistic") +
  xlab("Number of clusters K") +
  ylab("Gap statistic")

# Plot for annex : silhouette, k-means
plot.silhouette <- fviz_nbclust(data, kmeans,
                                nstart = 30,
                                method = "silhouette",
                                k.max = kMax) +
  theme_minimal() + ggtitle("The Silhouette Plot") + xlab("Number of clusters K")



#-----------------------------------------------------
# Function to plot GMM
#-----------------------------------------------------
fviz_mclust_bic_01<- function (object, model.names = NULL, shape = 19, color = "model", palette = NULL, legend = NULL, main = "Model selection", xlab = "Number of components", ylab = "BIC", ...)
{ if (!inherits(object, "Mclust")) stop("An object of class Mclust is required.")
  best_model <- object$modelName
  number_of_cluster <- object$G
  x <- object$BIC
  n <- ncol(x)
  dnx <- dimnames(x)
  x <- matrix(as.vector(x), ncol = n)
  dimnames(x) <- dnx
  x <- as.data.frame(x)
  if (is.null(model.names)) model.names <- dimnames(x)[[2]]
  x <- x[, model.names, drop = FALSE]
  x <- cbind.data.frame(cluster = rownames(x), x)
  x <- tidyr::gather_(x, key_col = "model", value_col = "BIC", gather_cols = colnames(x)[-1])
  x <- x[!is.na(x$BIC), , drop = FALSE]
  x$cluster<- factor(x$cluster, levels = dnx[[1]])
  x$model <- factor(x$model, levels = dnx[[2]])
  if (ggpubr:::.is_col_palette(palette)) palette <- ggpubr:::.get_pal(palette, k = length(model.names))
  ggline.opts <- list(data = x, x = "cluster", y = "BIC", group = "model", color = color, shape = shape, palette = palette, main = main, xlab = xlab, ylab = ylab, ...)
  p <- do.call(ggpubr::ggline, ggline.opts) + #labs(subtitle = paste0("Best model: ", best_model, " | Optimal clusters: n = ", number_of_cluster)) +
    geom_vline(xintercept = number_of_cluster, linetype = 2, color = "#4E84C4") + theme(legend.title = element_blank())
  if (missing(legend)) p + theme(legend.position = c(0.7, 0.2), legend.direction = "horizontal", legend.key.height = unit(0.5, "line")) + guides(color = guide_legend(nrow = 5, byrow = TRUE)) else p + theme(legend.position = legend)
}
#-----------------------------------------------------

GMM <- Mclust(data, modelNames = c("EII", "EEI", "VEI"), G = 1:kMax)
# The best model is VEI when all models are studied

plot.gmm <- fviz_mclust_bic_01(GMM, model.names = NULL, shape = 19,
                               color = "model",
                               palette = c( "#52854C","#4E84C4", "#293352"),
                               legend = NULL,
                               main = "GMM - Model selection",
                               xlab = "Number of clusters K",
                               ylab = "BIC") +
  theme_minimal()  # +
#geom_vline(xintercept = 5, colour="#4E84C4", linetype="dashed")


## Cluster stability
# ---------------------------------------------------------------------------

data <- as.data.frame(data)

nPropVec <- c(0.65, 0.75, 0.85)

NCIclustRstabNprop <- mclapply(nPropVec, function(i) {
  clustRstab(data= data,
             perturbedDataFun = subSample,
             kVec = 2:15,
             typeOfComp = "all",
             clCompScore = aricode::MARI,
             clAlgo = clAlgoKmeans,
             nsim = 500, baseLineCorrection = F,
             plot = F, nProp = i, pProp = 1)
}, mc.cores = 8)


nPropVec <- c(0.65, 0.75, 0.85)

NCIclustRstabNPprop <- mclapply(nPropVec, function(i) {
  clustRstab(data= data,
             perturbedDataFun = subSample,
             kVec = 2:15,
             typeOfComp = "all",
             clCompScore = aricode::MARI,
             clAlgo = clAlgoKmeans,
             nsim = 500, baseLineCorrection = F,
             plot = F, pProp = i, nProp = i)
}, mc.cores = 8)

save(nPropVec, NCIclustRstabNprop, NCIclustRstabNprop, NCIclustRstabNPprop, file = "inst/NCIclustRstab.RData")




# Prepare plot with diff nProp
plot.clustRstab <- do.call(what = rbind, args = lapply(1:length(nPropVec),
                                           function(i) cbind(t(NCIclustRstabNprop[[i]][1:2,]),
                                                             K = 2:kMax,
                                                             nProp = rep(nPropVec[i], length(2:kMax))
                                                             ))) %>%
  as.tibble() %>%
  mutate(nProp = as.factor(nProp)) %>%
  ggplot() +
  geom_line(aes(x=K, y=mean,
                                   group=nProp, colour = nProp)) +
  geom_ribbon(aes(x = K,
                                     ymin=mean-sd,
                                     ymax=mean+sd,
                                     group=nProp), alpha=0.2) +
  # geom_vline(xintercept = nb.grs, color = "black", linetype="dashed") +
  scale_x_continuous(breaks = 2:kMax) +
  ggtitle(paste("Cluster stability")) +
  ylab("Mean clustRstab") +
  xlab("Number of clusters K") +
  scale_color_manual(values=c( "#52854C","#4E84C4", "#293352")) +
  ylim(-0.05, 1.05) +
  theme_minimal()


plot.clustRstab


#
# clustRstabPlot <- do.call(what = rbind, args = lapply(1:length(nPropVec),
#                                                       function(i) cbind(t(NCIclustRstabPprop[[i]][1:2,]),
#                                                                         K = 2:kMax,
#                                                                         nProp = rep(nPropVec[i], length(2:kMax))
#                                                       ))) %>%
#   as.tibble() %>%
#   mutate(nProp = as.factor(nProp)) %>%
#   ggplot() +
#   geom_line(aes(x=K, y=mean,
#                 group=nProp, colour = nProp)) +
#   geom_ribbon(aes(x = K,
#                   ymin=mean-sd,
#                   ymax=mean+sd,
#                   group=nProp), alpha=0.2) +
#   # geom_vline(xintercept = nb.grs, color = "black", linetype="dashed") +
#   scale_x_continuous(breaks = 2:kMax) +
#   ggtitle(paste("Cluster stability")) +
#   ylab("Mean clustRstab") +
#   xlab("Number of clusters K") +
#   scale_color_manual(values=c( "#52854C","#4E84C4", "#293352")) +
#   ylim(-0.05, 1.05) +
#   theme_minimal()
#
#
# clustRstabPlot



# ----------------------------------------------------------------------------
# The Plot
# ---------------------------------------------------------------------------
library("gridExtra")
library(cowplot)
plotGrSelect <- ggdraw() +
  draw_plot(plot.elbow, 0, .5, 0.33, .5) +
  draw_plot(plot.silhouette, .33, .5, .33, .5) +
  draw_plot(plot.gap, .66, .5, .33, .5) +
  draw_plot(plot.clustRstab, 0, 0, .5, .5) +
  draw_plot(plot.gmm, .5, 0, .5, .5) +
  draw_plot_label(c("A", "B", "C", "D", "E" ), c(0, 0.33, 0.66, 0, 0.5), c(1, 1, 1, 0.5, 0.5), size = 15)
plotGrSelect

# ----------------------------------------------------------------------------
# Evolution MARI avec les classif kmeans / hclust vs labels
# ---------------------------------------------------------------------------

library(aricode)
library(clustRstab)

MARItoPlot <- cbind(MARI = c(sapply(2:kMax, function(i) MARI(type, clAlgoHCWard(data, i))),
                    sapply(2:kMax, function(i) MARI(type, clAlgoKmeans(data, i))),
                    sapply(2:kMax, function(i) MARI(type, clAlgoGmmEII(data, i)))),
                    K = rep(2:kMax, 3),
                    clAlgo = c(rep("HCward", length(2:kMax)),
                               rep("k-means", length(2:kMax)),
                               c(rep("GMM-EII", length(2:kMax))))) %>%
  as.tibble() %>%
  mutate(MARI = as.numeric(MARI), K = as.numeric(K), clAlgo = as.factor(clAlgo)) %>%
  ggplot() +
  geom_line(aes(x=K, y=MARI,
                group=clAlgo, colour = clAlgo)) +
  scale_x_continuous(breaks = 2:kMax) +
  #ggtitle(paste("Cluster comparison for cancer type classification")) +
  ylab("MARI") +
  xlab("Number of clusters K") +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  #scale_color_manual(values=c( "#52854C","#4E84C4", "#293352")) +
  ylim(0.0, 1.0) +
  theme_minimal()


MARItoPlot




# cbind(MARI = c(sapply(2:kMax, function(i) NID(type, clAlgoHCWard(data, i))),
#                              sapply(2:kMax, function(i) NID(type, clAlgoKmeans(data, i))),
#                              sapply(2:kMax, function(i) NID(type, clAlgoGmmEII(data, i)))),
#                     K = rep(2:kMax, 3),
#                     clAlgo = c(rep("HCward", length(2:kMax)),
#                                rep("k-means", length(2:kMax)),
#                                c(rep("GMM-", length(2:kMax))))) %>%
#   as.tibble() %>%
#   mutate(MARI = as.numeric(MARI), K = as.numeric(K), clAlgo = as.factor(clAlgo)) %>%
#   ggplot() +
#   geom_line(aes(x=K, y=MARI,
#                 group=clAlgo, colour = clAlgo)) +
#   scale_x_continuous(breaks = 2:kMax) +
#   ggtitle(paste("Cluster comparison for cancer type classification")) +
#   ylab("MARI") +
#   xlab("Number of clusters K") +
#   #scale_color_manual(values=c( "#52854C","#4E84C4", "#293352")) +
#   ylim(0.0, 1.0) +
#   theme_minimal()



# ----------------------------------------------------------------------------
# Evolution MARI avec les classif kmeans / hclust vs labels
# ---------------------------------------------------------------------------

table(type, clAlgoKmeans(data, 9))
table(type, clAlgoKmeans(data, 4));table(type, clAlgoHCWard(data, 6))
kable(table(type, clAlgoHCWard(data, 4)))

clustRstab(data,
           perturbedDataFun = subSample,
           nsim = 500,
           kVec = 2:15,
           clAlgo = clAlgoKmeans,
           clCompScore = MARI,
           typeOfComp = "all",
           baseLineCorrection = FALSE,
           plot = TRUE)
subSample(data, nProp = 0.8, pProp = 0.8)
noiseGaussian(data, noiseGaussianMean = 0, noiseGaussianSD = 1)
randProjData(data, randProjDim = 20, randProjMethod = "Haar")

clAlgoGmmEEI(data, k = 3)
clAlgoHCWard(data, k = 3)
clAlgoKmeans(data, k = 3)


#----------------------------------------------------------------
#----------------------------------------------------------------
# Illustration of clustRstab package
#----------------------------------------------------------------
#----------------------------------------------------------------

#----------------------------------------------------------------
#----------------------------------------------------------------
# Data perturbation
#----------------------------------------------------------------
#----------------------------------------------------------------




funsPerturbData <- list(
  subSample = clustRstab::subSample,
  noiseGaussian = clustRstab::noiseGaussian,
  randProjData = clustRstab::randProjData
)


paramsPerturbData <- list(
  nPropVec <- c(0.65, 0.75, 0.85),
  params <- cbind(c(0,1,1), c(1,2,3)),
  randProjVec <- c(20, 100, 500)
)


NCIclustRstabDataPert <- mclapply(1:3, function(f){
 lapply(nPropVec, function(i) {
  clustRstab(data= data,
             perturbedDataFun = funsPerturbData[f],
             kVec = 2:kMax,
             typeOfComp = "all",
             clCompScore = aricode::MARI,
             clAlgo = clustRstab::clAlgoKmeans,
             nsim = nsim, baseLineCorrection = F,
             plot = F,
             nProp = paramsPerturbData[[f]][i], pProp = paramsPerturbData[[f]][i],
             noiseGaussianMean = paramsPerturbData[[f]][i,1], noiseGaussianSD = paramsPerturbData[[f]][i,2],
             randProjMethod = "Haar", randProjDim = paramsPerturbData[[f]][i]
             )
})
  }, mc.cores = mc.cores)

nPropVec <- c(0.65, 0.75, 0.85)

NCIclustRstabNprop <- mclapply(nPropVec, function(i) {
  clustRstab(data= data,
             perturbedDataFun = clustRstab::subSample,
             kVec = 2:kMax,
             typeOfComp = "all",
             clCompScore = aricode::MARI,
             clAlgo = clustRstab::clAlgoKmeans,
             nsim = nsim, baseLineCorrection = F,
             plot = F, nProp = i, pProp = i)
}, mc.cores = mc.cores)


#----------------------------------------------------------------
# prop
#----------------------------------------------------------------
kMax = 15
nsim = 500
mc.cores = 14


nPropVec <- c(0.65, 0.75, 0.85)

NCIclustRstabNprop <- mclapply(nPropVec, function(i) {
  clustRstab(data= data,
             perturbedDataFun = clustRstab::subSample,
             kVec = 2:kMax,
             typeOfComp = "all",
             clCompScore = aricode::MARI,
             clAlgo = clustRstab::clAlgoKmeans,
             nsim = nsim, baseLineCorrection = F,
             plot = F, nProp = i, pProp = i)
}, mc.cores = mc.cores)


#----------------------------------------------------------------
# adding Noise
#----------------------------------------------------------------


noiseMeanVec <- c(0,1,1)
noiseSDVec <- c(1,2,3)
params <- cbind(noiseMeanVec, noiseSDVec)


NCIclustRstabNoise <- mclapply(1:3, function(i) {
  clustRstab(data= data,
             perturbedDataFun = clustRstab::noiseGaussian,
             kVec = 2:kMax,
             typeOfComp = "all",
             clCompScore = aricode::MARI,
             clAlgo = clAlgoKmeans,
             nsim = nsim, baseLineCorrection = F,
             plot = F, noiseGaussianMean = params[i,1],
             noiseGaussianSD = params[i,2] )
}, mc.cores = mc.cores)





#----------------------------------------------------------------
# Random projection
#----------------------------------------------------------------


randProjVec <- c(10, 100, 500)

NCIclustRstabRandProj <- mclapply(1:3, function(i) {
  clustRstab(data= data,
             perturbedDataFun = clustRstab::randProjData,
             kVec = 2:kMax,
             typeOfComp = "all",
             clCompScore = aricode::MARI,
             clAlgo = clAlgoKmeans,
             nsim = nsim, baseLineCorrection = F,
             plot = F, randProjMethod = "Haar", randProjDim = randProjVec[i])
}, mc.cores = mc.cores)


#----------------------------------------------------------------
# Clustering Function
#----------------------------------------------------------------

funsClust <- list(
  kMeans = clustRstab::clAlgoKmeans,
  HAC_Ward = clustRstab::clAlgoHCWard,
  GMM = clustRstab::clAlgoGmmEII
)

NCIclustRstabClustRFun <- mclapply(funsClust, function(i) {
  clustRstab(data = data,
             perturbedDataFun = clustRstab::subSample,
             kVec = 2:kMax,
             typeOfComp = "all",
             clCompScore = aricode::MARI,
             clAlgo = i,
             nsim = nsim, baseLineCorrection = F,
             plot = F)
}, mc.cores = mc.cores)



#----------------------------------------------------------------
# Type of comparison
#----------------------------------------------------------------


comp <- c("all", "random", "toInitial")

NCIclustRstabComp <- mclapply(1:3, function(i) {
  clustRstab(data= data,
             perturbedDataFun = clustRstab::subSample,
             kVec = 2:kMax,
             typeOfComp = comp[i],
             clCompScore = aricode::MARI,
             clAlgo = clAlgoKmeans,
             nsim = nsim, baseLineCorrection = F,
             plot = F)
}, mc.cores = mc.cores)


#----------------------------------------------------------------
#----------------------------------------------------------------
# Type of score
#----------------------------------------------------------------
#----------------------------------------------------------------
funsScore <- list(
  MARI = aricode::MARI,
  NID = aricode::NID,
  NID = aricode::NID
)

adjustement <- c(FALSE, TRUE, FALSE)


NCIclustRstabScore <- mclapply(1:3, function(i) {
  clustRstab(data= data,
             perturbedDataFun = clustRstab::subSample,
             kVec = 2:kMax,
             typeOfComp = "all",
             clCompScore = funsScore[[i]],
             clAlgo = clAlgoKmeans,
             nsim = nsim, baseLineCorrection = adjustement[i],
             plot = F)
}, mc.cores = mc.cores)




#----------------------------------------------------------------
#----------------------------------------------------------------
# Plots
#----------------------------------------------------------------
#----------------------------------------------------------------


#---------------------------------------------------------------
# Que les pert data fun
#---------------------------------------------------------------


# 1
plot.clustRstabProp <- do.call(what = rbind, args = lapply(1:length(nPropVec),
                                                           function(i) cbind(t(NCIclustRstabNprop[[i]][1:2,]),
                                                                             K = 2:kMax,
                                                                             nProp = rep(nPropVec[i], length(2:kMax))
                                                           ))) %>%
  as.tibble() %>%
  mutate(nProp = as.factor(nProp)) %>%
  ggplot() +
  geom_line(aes(x=K, y=mean,
                group=nProp, colour = nProp)) +
  geom_ribbon(aes(x = K,
                  ymin=mean-sd,
                  ymax=mean+sd,
                  group=nProp), alpha=0.2) +
  # geom_vline(xintercept = nb.grs, color = "black", linetype="dashed") +
  scale_x_continuous(breaks = 2:kMax) +
  #ggtitle(paste("Subsampling observations")) +
  ylab("clustRstab (mean MARI)") +
  xlab("Number of groups K") +
  scale_color_manual(values=c("#52854C","#4E84C4", "#293352")) +
  ylim(-0.05, 1.05) +
  theme_minimal()+
  labs(colour = "prop") +theme(text = element_text(size = 17))
  #theme(legend.key.size = unit(1.32, "cm"),
   #     axis.title.x=element_blank(),
    #    axis.text.x=element_blank(),
     #   axis.ticks.x=element_blank())


#2
plot.clustRstabNoise <- do.call(what = rbind, args = lapply(1:length(noiseMeanVec),
                                                            function(i) cbind(t(NCIclustRstabNoise[[i]][1:2,]),
                                                                              K = 2:kMax,
                                                                              GaussianNoise = rep(i, length(2:kMax))
                                                            ))) %>%
  as.tibble() %>%
  mutate(GaussianNoise = factor(GaussianNoise, levels = 1:3, labels = c(paste('\u03bc',"=0, ", '\u03C3', "=1"), "\u03bc=1, \u03C3=2","\u03bc=1, \u03C3=3"))) %>%
  ggplot() +
  geom_line(aes(x=K, y=mean,
                group=GaussianNoise, colour = GaussianNoise)) +
  geom_ribbon(aes(x = K,
                  ymin=mean-sd,
                  ymax=mean+sd,
                  group=GaussianNoise), alpha=0.2) +
  # geom_vline(xintercept = nb.grs, color = "black", linetype="dashed") +
  scale_x_continuous(breaks = 2:kMax) +
  #ggtitle(paste("Adding Gaussian Noise")) +
  ylab(" ")+
  xlab("Number of groups K") +
  scale_color_manual(values=c("#52854C","#4E84C4", "#293352")) +
  ylim(-0.05, 1.05) +
  theme_minimal() +
   theme( text = element_text(size = 17),
     axis.title.y=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())+
  labs(colour = "noiseGaussian")

#3
plot.clustRstabRandProj <- do.call(what = rbind, args = lapply(1:length(noiseMeanVec),
                                                               function(i) cbind(t(NCIclustRstabRandProj[[i]][1:2,]),
                                                                                 K = 2:kMax,
                                                                                 RandProjDim = rep(i, length(2:kMax))
                                                               ))) %>%
  as.tibble() %>%
  mutate(RandProjDim = factor(RandProjDim, levels = 1:3, labels = c(10, 100, 500))) %>%
  ggplot() +
  geom_line(aes(x=K, y=mean,
                group=RandProjDim, colour = RandProjDim)) +
  geom_ribbon(aes(x = K,
                  ymin=mean-sd,
                  ymax=mean+sd,
                  group=RandProjDim), alpha=0.2) +
  # geom_vline(xintercept = nb.grs, color = "black", linetype="dashed") +
  ylab(" ") +
  scale_x_continuous(breaks = 2:kMax) +
  #ggtitle("Random Prjection") +
  xlab("Number of groups K") +
  scale_color_manual(values=c( "#52854C","#4E84C4", "#293352")) +
  ylim(-0.05, 1.05) +
  theme_minimal() +
   theme(text = element_text(size = 17),
     #legend.text = element_text(size = 13),
     axis.title.y=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())+
  labs(colour = "randProjDim")


library(ggpubr)
ggarrange(
  plot.clustRstabProp,
  plot.clustRstabNoise,
  plot.clustRstabRandProj,
  labels = c("A", "B", "C"),
  common.legend = FALSE, legend = "bottom", ncol = 3,
  widths = c(1.1,0.95,0.95)
)

#---------------------------------------------------------------
# The others
#---------------------------------------------------------------
#4
plot.clustRstabClust <- do.call(what = rbind, args = lapply(1:3,
                                                            function(i) cbind(t(NCIclustRstabClustRFun[[i]][1:2,]),
                                                                              K = 2:kMax,
                                                                              typeOfComp = rep(i, length(2:kMax))
                                                            ))) %>%
  as.tibble() %>%
  mutate(typeOfComp = factor(typeOfComp, levels = 1:3, labels = c("kmeans", "HAC", "GMM"))) %>%
  ggplot() +
  geom_line(aes(x=K, y=mean,
                group=typeOfComp, colour = typeOfComp)) +
  geom_ribbon(aes(x = K,
                  ymin=mean-sd,
                  ymax=mean+sd,
                  group=typeOfComp), alpha=0.2) +
  # geom_vline(xintercept = nb.grs, color = "black", linetype="dashed") +
  scale_x_continuous(breaks = 2:kMax) +
  # ggtitle("Different comparisons") +
  ylab("clustRstab (mean MARI)")+
  xlab("Number of clusters K") +
  scale_color_manual(values=c( "#52854C","#4E84C4", "#293352")) +
  ylim(-0.05, 1.05) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = "clAlgo")+
 theme(text = element_text(size = 17))



#5
plot.clustRstabComp <- do.call(what = rbind, args = lapply(1:length(noiseMeanVec),
                                                           function(i) cbind(t(NCIclustRstabComp[[i]][1:2,]),
                                                                             K = 2:kMax,
                                                                             typeOfComp = rep(i, length(2:kMax))
                                                           ))) %>%
  as.tibble() %>%
  mutate(typeOfComp = factor(typeOfComp, levels = 1:3, labels = c("all", "random", "toInitial"))) %>%
  ggplot() +
  geom_line(aes(x=K, y=mean,
                group=typeOfComp, colour = typeOfComp)) +
  geom_ribbon(aes(x = K,
                  ymin=mean-sd,
                  ymax=mean+sd,
                  group=typeOfComp), alpha=0.2) +
  # geom_vline(xintercept = nb.grs, color = "black", linetype="dashed") +
  scale_x_continuous(breaks = 2:kMax) +
  # ggtitle("Different comparisons") +
  ylab("clustRstab (mean MARI)") +
  xlab("Number of clusters K") +
  scale_color_manual(values=c( "#52854C","#4E84C4", "#293352")) +
  ylim(-0.05, 1.05) +
  theme_minimal() +
  theme(text = element_text(size = 17))+
  labs(colour = "typeOfComp")


# 6
plot.clustRstabScore <- do.call(what = rbind, args = lapply(1:3,
                                                               function(i) cbind(t(NCIclustRstabScore[[i]][1:2,]),
                                                                                 K = 2:kMax,
                                                                                 Score = rep(i, length(2:kMax))
                                                               ))) %>%
  as.tibble() %>%
  mutate(Score = factor(Score, levels = 1:3, labels = c("MARI", "NID adj.", "NID non-adj."))) %>%
  ggplot() +
  geom_line(aes(x=K, y=mean,
                group=Score, colour = Score)) +
  geom_ribbon(aes(x = K,
                  ymin=mean-sd,
                  ymax=mean+sd,
                  group=Score), alpha=0.2) +
  # geom_vline(xintercept = nb.grs, color = "black", linetype="dashed") +
  scale_x_continuous(breaks = 2:kMax) +
  #ggtitle("Different comparisons") +
  ylab("clustRstab (mean MARI or NID)")+
  xlab("Number of clusters K") +
  scale_color_manual(values=c("#52854C","#4E84C4", "#293352")) +
  ylim(-0.05, 1.05) +
  theme_minimal()+
  labs(colour = "clCompSc") +
  theme(text = element_text(size = 17))
 #legenWidht = 1.5

library(ggpubr)
ggarrange(
  plot.clustRstabClust,
  plot.clustRstabComp,
  plot.clustRstabScore,
  labels = c("A", "B", "C"),
  common.legend = FALSE, legend = "bottom", ncol = 3, widths = c(0.99,0.99,1.02)
)

