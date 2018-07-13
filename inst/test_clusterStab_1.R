rm(list = ls())
library(clusterStability)

# Simulation cluster stability
data1 <- cbind(rnorm(n = 20, mean = 5, sd = 1), rnorm(n = 20, mean = 10, sd = 1))
data2 <- cbind(rnorm(n = 20, mean = 10, sd = 1), rnorm(n = 20, mean = 5, sd = 1))
data3 <- rbind(data1, data2)
head(data3)
tail(data3)

test <- instab.index2(data = data3, props = c(1,0.5), nb.grs = 1:10, score = "ARI", nsim = 50, mc.cores =  10, clustering.method = "Kmeans")
plot(1:10, test$instab.norm, type = "l")


test <- instab.index2(data = data3, props = c(1,0.5), nb.grs = 1:10, score = "NID", nsim = 50, mc.cores =  10, clustering.method = "Kmeans")
plot(1:10, test$instab.norm, type = "l")

# Il faut faire entrer les different paramteres dans les functions
# Il faut avoir un moyen de subsampler les obs.

