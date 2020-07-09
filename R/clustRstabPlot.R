# library(reshape2)
# library(ggplot2)

#------------------------------------------------------------------------------
clustRstabPlot <- function(clustRstabObj, kVec, nsim){

# Graphical representation
  toPlot <- data.frame(cbind(kVec, t(clustRstabObj)))
  colnames(toPlot) <- c("nbGrs", "clustRstab", "sd" )

# graphical representation empirical variance
  clustRstabPlot <- ggplot(toPlot) +
    aes(x=nbGrs, y=clustRstab) +
    geom_ribbon(aes(ymin = clustRstab - sd,
                    ymax = clustRstab + sd),
                fill = "grey70") +
    geom_smooth(stat="identity") +
    ylab("Mean clustRstab") +
    xlab("K") +
    scale_x_continuous(breaks = kVec) +
    #ylim(c(-0.3,1.3)) +
    ggtitle(paste('clustRstab, nsim = ', nsim))


  return(clustRstabPlot)
}

