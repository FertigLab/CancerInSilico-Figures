library(ggplot2)
library(methods)
library(CancerInSilico)
load('Figure_3_cleaned.RData') #fig3Data

len <- length(fig2Data[[1]]$numCells)
rawData <- lapply(fig2Data, function(mod)
    {
        cbind(1:len-1, mod$numCells,
            rep(mod$synced, len), rep(mod$initDensity, len),
            rep(mod$drugEffect, len), rep(mod$cycleLength, len))
    })
rawData <- do.call(rbind, rawData)
data <- data.frame(time=rawData[,1], numCells=rawData[,2],
    synced=rawData[,3]>0, initDen=rawData[,4], drugEffect=1/rawData[,5]/2,
    cycleLength=rawData[,6])

## Figure 3a - fit real data w/out sync

fig <- ggplot(subset(data, initDen==0.05 & cycleLength==24), 
    aes(x=time, y=numCells, col=drugEffect)) + geom_point()
ggsave(filename='fig3a.png', plot=fig)

## Figure 3b - fit real data w/ sync
