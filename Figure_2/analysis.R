library(methods)
library(CancerInSilico)
load('Figure_2_cleaned.RData') #fig2Data

## Figure 2a - Growth curves, w/ and w/out boundary

png(filename='fig2a.png')

iDen <- which(unname(sapply(fig2Data, function(l) l$initDensity))==0.2)
iBound <- which(unname(sapply(fig2Data, function(l) l$boundary))>0)
iNoBound <- which(unname(sapply(fig2Data, function(l) l$boundary))==0)
noBound <- intersect(iDen, iNoBound)
bound <- intersect(iDen, iBound)

yMax <- max(sapply(noBound, function(i) max(fig2Data[[i]]$numCells)))
xMax <- length(fig2Data[[1]]$numCells)
plot(NULL, xlim=c(0,xMax), ylim=c(80,yMax), log='y')
invisible(sapply(noBound, function(i) lines(fig2Data[[i]]$numCells)))
invisible(sapply(bound, function(i) lines(fig2Data[[i]]$numCells)))
invisible(dev.off())

## Figure 2b


## Figure 2c - Density over time, w/ and w/out boundary

png(filename='fig2c.png')
iCL <- which(unname(sapply(fig2Data, function(l) l$cycleLength))==24)
iBound <- which(unname(sapply(fig2Data, function(l) l$boundary))>0)
iNoBound <- which(unname(sapply(fig2Data, function(l) l$boundary))==0)
noBound <- intersect(iCL, iNoBound)
bound <- intersect(iCL, iBound)

xMax <- length(fig2Data[[1]]$numCells)
plot(NULL, xlim=c(0,xMax), ylim=c(1,4))
invisible(sapply(noBound, function(i) lines(1/fig2Data[[i]]$contactInhibition)))
invisible(sapply(bound, function(i) lines(1/fig2Data[[i]]$contactInhibition)))
invisible(dev.off())

