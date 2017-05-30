library(ggplot2)
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
plot.new()
plot.window(xlim=c(0,xMax), ylim=c(80,yMax), log='y')
axis(side=1)
axis(side=2)
invisible(sapply(noBound, function(i) lines(fig2Data[[i]]$numCells)))
invisible(sapply(bound, function(i) lines(fig2Data[[i]]$numCells)))
invisible(dev.off())

## Figure 2b - Contact Inhibition vs Density

rawData <- sapply(fig2Data, function(mod)
    {
        if (mod$initDensity == 0.05)
            cbind(mod$density, 1/mod$contactInhibition,
            rep(ifelse(mod$boundary > 0, T, F), length(mod$density)))
    })

rawData <- rawData[!sapply(rawData, is.null)]
rawData <- do.call(rbind, rawData)
df <- data.frame(den=rawData[,1], ci=rawData[,2], bd=rawData[,3]>0)
df$den <- round(df$den, 2)
df <- df[df$den < 0.7 | !df$bd,]
df <- df[df$den < 0.55 | df$bd,]

df$ci[df$bd] <- sapply(df$den[df$bd], function(d) mean(df$ci[df$den==d & df$bd]))
df$ci[!df$bd] <- sapply(df$den[!df$bd], function(d) mean(df$ci[df$den==d & !df$bd]))

fig <- ggplot(df, aes(x=den, y=ci, color=bd)) + geom_line()
ggsave(filename='fig2b.png', plot=fig)

## Figure 2c - Density over time, w/ and w/out boundary

png(filename='fig2c.png')

iCL <- which(unname(sapply(fig2Data, function(l) l$cycleLength))==24)
iBound <- which(unname(sapply(fig2Data, function(l) l$boundary))>0)
iNoBound <- which(unname(sapply(fig2Data, function(l) l$boundary))==0)
noBound <- intersect(iCL, iNoBound)
bound <- intersect(iCL, iBound)

xMax <- length(fig2Data[[1]]$numCells)
plot.new()
plot.window(xlim=c(0,xMax), ylim=c(0,1))
invisible(sapply(noBound, function(i) lines(fig2Data[[i]]$density)))
invisible(sapply(bound, function(i) lines(fig2Data[[i]]$density)))
invisible(dev.off())


