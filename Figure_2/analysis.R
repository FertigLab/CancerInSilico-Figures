library(ggplot2)
library(methods)
library(CancerInSilico)
load('Figure_2_cleaned.RData') #fig2Data

## Figure 2a - Growth curves, w/ and w/out boundary

rawData <- lapply(fig2Data, function(mod)
    {
        N <- length(mod$numCells)
        timePoints <- 1:N - 1
        if (mod$initDensity == 0.10 & (mod$cycleLength %in% c(12,20,28,36,44)))
        {
            cbind(timePoints, mod$numCells, rep(mod$cycleLength, N),
                rep(ifelse(mod$boundary > 0, TRUE, FALSE), N))
        }
    })
rawData <- rawData[!sapply(rawData, is.null)]
rawData <- do.call(rbind, rawData)
df <- data.frame(time=rawData[,1], nCells=rawData[,2], cycLength=rawData[,3],
    bd=rawData[,4]>0)

fig <- ggplot(df, aes(x=time, y=nCells, linetype=bd,
    group=interaction(cycLength, bd))) + geom_line() +
    scale_linetype_manual(values=c("dashed", "solid")) +
    scale_y_continuous(trans='log10', breaks=c(100,200,500,1000)) +
    labs(title = "Sensitivity to Growth Rate and Boundary Presence", 
        caption = "Figure 2a", x = "Time", y = "Number Of Cells (log scale)",
        linetype = "Boundary") 

ggsave(filename='fig2a.pdf', plot=fig)

## Figure 2b - Density over time, w/ and w/out boundary

# smooths out line
movAvg <- function(data, window)
{
    avg <- c()    
    for (i in 1:length(data))
    {
        mn <- max(1, i - window)
        mx <- min(length(data), i + window)
        avg <- c(avg, mean(data[mn:mx]))
    }
    return(avg)
}

rawData <- lapply(fig2Data, function(mod)
    {
        N <- length(mod$density)
        timePoints <- 1:N - 1
        if (mod$cycleLength == 24 & (mod$initDensity %in% c(0.05,0.1,0.2,0.3,0.4)))
        {
            cbind(timePoints, movAvg(mod$density, 4), rep(mod$initDensity, N),
                rep(ifelse(mod$boundary > 0, TRUE, FALSE), N))
        }
    })
rawData <- rawData[!sapply(rawData, is.null)]
rawData <- do.call(rbind, rawData)
df <- data.frame(time=rawData[,1], den=rawData[,2], iDen=rawData[,3],
    bd=rawData[,4]>0)

fig <- ggplot(df, aes(x=time, y=den, linetype=bd,
    group=interaction(iDen, bd))) + geom_line() +
    scale_linetype_manual(values=c("dashed", "solid")) +
    labs(title = "Boundary Presence Effects Maximum Population Density", 
        caption = "Figure 2b", x = "Time", y = "Population Density",
        linetype = "Boundary")

ggsave(filename='fig2b.pdf', plot=fig)

