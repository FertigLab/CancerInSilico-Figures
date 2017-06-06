library(ggplot2)
library(methods)
library(CancerInSilico)
load('Figure_3_cleaned.RData') #fig3Data
load('tangData.RData') #tangData

# useful
l2norm <- function(a,b)
{
  sqrt(sum((a-b)^2))
}

## Figure 3a - fit real data w/out sync (days 2-7)

### fit no drug data

noDrugSims <- fig3Data[sapply(fig3Data, function(d) d$drugEffect==1)]
noDrugData <- subset(tangData, dosage==0 & day>1)$numCells

### fit drug data

timePoints <- seq(25,145,24)
l2 <- sapply(noDrugSims, function(sim)
    {
        d <- noDrugData * sim$numCells[timePoints[1]] / noDrugData[1]
        l2norm(sim$numCells[timePoints], d)
    })

noDrugFit <- noDrugSims[[which(l2==min(l2))]]
noDrugFit$numCells <- noDrugFit$numCells * noDrugData[1] / noDrugFit$numCells[25]

drugSims <- fig3Data[sapply(fig3Data, function(d)
    d$initDensity==noDrugFit$initDensity & d$cycleLength==noDrugFit$cycleLength)]

doses <- c(0.01, 0.1, 1, 10, 100)
drugFits <- lapply(doses, function(d)
    {
        drugData <- subset(tangData, dosage==d & day>1)$numCells
        l2 <- sapply(drugSims, function(sim)
            {
                d <- drugData * sim$numCells[timePoints[1]] / drugData[1]
                l2norm(sim$numCells[timePoints], d)
            })
        fit <- drugSims[[which(l2==min(l2))]]
        fit$numCells <- fit$numCells * drugData[1] / fit$numCells[25]
        return(fit)
    })

plotData <- data.frame(day=seq(2,7), real1=noDrugData, fit1=noDrugFit$numCells[timePoints],
    real2=subset(tangData, dosage==10 & day>1)$numCells, fit2=drugFits[[4]]$numCells[timePoints])

fig <- ggplot(plotData, aes(x=day)) + geom_point(aes(y=real1)) + geom_line(aes(y=fit1)) + 
    geom_point(aes(y=real2, color='10ug')) + geom_line(aes(y=fit2, color='10ug'))

ggsave(filename='fig3a.png', plot=fig)

#fig <- ggplot(tangData, aes(x=day, y=numCells, col=log(dosage), group=dosage)) + geom_line()
#ggsave(filename='rawData.png')


