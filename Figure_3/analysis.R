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

noDrugData <- subset(tangData, dosage==0 & day>1)$numCells
noDrugData <- noDrugData * 80 / noDrugData[1]
noDrugSim <- fig3Data[sapply(fig3Data, function(d) !d$synced & d$drugEffect==1)]

t <- seq(1,121,24)
l2 <- sapply(noDrugSim, function(sim) l2norm(sim$numCells[t], noDrugData))
noDrugFit <- noDrugSim[[which(l2==min(l2))]]

validSim <- fig3Data[sapply(fig3Data, function(d) !d$synced &
    d$initDensity==noDrugFit$initDensity & d$cycleLength==noDrugFit$cycleLength)]

doses <- c(0.01, 0.1, 1, 10, 100)
fit <- lapply(doses, function(d)
    {
        noDrugData <- subset(tangData, dosage==d & day>1)$numCells
        noDrugData <- noDrugData * 80 / noDrugData[1]
        l2 <- sapply(validSim, function(sim) l2norm(sim$numCells[t], noDrugData))
        return(validSim[[which(l2==min(l2))]])
    })

plotData <- data.frame(day=seq(2,7), real1=noDrugData, fit1=noDrugFit$numCells[t],
    real2=subset(tangData, dosage==10 & day>1)$numCells, fit2=fit[[1]]$numCells[t])
fig <- ggplot(plotData, aes(x=day)) + geom_point(aes(y=real1)) + geom_line(aes(y=fit1)) +
    geom_point(aes(y=real2, color='red')) + geom_line(aes(y=fit2))
ggsave(filename='fig3a.png', plot=fig)

#fig <- ggplot(tangData, aes(x=day, y=numCells, col=log(dosage), group=dosage)) + geom_line()
#ggsave(filename='rawData.png')


