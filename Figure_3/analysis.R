library(ggplot2)
library(methods)
library(CancerInSilico)
load('Figure_3_cleaned.RData') #fig3Data
load('tangData.RData') #tangData

# useful
l2norm <- function(a,b) sqrt(sum((a-b)^2))

## Figure 3 - fit real data w/out sync (days 2-7)

noDrugData <- subset(tangData, dosage==0 & day>1)$numCells
timePoints <- seq(24,144,24)

### fit no drug data

noDrugSims <- fig3Data[sapply(fig3Data, function(d) d$drugEffect==1)]
l2 <- sapply(noDrugSims, function(sim)
    {
        nCells <- sim$numCells[timePoints]
        nCells <- nCells * noDrugData[1] / nCells[1]
        l2norm(nCells, noDrugData)
    })

noDrugFit <- noDrugSims[[which(l2==min(l2))]]
noDrugFit$numCells <- noDrugFit$numCells * noDrugData[1] /
    noDrugFit$numCells[timePoints[1]]

print(paste('init density:', noDrugFit$initDensity))
print(paste('mean cycle length:', noDrugFit$cycleLength))

### fit drug data

drugSims <- fig3Data[sapply(fig3Data, function(d)
    d$initDensity==noDrugFit$initDensity & d$cycleLength==noDrugFit$cycleLength)]

doses <- c(0.01, 0.1, 1, 10, 100)
drugFits <- lapply(doses, function(d)
    {
        drugData <- subset(tangData, dosage==d & day>1)$numCells
        l2 <- sapply(drugSims, function(sim)
            {
                nCells <- sim$numCells[timePoints]
                nCells <- nCells  * drugData[1] / nCells[1]
                l2norm(nCells, drugData)
            })
        fit <- drugSims[[which(l2==min(l2))]]
        fit$numCells <- fit$numCells * drugData[1] / 
            fit$numCells[timePoints[1]]
        return(fit)
    })

plotData <- data.frame(day=seq(2,7),
    noDrugReal=noDrugData, 
    noDrugSim=noDrugFit$numCells[timePoints],
    drugReal=subset(tangData, dosage==10 & day>1)$numCells,
    drugSim=drugFits[[which(doses==10)]]$numCells[timePoints])

fig <- ggplot(plotData, aes(x=day)) + 
    geom_point(aes(y=noDrugReal)) + 
    geom_line(aes(y=noDrugSim)) + 
    geom_point(aes(y=drugReal, color='10ug')) +
    geom_line(aes(y=drugSim, color='10ug'))

ggsave(filename='fig3.png', plot=fig)

