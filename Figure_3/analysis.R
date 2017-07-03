library(ggplot2)
library(methods)
library(CancerInSilico)
load('Figure_3_cleaned.RData') #fig3Data
load('tangData.RData') #tangData

# useful
l2norm <- function(a,b) sqrt(sum((a-b)^2))

## Figure 3 - fit real data w/out sync (days 2-7)

noDrugData <- subset(tangData, dosage==0)$numCells
timePoints <- seq(1,145,24)

### fit no drug data

noDrugSims <- fig3Data[sapply(fig3Data, function(d) d$drugEffect==1)]
l2 <- sapply(noDrugSims, function(sim)
    {
        nCells <- sim$numCells[timePoints]
        nCells <- nCells * noDrugData[2] / nCells[2]
        l2norm(nCells[2:7], noDrugData[2:7])
    })

noDrugFit <- noDrugSims[[which(l2==min(l2))]]
noDrugFit$numCells <- noDrugFit$numCells * noDrugData[2] /
    noDrugFit$numCells[timePoints[2]]

print(paste('init density:', noDrugFit$initDensity))
print(paste('mean cycle length:', noDrugFit$cycleLength))

### fit drug data

drugSims <- fig3Data[sapply(fig3Data, function(d)
    d$initDensity==noDrugFit$initDensity & d$cycleLength==noDrugFit$cycleLength)]

doses <- c(0.01, 0.1, 1, 10, 100)
drugFits <- lapply(doses, function(d)
    {
        drugData <- subset(tangData, dosage==d)$numCells
        ndx <- 2:7
        if (d == 1) ndx <- setdiff(2:7, 4)
        l2 <- sapply(drugSims, function(sim)
            {
                nCells <- sim$numCells[timePoints]
                nCells <- nCells  * drugData[2] / nCells[2]
                l2norm(nCells[ndx], drugData[ndx])
            })
        fit <- drugSims[[which(l2==min(l2))]]
        fit$numCells <- fit$numCells * drugData[2] / 
            fit$numCells[timePoints[2]]
        print(fit$drugEffect)
        return(fit)
    })

plotData <- data.frame(day=seq(1,7),
    noDrugReal=noDrugData, 
    noDrugSim=noDrugFit$numCells[timePoints],
    drugReal10=subset(tangData, dosage==10)$numCells,
    drugSim10=drugFits[[which(doses==10)]]$numCells[timePoints],
    drugReal100=subset(tangData, dosage==100)$numCells,
    drugSim100=drugFits[[which(doses==100)]]$numCells[timePoints])

fig <- ggplot(plotData, aes(x=day)) + 
    geom_point(aes(y=noDrugReal), size=2) + 
    geom_line(aes(y=noDrugSim)) + 
    geom_point(aes(y=drugReal10, color='10ug'), size=2) +
    geom_line(aes(y=drugSim10, color='10ug')) +
    geom_point(aes(y=drugReal100, color='100ug'), size=2) +
    geom_line(aes(y=drugSim100, color='100ug')) + 
    labs(title = "Fitting Real Data with CancerInSilico", 
        caption = "Figure 3", x = "Day", y = "Number Of Cells",
        color = "Dosage") 

ggsave(filename='fig3.png', plot=fig)

