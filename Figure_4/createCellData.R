library(CancerInSilico)
library(foreach)
library(doParallel)

# load the fitted simulations from figure 3
load("../Figure_3/Fig3_Fitted_Simulations.RData")

# run cell model
runModel <- function(sim, time, addVariance)
{
    if (addVariance)
    {
        cellTypes <- c(new('CellType', name='DEFAULT', minCycle=sim$cycleLength-6,
            cycleLength=function() sim$cycleLength - 6 + rexp(1,1/6)))
    }
    else
    {
        cellTypes <- c(new('CellType', name='DEFAULT', minCycle=sim$cycleLength,
            cycleLength=function() sim$cycleLength))
    }

    if (sim$drugEffect == 1)
    {
        drugs <- list()
    }
    else
    {
        drugs <- c(new('Drug', name='DEFAULT', timeAdded=0,
        cycleLengthEffect=function(a,b) rnorm(n=1, mean=b*sim$drugEffect, sd=4)))
    }

    inSilicoCellModel(initialNum=60, runTime=time, density=sim$initDensity,
        boundary=1, syncCycle=FALSE, randSeed=0, outputIncrement=4,
        recordIncrement=0.25, timeIncrement=0.001, cellTypes=cellTypes,
        cellTypeInitFreq=c(1), drugs=drugs, maxTranslation=0.1,
        maxRotation=0.3, nG=28, epsilon=10, delta=0.2)
}

# call runModel in parallel
cl <- makeCluster(6)
registerDoParallel(cl)
cellModels <- foreach(i = 1:16, .packages="CancerInSilico") %dopar%
{
    if (i <= 4)
    {
        j <- ifelse(i %in% c(2,4), 3, 1)
        addVar <- (i > 2)
        runModel(tangFit[[j]], 168, addVar)
    }
    else
    {
        j <- ifelse(i > 10, i - 10, i - 4)
        addVar <- (i > 10)
        runModel(kagoharaFit[[j]], 144, addVar)
    }
}
stopCluster(cl)

names(cellModels) <- c(rep(names(tangFit[c(1,3)]), 2), rep(names(kagoharaFit), 2))

# save data
save(cellModels, file="Fig4_RawData.RData")
