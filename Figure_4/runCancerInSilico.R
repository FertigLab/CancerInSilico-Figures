library(CancerInSilico)
library(foreach)
library(doParallel)

# load the fitted simulations from figure 3
load("../Figure_3/Fig3_Fitted_Simulations.RData")

# run cell model
runModel <- function(sim)
{
    type <- new('CellType', name='DEFAULT', minCycle=sim$cycleLength,
        cycleLength=function() sim$cycleLength)
    drug <- new('Drug', name='DEFAULT', timeAdded=24,
        cycleLengthEffect=function(a,b) rnorm(n=1, mean=b*sim$drugEffect, sd=4))
    inSilicoCellModel(initialNum=100, runTime=168, density=sim$initDensity,
        boundary=1, syncCycle=FALSE, randSeed=0, outputIncrement=4,
        recordIncrement=0.25, timeIncrement=0.001, cellTypes=c(type),
        cellTypeInitFreq=c(1), drugs=c(drug), maxDeformation=0.1,
        maxTranslation=0.1, maxRotation=0.3, nG=28, epsilon=10, delta=0.2)
}

# call runModel in parallel
cl <- makeCluster(6)
registerDoParallel(cl)
cellModels <- foreach(i = 1:9, .packages="CancerInSilico") %dopar%
{
    if (i <= 3)
    {
        runModel(tangFit[[i]])
    }
    else
    {
        runModel(kagoharaFit[[i-3]])
    }
}
stopCluster(cl)

names(cellModels) <- c(names(tangFit),names(kagoharaFit))
names(cellModels)

# save data
save(cellModels, file="Fig4_RawData.RData")
