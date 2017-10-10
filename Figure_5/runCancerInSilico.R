library('CancerInSilico')
library(methods)

args <- commandArgs(TRUE)
arrayNum <- as.integer(args[1])
returnSize <- as.integer(args[2])

#### Set Defaults ####

initialNum <- 100
runTime <- 168
density <- 0.05
boundary <- 1
syncCycles <- FALSE
randSeed <- 0
outputIncrement <- 4
recordIncrement <- 0.25
timeIncrement <- 0.001
cellTypes <- c(new('CellType', name='DEFAULT'))
cellTypeInitFreq <- c(1)
drugs <- list()
maxDeformation <- 0.1
maxTranslation <- 0.1
maxRotation <- 0.3
nG <- 28
epsilon <- 10.0
delta <- 0.2

#### Set Custom Values ####

type1 <- new('CellType', name='A', minCycle=10-2,
    cycleLength=function() max(10-2, rnorm(1,10,1)))
type2 <- new('CellType', name='B', minCycle=16-2,
    cycleLength=function() max(16-2, rnorm(1,16,1)))
type3 <- new('CellType', name='C', minCycle=22-2,
    cycleLength=function() max(22-2, rnorm(1,22,1)))
type4 <- new('CellType', name='D', minCycle=28-2,
    cycleLength=function() max(28-2, rnorm(1,28,1)))
type5 <- new('CellType', name='E', minCycle=34-2,
    cycleLength=function() max(34-2, rnorm(1,34,1)))

numReplicates <- 300
totalRuns <- 3 * numReplicates

if (arrayNum <= numReplicates)
{
    cellTypes <- c(type1, type5)
    cellTypeInitFreq <- c(rmultinom(1, 100, rep(1/2, 2))) / 100

} else if (arrayNum <= 2 * numReplicates) {

    cellTypes <- c(type1, type3, type5)
    cellTypeInitFreq <- c(rmultinom(1, 100, rep(1/3, 3))) / 100

} else {

    cellTypes <- c(type1, type2, type3, type4, type5)
    cellTypeInitFreq <- c(rmultinom(1, 100, rep(1/5, 5))) / 100

}

#### Run Simulation ####

if (!is.na(returnSize)) {

    cat(as.numeric(totalRuns))

} else {

    output <- inSilicoCellModel(initialNum=initialNum,
        runTime=runTime,
        density=density,
        boundary=boundary,
        syncCycles=syncCycles,
        randSeed=randSeed,
        modelType=modelType,
        outputIncrement=outputIncrement,
        recordIncrement=recordIncrement,
        timeIncrement=timeIncrement,
        cellTypes=cellTypes,
        cellTypeInitFreq=cellTypeInitFreq,
        drugs=drugs,
        maxDeformation=maxDeformation,
        maxTranslation=maxTranslation,
        maxRotation=maxRotation,
        nG=nG,
        epsilon=epsilon,
        delta=delta
    )

    save(output, file=paste("output_", arrayNum, ".RData", sep=""))

}



