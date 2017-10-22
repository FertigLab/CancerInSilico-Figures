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

#Main source of variation is initial proportion
#Three parameters: initial proportion, number of cell types, difference in cell types
#A – show variance with one type is small
#B – show variance is large with two cell types and varying prop
#C – show adding cell types does not increase this variance
#D – show increasing variance with increase variance in prop
#Each fig has 3 components: total over time, prop dist, total dist at 48 hours

singleType <- new('CellType', name='DEFAULT', minCycle=24-20,
    cycleLength=function() max(24-20, rnorm(1,24,10)))

type1 <- new('CellType', name='A', minCycle=12-2,
    cycleLength=function() max(12-2, rnorm(1,12,1)))

type2 <- new('CellType', name='B', minCycle=16-2,
    cycleLength=function() max(16-2, rnorm(1,16,1)))

type3 <- new('CellType', name='C', minCycle=20-2,
    cycleLength=function() max(20-2, rnorm(1,20,1)))

type4 <- new('CellType', name='D', minCycle=28-2,
    cycleLength=function() max(28-2, rnorm(1,28,1)))

type5 <- new('CellType', name='E', minCycle=32-2,
    cycleLength=function() max(32-2, rnorm(1,32,1)))

type6 <- new('CellType', name='E', minCycle=36-2,
    cycleLength=function() max(36-2, rnorm(1,36,1)))

numReplicates <- 200
totalRuns <- 4 * numReplicates

if (arrayNum <= numReplicates)
{
    cellTypes <- c(singleType)
    cellTypeInitFreq <- c(1)

} else if (arrayNum <= 4 * numReplicates) {

    cellTypes <- c(type1, type6)
    aProp <- runif(1,0,1)
    cellTypeInitFreq <- c(aProp, 1 - aProp)
}

#} else if (arrayNum <= 5 * numReplicates) {
#
#    cellTypes <- c(type1, type2, type5, type6)
#    cellTypeInitFreq <- c(rmultinom(1, 100, rep(1/4, 4))) / 100
#
#} else {
#
#    cellTypes <- c(type1, type2, type3, type4, type5, type6)
#    cellTypeInitFreq <- c(rmultinom(1, 100, rep(1/6, 6))) / 100
#
#}

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



