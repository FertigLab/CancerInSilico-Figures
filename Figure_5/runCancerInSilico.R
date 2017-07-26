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

cellTypeA <- new('CellType', name='A', cycleLength=function() 32, minCycle=32)

allCellTypeBFreq <- seq(0,0.5,0.05)
allCellTypeB <- lapply(seq(8,48,6), function(l) new('CellType',
    name='B', minCycle=l, cycleLength=function() l))
allDensities <- c(1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 0.01, 0.05)

dim <- c(length(allCellTypeBFreq), length(allCellTypeB), length(allDensities))
indexArray <- array(1:prod(dim), dim)
index <- which(indexArray==arrayNum, arr.ind=TRUE)

cellTypeBFreq <- allCellTypeBFreq[index[1]]
cellTypeB <- allCellTypeB[index[2]]
density <- allDensities[index[3]]

cellTypes <- c(cellTypeA, cellTypeB)
cellTypeInitFreq <- c(1 - cellTypeBFreq, cellTypeBFreq)

#### Run Simulation ####

if (!is.na(returnSize)) {

    cat(as.numeric(prod(dim)))

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



