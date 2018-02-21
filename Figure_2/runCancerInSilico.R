library('CancerInSilico')
library(methods)

args <- commandArgs(TRUE)
arrayNum <- as.integer(args[1])
returnSize <- as.integer(args[2])

#### Set Defaults ####

initialNum <- 80
runTime <- 168
density <- 0.2
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

allDensities <- seq(0.05,0.45,0.05)
allBoundaries <- c(0,1)
allCellTypes <- lapply(seq(12,48,4), function(l) new('CellType',
    name='DEFAULT', minCycle=l, cycleLength=function() l))

dim <- c(length(allDensities), length(allBoundaries), length(allCellTypes))
indexArray <- array(1:prod(dim), dim)
index <- which(indexArray==arrayNum, arr.ind=TRUE)

density <- allDensities[index[1]]
boundary <- allBoundaries[index[2]]
cellTypes <- c(allCellTypes[index[3]])

# 2a - cycleLength = {12,20,28,36,44}
# 2b - density = {0.05, 0.10, 0.20, 0.30, 0.40}

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
