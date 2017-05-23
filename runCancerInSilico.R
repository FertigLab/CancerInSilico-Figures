library('CancerInSilico')
library(methods)

args <- commandArgs(TRUE)
arrayNum <- as.integer(args[1])
jobName <- args[2]

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

#allCellTypes <- lapply(seq(12,48,4), function(l) new('CellType',
#    name='B', minCycle=l, cycleLength=function() l))

#ctA <- new('CellType', name='A', cycleLength=function() 24, minCycle=24)
#ctB <- allCellTypes[arrayNum]
#cellTypes <- c(ctA, ctB)
#cellTypeInitFreq <- c(0.55, 0.45)

#### Run Simulation ####

repeat
{
    output <- runCellSimulation(initialNum=initialNum,
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
    randSeed <- randSeed + 100

    if (length(output@cells) > output@runTime / output@recordIncrement)
        break
}

save(output, file=paste("output_", jobName, "_", arrayNum, ".RData", sep=""))
