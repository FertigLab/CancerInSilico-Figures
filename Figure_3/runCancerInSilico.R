library('CancerInSilico')
library(methods)

args <- commandArgs(TRUE)
arrayNum <- as.integer(args[1])
returnSize <- as.integer(args[2])

#### Set Defaults ####

initialNum <- 100
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

allDensities <- seq(0.001, 0.1, length=10)
allCellTypes <- lapply(seq(24,36,6), function(l) new('CellType',
    name='DEFAULT', minCycle=l, cycleLength=function() l))
allDrugs <- lapply(seq(1.0, 2.0, 0.05), function(l) new('Drug',
    name='DEFAULT', timeAdded=24, cycleLengthEffect=function(a,b)
    rnorm(n=1, mean=b*l, sd=4)))

dim <- c(length(allDensities), length(allCellTypes), length(allDrugs))
indexArray <- array(1:prod(dim), dim)
index <- which(indexArray==arrayNum, arr.ind=TRUE)

density <- allDensities[index[1]]
cellTypes <- c(allCellTypes[index[2]])
drugs <- c(allDrugs[index[3]])

#### Run Simulation ####

if (!is.na(returnSize)) {

    cat(as.numeric(prod(dim)))

} else {

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

    save(output, file=paste("output_", arrayNum, ".RData", sep=""))

}
