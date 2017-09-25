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

sd <- c(0.1, 0.5, 1, 2, 4, 6, 8, 10)
freqA <- seq(0, 1.0, 0.1)

typeDEF <- lapply(sd, function(stddev) new('CellType', name='DEF', minCycle=24-2*stddev,
    cycleLength=function() max(24-2*stddev, rnorm(1,24,stddev))))
typeA <- lapply(c(10,14,18,22), function(mean) new('CellType', name='A', minCycle=mean-2,
    cycleLength=function() max(mean-2, rnorm(1,mean,1))))
typeB <- lapply(c(38,34,30,26), function(mean) new('CellType', name='B', minCycle=mean-2,
	cycleLength=function() max(mean-2, rnorm(1,mean,1))))

total <- length(sd) + length(typeA) * length(freqA)

if (arrayNum <= length(sd))
{
    cellTypeA <- typeDEF[[arrayNum]]
    cellTypeB <- typeDEF[[arrayNum]]
    cellTypeInitFreq <- c(1,0)

} else {

    n <- arrayNum - length(sd) - 1
    i1 <- (n %% length(typeA)) + 1
    i2 <- floor(n / length(typeA)) + 1

    cellTypeA <- typeA[[i1]]
    cellTypeB <- typeB[[i1]]
    cellTypeInitFreq <- c(freqA[i2], 1 - freqA[i2])
}

cellTypes <- c(cellTypeA, cellTypeB)

#### Run Simulation ####

if (!is.na(returnSize)) {

    cat(as.numeric(total))

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



