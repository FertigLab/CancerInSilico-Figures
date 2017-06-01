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

initialNum <- 1000
runTime <- 72
allDensities <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
dim <- c(length(allDensities))
density <- allDensities[arrayNum]
recordIncrement <- 4

#### Run Simulation ####

if (!is.na(returnSize)) {

    cat(as.numeric(prod(dim)))

} else {

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

    save(output, file=paste("output_", arrayNum, ".RData", sep=""))
}
