library('CancerInSilico')
library(methods)

args <- commandArgs(TRUE)
arrayNum <- as.integer(args[1])
returnSize <- as.integer(args[2])

#### Set Defaults ####

initialNum <- 100
runTime <- 168
density <- 0.1
boundary <- 1
syncCycles <- FALSE
randSeed <- 0
outputIncrement <- 4
recordIncrement <- 1
timeIncrement <- 0.001
cellTypes <- c(new('CellType', name='DEFAULT', minCycle=40, cycleLength=function() 40))
cellTypeInitFreq <- c(1)
drugs <- list()
maxDeformation <- 0.1
maxTranslation <- 0.1
maxRotation <- 0.3
nG <- 28
epsilon <- 10.0
delta <- 0.2

#### Set Custom Values ####

allDrugs <- lapply(seq(1.0, 2.0, 0.5), function(l) new('Drug',
    name='DEFAULT', timeAdded=24, cycleLengthEffect=function(a,b) b*l))

dim <- c(length(allDrugs))
indexArray <- array(1:prod(dim), dim)
index <- which(indexArray==arrayNum, arr.ind=TRUE)

drugs <- c(allDrugs[index[1]])

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
