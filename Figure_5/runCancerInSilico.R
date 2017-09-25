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

# first sweep, single cell type, fixed mean change variance
mean <- 24
sd <- c(0.1, 0.5, 1, 2, 4, 6, 8, 10)

# second sweep, two cell types, fix variance, change mean and freq
freqA <- seq(0, 1.0, 0.1)
freqB <- 1 - freqA
meanA <- c(10, 14, 18, 22)
meanB <- c(38, 34, 30, 26)
sdA <- 1
sdB <- 1

# | Afreq | Bfreq | Amean | Bmean | Asd | Bsd |
scenario1 <- matrix(nrow = length(sd), ncol = 6)
scenario1[,1] = 1
scenario1[,2] = 0
scenario1[,3] = rep(mean, nrow(scenario1))
scenario1[,4] = rep(mean, nrow(scenario1))
scenario1[,5] = sd
scenario1[,6] = rep(0, nrow(scenario1))

# | Afreq | Bfreq | Amean | Bmean | Asd | Bsd |
scenario2 <- matrix(nrow = length(meanA) * length(freqA), ncol = 6)
scenario2[,1] = rep(freqA, each=length(meanA))
scenario2[,2] = 1 - scenario2[,1]
scenario2[,3] = rep(meanA, times=length(freqA))
scenario2[,4] = rep(meanB, times=length(freqA))
scenario2[,5] = rep(sdA, times=nrow(scenario2))
scenario2[,6] = rep(sdB, times=nrow(scenario2))

allScenarios <- rbind(scenario1, scenario2) 

# parameters used
meanA <- allScenarios[arrayNum, 3]
meanB <- allScenarios[arrayNum, 4]
sdA <- allScenarios[arrayNum, 5]
sdB <- allScenarios[arrayNum, 6]

cellTypeA <- new('CellType', name='A', minCycle=meanA - 2*sdA,
    cycleLength=function() max(meanA - 2*sdA, rnorm(1,meanA,sdA)))
cellTypeB <- new('CellType', name='B', minCycle=meanB - 2*sdB,
    cycleLength=function() max(meanB - 2*sdB, rnorm(1,meanB,sdB)))

cellTypes <- c(cellTypeA, cellTypeB)
cellTypeInitFreq <- c(allScenarios[arrayNum, 1], allScenarios[arrayNum, 2])

#### Run Simulation ####

if (!is.na(returnSize)) {

    cat(as.numeric(nrow(allScenarios)))

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



