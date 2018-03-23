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

typeA <- new('CellType', name='A', minCycle=8,
    cycleLength=function() 8 + rexp(1,1/4))

typeB <- new('CellType', name='B', minCycle=32,
    cycleLength=function() 32 + rexp(1,1/4))

nReplicates <- 40 # number of runs per variance
nVars <- 10 # number of x-axis variances

allSingleCellTypes <- lapply(1:nVars, function(v) new('CellType', name=paste('DEFAULT_', v, sep=""),
	minCycle=24 - v, cycleLength=function() 24 - v + rexp(1,1/v)))

if (arrayNum <= nVars * nReplicates) { # one type
	
	cellTypes <- c(allSingleCellTypes[floor((arrayNum - 1) / nReplicates) + 1])

} else { # two types

	range <- (floor((arrayNum - (nVars * nReplicates + 1)) / nReplicates) + 1) * 0.05
	aProp <- runif(1,0.5-range,0.5+range)
	cellTypeInitFreq <- c(aProp, 1 - aProp)
	cellTypes <- c(typeA, typeB)
}

boundary <- 0
randSeed <- arrayNum

#### Run Simulation ####

if (!is.na(returnSize)) {

    cat(as.numeric(2 * nVars * nReplicates))

} else {

    output <- inSilicoCellModel(initialNum=initialNum,
        runTime=runTime,
        density=density,
        boundary=boundary,
        syncCycles=syncCycles,
        randSeed=randSeed,
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



