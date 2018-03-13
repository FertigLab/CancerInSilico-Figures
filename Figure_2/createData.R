# load libraries

library(CancerInSilico)
library(foreach)
library(doParallel)

# command line arguments

args <- commandArgs(TRUE)
nCores <- ifelse(length(args) == 0, 1, as.integer(args[1]))
print(paste("running on", nCores, "cores"))

# simulate data

runModel <- function(cycleLength, density, boundary)
{
    type <- new('CellType', name='DEFAULT', minCycle=cycleLength,
        cycleLength=function() cycleLength)
    output <- inSilicoCellModel(initialNum=100, runTime=168, density=density,
        boundary=boundary, syncCycle=FALSE, randSeed=123, outputIncrement=4,
        recordIncrement=0.25, timeIncrement=0.001, cellTypes=c(type),
        cellTypeInitFreq=c(1), maxDeformation=0.1, maxTranslation=0.1,
        maxRotation=0.3, nG=28, epsilon=10, delta=0.2)

    time <- 0:output@runTime
    N <- length(time)
    den <- sapply(time, getDensity, model=output)
    nCells <- sapply(time, getNumberOfCells, model=output)
    cycleLength <- rep(cycleLength, N)
    density <- rep(density, N)
    boundary <- rep(boundary, N)
    
    cbind(time, density, cycleLength, boundary, nCells, den)    
}

cycleLength <- c(12,24,48)
density <- c(0.1, 0.2, 0.3)
boundary <- c(0,1)

dim <- c(length(cycleLength), length(density), length(boundary))
indexArray <- array(1:prod(dim), dim)

cl <- makeCluster(nCores)
registerDoParallel(cl)
data <- foreach(i = 1:prod(dim), .packages="CancerInSilico") %dopar%
{
    ndx <- which(indexArray==i, arr.ind=TRUE)
    runModel(cycleLength[ndx[1]], density[ndx[2]], boundary[ndx[3]])
}
stopCluster(cl)

mat <- do.call(rbind, data)
df <- data.frame(time=mat[,1], initialDensity=mat[,2], cycleLength=mat[,3],
    boundary=ifelse(mat[,4] > 0, TRUE, FALSE), nCells=mat[,5], density=mat[,6])
save(df, file="fig2.RData")
