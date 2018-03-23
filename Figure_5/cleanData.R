library(CancerInSilico)

args <- commandArgs(TRUE)
dir <- args[1]

## read data, extract neccesary info
allFiles <- list.files(path=dir, full.names = TRUE, recursive = TRUE, pattern = "*.RData")

getMean <- function(cellType)
{
    round(mean(sapply(1:1e4, function(dummy) cellType@cycleLength())))
}

getSD <- function(cellType)
{
    (getMean(cellType) - cellType@minCycle) / 2
}

getTypeAProportion <- function(time, mod)
{
    N <- getNumberOfCells(mod, time)
    sum(sapply(1:N, function(i) getCellType(mod, time, i)==1)) / N
}

fileNo <- 1
pb <- txtProgressBar(min=1, max=length(allFiles), style=3)
fig5Data <- list()
for (file in allFiles)
{
    load(file)
   
    fig5Data[[file]] <- list(
        'numCells'    = sapply(0:output@runTime, getNumberOfCells, model=output),
        'numTypes'    = length(output@cellTypes),
	'typeName'    = output@cellTypes[[1]]@name,
        'typeFreq' = output@cellTypeInitFreq,
	'typeAProp' = getTypeAProportion(output@runTime, output)
    )
    
    fileNo <- fileNo + 1
    setTxtProgressBar(pb, fileNo)
}
close(pb)
print('saving...')
save(fig5Data, file='Figure_5_cleaned.RData')


