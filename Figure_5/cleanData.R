library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='~/data/figure_data/Figure_5', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

getMean <- function(cellType)
{
    round(mean(sapply(1:1e4, function(dummy) cellType@cycleLength())))
}

getSD <- function(cellType)
{
    (getMean(cellType) - cellType@minCycle) / 2
}

fileNo <- 1
pb <- txtProgressBar(min=1, max=length(allFiles), style=3)
fig5Data <- list()
for (file in allFiles)
{
    load(file)
    nCells <- sapply(0:output@runTime, getNumberOfCells, model=output)

    fig5Data[[file]] <- list(
        'numCells'    = nCells,
        'meanA' = getMean(output@cellTypes[[1]]),
        'meanB' = getMean(output@cellTypes[[2]]),
        'sdA' = getSD(output@cellTypes[[1]]),
        'sdB' = getSD(output@cellTypes[[2]]),
        'freqA' = output@cellTypeInitFreq[1]
    )
    
    fileNo <- fileNo + 1
    setTxtProgressBar(pb, fileNo)
}
close(pb)
print('saving...')
save(fig5Data, file='Figure_5_cleaned.RData')


