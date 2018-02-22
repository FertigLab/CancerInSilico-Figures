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
   
    fig5Data[[file]] <- list(
        'density'    = sapply(0:output@runTime, getDensity, model=output),
        'numTypes'    = length(output@cellTypes),
	'typeName'    = output@cellTypes[[1]]@name,
        'typeFreq' = output@cellTypeInitFreq
    )
    
    fileNo <- fileNo + 1
    setTxtProgressBar(pb, fileNo)
}
close(pb)
print('saving...')
save(fig5Data, file='Figure_5_cleaned.RData')


