library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='../Data/Figure_3', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

fileNo <- 1
pb <- txtProgressBar(min=1, max=length(allFiles), style=3)
fig3Data <- list()
for (file in allFiles)
{
    load(file)
    nCells <- sapply(0:output@runTime, getNumberOfCells, model=output)

    fig3Data[[file]] <- list(
        'initDensity' = output@density,
        'numCells'    = nCells,
        'drugEffect'  = output@drugs[[1]]@cycleLengthEffect(0,1),
        'cycleLength' = output@cellTypes[[1]]@minCycle)

    fileNo <- fileNo + 1
    setTxtProgressBar(pb, fileNo)
}
close(pb)
print('saving...')
save(fig3Data, file='Figure_3_cleaned.RData')


