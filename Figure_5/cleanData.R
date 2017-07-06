library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='../Data/Figure_5', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

fileNo <- 1
pb <- txtProgressBar(min=1, max=length(allFiles), style=3)
final_proportion_data <- list()
for (file in allFiles)
{
    load(file)

    cellTypeBFreq <- output@cellTypeInitFreq[2]
    cellTypeBCycleLength <- output@cellTypes[[2]]@minCycle

    finalCellTypes <- sapply(1:getNumberOfCells(output, output@runtime),
        getCellType, model=output, time=output@runtime)

    cellTypeBFinalFreq <- sum(finalCellTypes==2) / length(finalCellTypes)

    fig5data[[file]] <- c(cellTypeBFreq, cellTypeBCycleLength, cellTypeBFinalFreq)

    fileNo <- fileNo + 1
    setTxtProgressBar(pb, fileNo)
}
close(pb)
print('saving...')
save(fig5data, file='Figure_5_cleaned.RData')


