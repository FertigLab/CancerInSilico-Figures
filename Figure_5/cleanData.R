library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='../Data/Figure_5', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

fileNo <- 1
pb <- txtProgressBar(min=1, max=length(allFiles), style=3)
cellTypeBFreq <- c()
cycleLength <- c()
density <- c()
for (file in allFiles)
{
    load(file)

    getCellTypeBFreq <- function(time)
    {
        cellTypes <- sapply(1:getNumberOfCells(output, time),
            getCellType, model=output, time=time)
        freq <- sum(cellTypes==2) / length(cellTypes)
    } 

    N <- output@runTime
    cellTypeBFreq <- sapply(0:N, getCellTypeBFreq)

    cycleLength <- rep(output@cellTypes[[2]]@minCycle, N + 1)
    density <- rep(output@density, N + 1)

    fileNo <- fileNo + 1
    setTxtProgressBar(pb, fileNo)
}
fig5data <- data.frame(cellTypeBFreq = cellTypeBFreq,
    cycleLength = cycleLength, density = density)
close(pb)
print('saving...')
save(fig5data, file='Figure_5_cleaned.RData')


