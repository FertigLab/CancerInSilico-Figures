library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='../Data/Figure_5', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

## extract cell type proportion
getCellTypeBFreq <- function(time, mod)
{
    cellTypes <- sapply(1:getNumberOfCells(mod, time),
        getCellType, model=mod, time=time)
    freq <- sum(cellTypes==2) / length(cellTypes)
} 

## setup progress bar
fileNo <- 1
pb <- txtProgressBar(min=1, max=length(allFiles), style=3)

## record cell info
cellTypeBFreq <- c()
cycleLength <- c()
density <- c()
time <- c()
for (file in allFiles)
{
    load(file)

    t <- 0:output@runTime

    time <- c(time, t)
    density <- c(density, rep(output@density, length(t)))
    cycleLength <- c(cycleLength, rep(output@cellTypes[[2]]@minCycle, length(t)))
    cellTypeBFreq <- c(cellTypeBFreq, sapply(t, getCellTypeBFreq, mod=output))

    fileNo <- fileNo + 1
    setTxtProgressBar(pb, fileNo)
}

## store cell info in data frame
fig5data <- data.frame(cellTypeBFreq = cellTypeBFreq,
    cycleLength = cycleLength, density = density, time = time)

## close progress bar and save data
close(pb)
print('saving...')
save(fig5data, file='Figure_5_cleaned.RData')


