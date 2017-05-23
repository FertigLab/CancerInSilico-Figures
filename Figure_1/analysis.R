## read data
allFiles <- list.files(path='../Data/Figure_1', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

modelOutput <- list()
for (file in allFiles)
{
    load(file)
    totalCells <- sapply(0:output@runTime, getNumberOfCells, model=output)
    modelOutput[[file]] <- c(output@boundary, output@density,
        output@cellTypes[[1]]@minCycle, totalCells)
}


