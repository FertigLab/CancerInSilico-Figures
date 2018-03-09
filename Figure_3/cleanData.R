library(CancerInSilico)

args <- commandArgs(TRUE)
dir <- args[1]

## read data, extract neccesary info
allFiles <- list.files(path=dir, full.names=TRUE, recursive=TRUE,
    pattern="*.RData")

getDrugEffect <- function(d)
{
    effect <- sapply(1:10000, function(dummy) d@cycleLengthEffect(0,1000))
    corrected <- 5 * round(mean(effect) / 5)
    return(corrected / 1000)
}

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
        'drugEffect'  = getDrugEffect(output@drugs[[1]]),
        'cycleLength' = output@cellTypes[[1]]@minCycle)

    fileNo <- fileNo + 1
    setTxtProgressBar(pb, fileNo)
}
close(pb)
print('saving...')
save(fig3Data, file='Figure_3_cleaned.RData')


