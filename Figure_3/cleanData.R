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

getCellArea <- function(model, time)
{
    N <- getNumberOfCells(model, time)
    radius <- sapply(1:N, getRadius, model=model, time=time)
    axis_length <- sapply(1:N, getAxisLength, model=model, time=time)
    interphase <- axis_length <= 2 * radius
    N_mitosis <- length(!interphase)
    return(sum(pi * radius[interphase]^2) + 2 * pi * N_mitosis)
}

fileNo <- 1
pb <- txtProgressBar(min=1, max=length(allFiles), style=3)
fig3Data <- list()
for (file in allFiles)
{
    load(file)
    cellArea <- sapply(0:output@runTime, getCellArea, model=output)

    fig3Data[[file]] <- list(
        'initDensity' = output@density,
        'cellArea'    = cellArea,
        'drugEffect'  = getDrugEffect(output@drugs[[1]]),
        'cycleLength' = output@cellTypes[[1]]@minCycle)

    fileNo <- fileNo + 1
    setTxtProgressBar(pb, fileNo)
}
close(pb)
print('saving...')
save(fig3Data, file='Figure_3_cleaned.RData')


