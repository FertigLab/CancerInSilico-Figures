library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='../Data/Figure_2', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

fileNo <- 1
pb <- txtProgressBar(min=1, max=length(allFiles), style=3)
fig2Data <- list()
for (file in allFiles)
{
    load(file)
    nCells <- sapply(0:output@runTime, getNumberOfCells, model=output)
    contact <- sapply(0:output@runTime, function(t)
        mean(getContactInhibition(output, t)))
    den <- sapply(0:output@runTime, getDensity, model=output)

    fig2Data[[file]] <- list('boundary'=output@boundary,
        'initDensity'=output@density, 'numCells'=nCells,
        'cycleLength'=output@cellTypes[[1]]@minCycle,
        'contactInhibition'=contact, 'density'=den)

    fileNo <- fileNo + 1
    setTxtProgressBar(pb, fileNo)
}
close(pb)
print('saving...')
save(fig2Data, file='Figure_2_cleaned.RData')


