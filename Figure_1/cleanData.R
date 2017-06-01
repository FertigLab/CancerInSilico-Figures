library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='../Data/Figure_1', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

fileNo <- 1
pb <- txtProgressBar(min=1, max=length(allFiles), style=3)
fig1Data <- list()
for (file in allFiles)
{
    load(file)
    fig1Data[[file]] <- output
}
close(pb)
print('saving...')
save(fig1Data, file='Figure_1_cleaned.RData')


