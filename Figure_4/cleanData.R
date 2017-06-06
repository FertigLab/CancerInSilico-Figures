library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='../Data/Figure_4', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

fileNo <- 1
pb <- txtProgressBar(min=1, max=length(allFiles), style=3)
fig4Data <- list()
for (file in allFiles)
{
    load(file)
    fig4Data[[file]] <- output
    
    fileNo <- fileNo + 1
    setTxtProgressBar(pb, fileNo)
}
close(pb)
print('saving...')
save(fig4Data, file='Figure_4_cleaned.RData')


