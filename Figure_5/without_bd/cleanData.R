library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='~/CancerInSilico-Figures/Data/Figure_5_no_bd/', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

fileNo <- 1
pb <- txtProgressBar(min=1, max=length(allFiles), style=3)
final_proportion_data <- list()
for (file in allFiles)
{
    load(file)

    init_proportion_typeA <- output@cellTypeInitFreq[[1]]
    cycle_length_typeB <- output@cellTypes[[2]]@minCycle

    finalTime <- output@runTime
    finalCellType_list <- getCellTypes(output, finalTime)

    celltype_counts <- c(0,0)
    celltype_counts[1] <- sum(finalCellType_list == 0)
    celltype_counts[2] <- sum(finalCellType_list == 1)
    finalProportions <- celltype_counts[1]/sum(celltype_counts)

    final_proportion_data[[file]] <- c(init_proportion_typeA, cycle_length_typeB, finalProportions)

    fileNo <- fileNo + 1
    setTxtProgressBar(pb, fileNo)
}
close(pb)
print('saving...')
save(final_proportion_data, file='Figure_5b_cleaned.RData')


