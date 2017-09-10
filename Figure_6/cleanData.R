library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='~/data/figure_data/Figure_4', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

print('saving...')
save(fig6data, file='Figure_6_cleaned.RData')


