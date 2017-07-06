library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='../Data/Figure_4', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

fig4Data <- list()
for (file in allFiles)
{
    load(file)
    fig4Data[[file]] <- output
}

data(SamplePathways)

sampFreq <- 4
nCells <- 50
hours <- seq(0,144,sampFreq)
col_ndx <- 1:length(hours)

pwyGrowth <- calibratePathway(pwyGrowth, referenceGeneExpression)
pwyMitosis <- calibratePathway(pwyMitosis, referenceGeneExpression)
pwySPhase <- calibratePathway(pwySPhase, referenceGeneExpression)
pwyContactInhibition <- calibratePathway(pwyContactInhibition,
    referenceGeneExpression)
allPwys <- c(pwyGrowth, pwyMitosis, pwySPhase, pwyContactInhibition)

ge_pbs <- inSilicoGeneExpression(fig4Data[[1]], allPwys, nCells=nCells,
    sampFreq=sampFreq, perError=0.05)

ge_10ug <- inSilicoGeneExpression(fig4Data[[2]], allPwys, nCells=nCells,
    sampFreq=sampFreq, perError=0.05)

ge_100ug <- inSilicoGeneExpression(fig4Data[[3]], allPwys, nCells=nCells,
    sampFreq=sampFreq, perError=0.05)

ge_pbs$expression <- ge_pbs$expression[,col_ndx]
ge_10ug$expression <- ge_10ug$expression[,col_ndx]
ge_100ug$expression <- ge_100ug$expression[,col_ndx]

save(ge_pbs, ge_10ug, ge_100ug, file='Figure_4_cleaned.RData')

