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

movAvg <- function(data)
{
    window <- 24 / sampFreq
    avg <- c()    
    for (i in 1:length(data))
    {
        mn <- max(1, i - window)
        mx <- min(length(data), i + window)
        avg <- c(avg, mean(data[mn:mx]))
    }
    return(avg)
}

pwyActivity <- data.frame(hour=hours,
    activationGrowth_pbs    =   ge_pbs$pathways[[1]][col_ndx],
    mitosis_pbs             =   movAvg(ge_pbs$pathways[[2]][col_ndx]),
    GtoS_pbs                =   movAvg(ge_pbs$pathways[[3]][col_ndx]),
    contactInhibition_pbs   =   ge_pbs$pathways[[4]][col_ndx],

    activationGrowth_10ug   =   ge_10ug$pathways[[1]][col_ndx],
    mitosis_10ug            =   movAvg(ge_10ug$pathways[[2]][col_ndx]),
    GtoS_10ug               =   movAvg(ge_10ug$pathways[[3]][col_ndx]),
    contactInhibition_10ug  =   ge_10ug$pathways[[4]][col_ndx],

    activationGrowth_100ug  =   ge_100ug$pathways[[1]][col_ndx],
    mitosis_100ug           =   movAvg(ge_100ug$pathways[[2]][col_ndx]),
    GtoS_100ug              =   movAvg(ge_100ug$pathways[[3]][col_ndx]),
    contactInhibition_100ug =   ge_100ug$pathways[[4]][col_ndx]
)

save(hours, ge_pbs, ge_10ug, ge_100ug, pwyActivity,
    file='Figure_4_cleaned.RData')

