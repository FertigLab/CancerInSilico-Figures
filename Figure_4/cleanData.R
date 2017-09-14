library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='~/data/figure_data/Figure_4', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

fig4Data <- list()
for (file in allFiles)
{
    load(file)
    fig4Data[[file]] <- output
}

params <- new('GeneExpressionParams')
params@RNAseq <- FALSE
params@singleCell <- FALSE
params@nCells <- 50
params@sampleFreq <- 4
params@perError <- 0.05

hours <- seq(0,144,params@sampleFreq)
colNdx <- 1:length(hours)

data(SamplePathways)

pwyGrowth <- calibratePathway(pwyGrowth, referenceGeneExpression)
pwyMitosis <- calibratePathway(pwyMitosis, referenceGeneExpression)
pwySPhase <- calibratePathway(pwySPhase, referenceGeneExpression)
pwyContactInhibition <- calibratePathway(pwyContactInhibition,
    referenceGeneExpression)
allPwys <- c(pwyGrowth, pwyMitosis, pwySPhase, pwyContactInhibition)

ge_pbs <- inSilicoGeneExpression(fig4Data[[1]], allPwys, params)
ge_10ug <- inSilicoGeneExpression(fig4Data[[2]], allPwys, params)
ge_100ug <- inSilicoGeneExpression(fig4Data[[3]], allPwys, params)

ge_pbs$expression <- ge_pbs$expression[,colNdx]
ge_10ug$expression <- ge_10ug$expression[,colNdx]
ge_100ug$expression <- ge_100ug$expression[,colNdx]

movAvg <- function(data)
{
    window <- 24 / params@sampleFreq
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
    activationGrowth_pbs    =   ge_pbs$pathways[[1]][colNdx],
    mitosis_pbs             =   movAvg(ge_pbs$pathways[[2]][colNdx]),
    GtoS_pbs                =   movAvg(ge_pbs$pathways[[3]][colNdx]),
    contactInhibition_pbs   =   ge_pbs$pathways[[4]][colNdx],

    activationGrowth_10ug   =   ge_10ug$pathways[[1]][colNdx],
    mitosis_10ug            =   movAvg(ge_10ug$pathways[[2]][colNdx]),
    GtoS_10ug               =   movAvg(ge_10ug$pathways[[3]][colNdx]),
    contactInhibition_10ug  =   ge_10ug$pathways[[4]][colNdx],

    activationGrowth_100ug  =   ge_100ug$pathways[[1]][colNdx],
    mitosis_100ug           =   movAvg(ge_100ug$pathways[[2]][colNdx]),
    GtoS_100ug              =   movAvg(ge_100ug$pathways[[3]][colNdx]),
    contactInhibition_100ug =   ge_100ug$pathways[[4]][colNdx]
)

save(ge_pbs, ge_10ug, ge_100ug, pwyActivity, file='Figure_4_cleaned.RData')

