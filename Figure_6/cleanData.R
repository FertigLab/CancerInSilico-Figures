library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='~/data/figure_data/Figure_4', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

fig6Data <- list()
for (file in allFiles)
{
    load(file)
    fig6Data[[file]] <- output
}

params <- new('GeneExpressionParams')
params@RNAseq <- TRUE
params@singleCell <- TRUE
params@nCells <- 50
params@sampleFreq <- 4
params@splatParams <- splatter::setParam(params@splatParams, "dropout.present", TRUE)

hours <- seq(0,144,params@sampleFreq)
colNdx <- 1:length(hours)

data(SamplePathways)

refCountData <- round(2 ^ referenceGeneExpression - 1)

pwyGrowth <- calibratePathway(pwyGrowth, refCountData)
pwyMitosis <- calibratePathway(pwyMitosis, refCountData)
pwySPhase <- calibratePathway(pwySPhase, refCountData)
pwyContactInhibition <- calibratePathway(pwyContactInhibition, refCountData)
allPwys <- c(pwyGrowth, pwyMitosis, pwySPhase, pwyContactInhibition)

ge_pbs <- inSilicoGeneExpression(fig6Data[[1]], allPwys, params)
ge_10ug <- inSilicoGeneExpression(fig6Data[[2]], allPwys, params)
ge_100ug <- inSilicoGeneExpression(fig6Data[[3]], allPwys, params)

params@singleCell <- FALSE

ge_pbs_bulk <- inSilicoGeneExpression(fig6Data[[1]], allPwys, params)
ge_10ug_bulk <- inSilicoGeneExpression(fig6Data[[2]], allPwys, params)
ge_100ug_bulk <- inSilicoGeneExpression(fig6Data[[3]], allPwys, params)

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

save(ge_pbs, ge_10ug, ge_100ug, ge_pbs_bulk,
    ge_10ug_bulk, ge_100ug_bulk, pwyActivity, file='Figure_6_cleaned.RData')

