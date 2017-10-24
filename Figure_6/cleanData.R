library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='~/data/figure_data/Figure_5', full.names = TRUE,
    recursive = TRUE, pattern = "*.RData")

fig6Data <- list()
for (file in allFiles)
{
    load(file)
    if (output@cellTypeInitFreq[1] > 0.35 & output@cellTypeInitFreq[1] < 0.45)
    {
        fig6Data[[1]] <- output
        break
    }
}

params <- new('GeneExpressionParams')
params@RNAseq <- TRUE
params@singleCell <- TRUE
params@nCells <- 96
params@sampleFreq <- 4
params@splatParams <- splatter::setParam(params@splatParams, "dropout.present", TRUE)

hours <- seq(0,144,params@sampleFreq)
colNdx <- 1:length(hours)

data(SamplePathways)

refCountData <- round(2 ^ referenceGeneExpression - 1)

pwyCellTypeA <- new('Pathway', genes = pwyGrowth@genes[1:241],
    transformMidpoint = 0.05, transformSlope = 5 / 0.1,
    expressionScale = function(model, cell, time)
    {
        type <- getCellType(model, time, cell)
        if (type == 1)
            return(1)
        else
            return(0)
    }
)

pwyCellTypeB <- new('Pathway', genes = pwyGrowth@genes[242:482],
    transformMidpoint = 0.05, transformSlope = 5 / 0.1,
    expressionScale = function(model, cell, time)
    {
        type <- getCellType(model, time, cell)
        if (type == 2)
            return(1)
        else
            return(0)
    }
)

pwyCellTypeA <- calibratePathway(pwyCellTypeA, refCountData)
pwyCellTypeB <- calibratePathway(pwyCellTypeB, refCountData)
pwyMitosis <- calibratePathway(pwyMitosis, refCountData)
pwySPhase <- calibratePathway(pwySPhase, refCountData)
allPwys <- c(pwyCellTypeA, pwyCellTypeB, pwyMitosis, pwySPhase)

ge <- inSilicoGeneExpression(fig6Data[[1]], allPwys, params)

params@singleCell <- FALSE

ge_bulk <- inSilicoGeneExpression(fig6Data[[1]], allPwys, params)

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
    cellTypeA = ge_bulk$pathways[[1]][colNdx],
    cellTypeB = ge_bulk$pathways[[2]][colNdx],
    GtoM      = movAvg(ge_bulk$pathways[[3]][colNdx]),
    GtoS      = movAvg(ge_bulk$pathways[[4]][colNdx])
)

# get cell phase/type info
cellType <- c()
nCells <- getNumberOfCells(fig6Data[[1]], 144)
cellPhase <- matrix(nrow=nCells, ncol=145)

SPhaseExp <- function(model, cell, time)
{
    window <- c(max(time - 1, 0), min(time + 1, model@runTime))

    r1 <- getRadius(model, window[1], cell)
    r2 <- getRadius(model, window[2], cell)

    type <- model@cellTypes[[getCellType(model, time, cell)]]

    return(ifelse(r1 < sqrt(1.5 * type@size) & r2 > sqrt(1.5 * type@size),
        1, 0))
}

for (c in 1:nCells)
{
    cellType[c] <- getCellType(fig6Data[[1]], 144, c)
    for (t in 0:144)
    {
        if (SPhaseExp(fig6Data[[1]], c, t)
            cellPhase[c,t+1] <- 'S'
        else
            cellPhase[c,t+1] <- getCellPhase(fig6Data[[1]], t, c)
    }
}

save(pwyActivity, ge, ge_bulk, cellPhase, cellType, file='Figure_6_cleaned.RData')

