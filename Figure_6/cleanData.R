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

pwyActivity <- data.frame(cellTypeA = ge_bulk$pathways[[1]],
    cellTypeB = ge_bulk$pathways[[2]],
    GtoM      = movAvg(ge_bulk$pathways[[3]]),
    GtoS      = movAvg(ge_bulk$pathways[[4]])
)

# get cell phase/type info
cellType <- c()
nCells <- getNumberOfCells(fig6Data[[1]], 168)
cellPhase <- matrix('NoCell', nrow=nCells, ncol=169)

SPhaseExp <- function(model, cell, time)
{
    window <- c(max(time - 1, 0), min(time + 1, model@runTime))

    r1 <- getRadius(model, window[1], cell)
    r2 <- getRadius(model, window[2], cell)

    if (is.na(r1)) return(FALSE)
    return(r1 < sqrt(1.5) & r2 > sqrt(1.5))
}

for (c in 1:nCells)
{
    cellType[c] <- getCellType(fig6Data[[1]], 168, c)
}

for (t in 0:168)
{
    nCells <- getNumberOfCells(fig6Data[[1]], t)

    for (c in 1:nCells)
    {
        if (SPhaseExp(fig6Data[[1]], c, t)) {
            cellPhase[c,t+1] <- 'S'
        } else {
            cellPhase[c,t+1] <- getCellPhase(fig6Data[[1]], t, c)
        }
    }
}

save(pwyActivity, ge, ge_bulk, cellPhase, cellType, file='Figure_6_cleaned.RData')

