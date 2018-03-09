library(CancerInSilico)
library(ggplot2)
set.seed(123)

# simulate data

typeA <- new('CellType', name='A', minCycle=24-4,
    cycleLength=function() max(24-4, rnorm(1,24,1)))

typeB <- new('CellType', name='B', minCycle=36-4,
    cycleLength=function() max(36-4, rnorm(1,36,1)))

twoTypesModel <- inSilicoCellModel(initialNum=100, runTime=168, density=0.05,
    boundary=1, syncCycle=FALSE, randSeed=123, outputIncrement=4,
    recordIncrement=0.25, timeIncrement=0.001, cellTypes=c(typeA, typeB),
    cellTypeInitFreq=c(0.5,0.5), maxDeformation=0.1,
    maxTranslation=0.1, maxRotation=0.3, nG=28, epsilon=10, delta=0.2)

# generate gene expression data

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
pwyContactInhibition <- calibratePathway(pwyContactInhibition, refCountData)
allPwys <- c(pwyCellTypeA, pwyCellTypeB, pwyMitosis, pwySPhase,
    pwyContactInhibition)

ge <- inSilicoGeneExpression(twoTypesModel, allPwys, params)

cellType <- c()
nCells <- getNumberOfCells(twoTypesModel, 168)
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
    cellType[c] <- getCellType(twoTypesModel, 168, c)
}

for (t in 0:168)
{
    nCells <- getNumberOfCells(twoTypesModel, t)
    for (c in 1:nCells)
    {
        cellPhase[c,t+1] <- ifelse(SPhaseExp(twoTypesModel, c, t), 'S',
            getCellPhase(twoTypesModel, t, c))
    }
}

ge <- ge$expression
save(ge, cellType, cellPhase, file="fig6.RData")
