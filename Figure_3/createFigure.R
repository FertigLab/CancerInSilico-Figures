############################# load libraries ###################################

library(CancerInSilico)
library(Rtsne)
library(viridis)
library(rgl)
set.seed(42)

########################## adjustable parameters ###############################

totalTime <- 168
initNumCells <- 100
initDensity <- 0.05
nGenes <- 100

######################## define multiple cell types ############################

cellTypeA <- new("CellType", name="type.A", size=1, minCycle=20,
    cycleLength=function() 20 + rexp(1, 1/4))

cellTypeB <- new("CellType", name="type.B", size=1, minCycle=32,
    cycleLength=function() 32 + rexp(1, 1/4))

allCellTypes <- list(cellTypeA, cellTypeB)

########################### simulate cell growth ###############################

mod <- inSilicoCellModel(initialNum=initNumCells, runTime=totalTime,
    density=initDensity, outputIncrement=12, cellTypes=allCellTypes,
    cellTypeInitFreq=c(0.5,0.5), syncCycles=FALSE, randSeed=42)

######################### get cell types by index ##############################

types <- sapply(allCellTypes, function(t) t@name)
nCells <- getNumberOfCells(mod, totalTime)
cellType <- c()
for (c in 1:nCells)
{
    cellType[c] <- types[getCellType(mod, totalTime, c)]
}

extractCellIndex <- function(str)
{
    as.numeric(substring(strsplit(str, "_")[[1]][1], 2))
}

extractTimePoint <- function(str)
{
    as.numeric(gsub("t", "", strsplit(str, "_")[[1]][2]))
}

convertNamesToCellTypes <- function(cnames)
{
    ndx <- unname(sapply(cnames, extractCellIndex))
    return(cellType[ndx])
}

convertNamesToTimes <- function(cnames)
{
    unname(sapply(cnames, function(x) extractTimePoint(x)))
}

sphaseActivityFunction <- function(model, cell, time)
{
    r1 <- getRadius(model, max(time - 1, 0), cell)
    r2 <- getRadius(model, time, cell)
    if (is.na(r1))
        return(0)
    else
        return(as.numeric(r1 < sqrt(1.5) & r2 > sqrt(1.5)))
}

convertNamesToPhases <- function(model, cnames)
{
    unname(sapply(cnames, function(str)
    {
        ndx <- extractCellIndex(str)
        time <- extractTimePoint(str)
        act <- sphaseActivityFunction(model, ndx, time)
        return(ifelse(act == 1, "S", getCellPhase(model, time, ndx)))
    }))
}

############################### define pathways ################################

typeAPwy <- new("Pathway", genes=paste("typeA.", 1:nGenes, sep=""),
    expressionScale=function(model, cell, time)
    {
        as.numeric(getCellType(model, time, cell) == 1)
    }
)

typeBPwy <- new("Pathway", genes=paste("typeB.", 1:nGenes, sep=""),
    expressionScale=function(model, cell, time)
    {
        as.numeric(getCellType(model, time, cell) == 2)
    }
)

sphasePwy <- new("Pathway", genes=paste("sphase.", 1:nGenes, sep=""),
    expressionScale=sphaseActivityFunction)

mitosisPwy <- new("Pathway", genes=paste("mitosis.", 1:nGenes, sep=""),
    expressionScale=function(model, cell, time)
    {
        a1 <- getAxisLength(model, max(time - 1, 0), cell)
        a2 <- getAxisLength(model, time, cell)
        if (is.na(a1))
            return(0)
        else
            return(as.numeric(a2 < a1))
    }
)

contactInhibitionPwy <- new("Pathway", genes=paste("ci.", 1:nGenes, sep=""),
    expressionScale=function(model, cell, time)
    {
        getLocalDensity(model, time, cell, 3.3)
    }
)

######################## set gene expression parameters ########################

params <- new("GeneExpressionParams")
params@randSeed <- 42
params@nCells <- 100
params@sampleFreq <- 4
params@RNAseq <- TRUE
params@singleCell <- TRUE
params@dropoutPresent <- TRUE
params@bcvCommon <- 0.001
params@combineFUN <- max

calLambda <- 10
calStdDev <- 2

#################################### figure ####################################

allPwys <- c(typeAPwy, typeBPwy, sphasePwy, mitosisPwy, contactInhibitionPwy)
allPwys <- calibratePathways(allPwys, lambda=calLambda, stddev=calStdDev)
res <- inSilicoGeneExpression(mod, allPwys, params)

tsne_out <- Rtsne(t(res$expression), dims=3, initial_dims=3, perplexity=30,
    theta=0.5, pca=TRUE, max_iter=1000, verbose=TRUE)

# part A, color by cell type

sampledTypes <- convertNamesToCellTypes(colnames(res$expression))
typeColor <- c("red", "blue")[1 + as.numeric(sampledTypes=="type.A")]

plot3d(tsne_out$Y, col=typeColor, xlab="", ylab="", zlab="")
legend3d("topright", legend=c("Type A", "Type B"), pch=16, cex=1,
    col=c("red", "purple"), inset=c(0.04))
rgl.postscript("fig3A_CELL_TYPE.pdf", fmt="pdf")

# part B, color by time

myColorRamp <- function(values, palette=viridis(255))
{
    if (min(values) < 0)
    {  
       values <- values + abs(min(values))
    }
    v <- (values - min(values)) / diff(range(values))
    x <- colorRamp(palette)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
}

sampledTimes <- convertNamesToTimes(colnames(res$expression))
timeColor <- myColorRamp(sampledTimes)

plot3d(tsne_out$Y, col=timeColor, xlab="", ylab="", zlab="")
legend3d("topright", legend=c("0 hours", "168 hours"), pch=16, cex=1,
    col=c(min(timeColor), max(timeColor)), inset=c(0.02))
rgl.postscript("fig3B_TIME.pdf", fmt="pdf")

# part C, color by phase

sampledPhases <- convertNamesToPhases(mod, colnames(res$expression))
phaseColor <- rep("blue", length(sampledPhases))
phaseColor[sampledPhases == "S"] <- "orange"
phaseColor[sampledPhases == "M"] <- "red"

plot3d(tsne_out$Y, col=phaseColor, xlab="", ylab="", zlab="")
legend3d("topright", legend=c("Interphase", "SPhase", "Mitosis"), pch=16, cex=1,
    col=c("blue", "orange", "red"), inset=c(0.02))
rgl.postscript("fig3C_CELL_PHASE.pdf", fmt="pdf", drawText=TRUE)