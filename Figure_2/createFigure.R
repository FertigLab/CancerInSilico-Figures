############################# load libraries ###################################

library(CancerInSilico)
library(gplots)
library(viridis)
library(ggplot2)
set.seed(42)

########################## adjustable parameters ###############################

totalTime <- 168
initNumCells <- 100
initDensity <- 0.05
nGenes <- 50

######################## define multiple cell types ############################

allCellTypes <- list(new("CellType", name="type.A", size=1, minCycle=24,
    cycleLength=function() 24 + rexp(1, 1/4)))

############################## define a drug ###################################

allDrugs <- list(new("Drug", name="GrowthSuppressor", timeAdded=0, 
    cycleLengthEffect=function(a,b) rnorm(n=1, mean=2*b, sd=4)))

########################### simulate cell growth ###############################

mod <- inSilicoCellModel(initialNum=initNumCells, runTime=totalTime,
    density=initDensity, outputIncrement=12, cellTypes=allCellTypes,
    cellTypeInitFreq=c(1), syncCycles=FALSE, randSeed=42)

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

############################### define pathways ################################

mitosisPwy <- new("Pathway", genes=paste("mitosis.", 1:nGenes, sep=""),
    expressionScale=function(model, cell, time)
    {
        a1 <- getAxisLength(model, max(time - 1, 0), cell)
        a2 <- getAxisLength(model, min(time + 1, model@runTime), cell)
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

sphasePwy <- new("Pathway", genes=paste("sphase.", 1:nGenes, sep=""),
    expressionScale=function(model, cell, time)
    {
        r1 <- getRadius(model, max(time - 1, 0), cell)
        r2 <- getRadius(model, min(time + 1, model@runTime), cell)
        if (is.na(r1))
            return(0)
        else
            return(as.numeric(r1 < sqrt(1.5) & r2 > sqrt(1.5)))
    }
)

######################## set gene expression parameters ########################

calLambda <- 10
calStdDev <- 2

params <- new("GeneExpressionParams")
params@randSeed <- 42
params@nCells <- 50
params@sampleFreq <- 4
params@singleCell <- FALSE
params@bcvCommon <- 0.001
params@combineFUN <- max
params@nDummyGenes <- 0 * nGenes
params@dummyDist <- function(N) rnorm(N, mean=rexp(1, 1/calLambda), sd=calStdDev)
params@perError <- 0.2

#################################### figure ####################################

allPwys <- c(mitosisPwy, contactInhibitionPwy, sphasePwy)
allPwys <- calibratePathways(allPwys, lambda=calLambda, stddev=calStdDev)

# part A, bulk microarray visualization

params@RNAseq <- FALSE
microarray <- inSilicoGeneExpression(mod, allPwys, params)

pdf("Fig2a.pdf")
ndx <- apply(microarray$expression, 1, var) == 0 # remove zero variance rows
plotData <- microarray$expression[!ndx,]
plotData <- plotData[sample(1:nrow(plotData)),]
heatmap.2(plotData, 
    col=greenred, scale="row",
    trace="none", hclust=function(x) hclust(x,method="complete"),
    distfun=function(x) as.dist((1-cor(t(x)))/2), 
    Colv=FALSE, dendrogram="row",
    labRow = FALSE, labCol = FALSE,
    main="Bulk Gene Expression",
    xlab="time")
dev.off()

# part B, PCA colored by time

times <- seq(0, mod@runTime, params@sampleFreq)
timeColor <- myColorRamp(times)

pca <- prcomp(microarray$expression, center=FALSE, scale.=FALSE)
dfB <- data.frame(time=times, mitosis=unname(microarray$pathways[[1]]),
    sphase=unname(microarray$pathways[[3]]),
    pc1=pca$rotation[,1], pc2=pca$rotation[,2], pc3=pca$rotation[,3])

fig <- ggplot(dfB, aes(x=pc1, y=pc2, color=time, size=1.5)) + 
    geom_point() + 
    labs(title="PCA by Time", x="PC1", y="PC2") +
    theme_classic() + 
    scale_colour_gradientn(colours=timeColor)
ggsave(filename="Fig2b_TIME.pdf", plot=fig)

# part C, PCA colored by phase

phase <- rep("Interphase", length(times))
phase[dfB$mitosis > dfB$sphase] <- "Mitosis"
phase[dfB$mitosis < dfB$sphase] <- "SPhase"

dfC <- data.frame(time=times, phase=phase, pc1=pca$rotation[,1],
    pc2=pca$rotation[,2])

fig <- ggplot(dfC, aes(x=pc1, y=pc2, color=phase, size=1.5)) + 
    geom_point() + 
    labs(title="PCA by Phase", x="PC1", y="PC2") +
    theme_classic() +
    scale_color_manual(values=c("blue", "red", "orange"))
ggsave(filename="Fig2c_PHASE.pdf", plot=fig)
