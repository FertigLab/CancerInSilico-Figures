# TODO plot vs time
# TODO add time course pattern
# color pc 1 & 2 together

############################# load libraries ###################################

library(CancerInSilico)
library(fastICA)
library(viridis)
library(rgl)
library(ggplot2)
set.seed(42)

########################## adjustable parameters ###############################

totalTime <- 168
initNumCells <- 100
initDensity <- 0.05
nGenes <- 50

######################## define multiple cell types ############################

cellTypeA <- new("CellType", name="type.A", size=1, minCycle=16,
    cycleLength=function() 16 + rexp(1, 1/4))

cellTypeB <- new("CellType", name="type.B", size=1, minCycle=16,
    cycleLength=function() 16 + rexp(1, 1/4))

allCellTypes <- list(cellTypeA, cellTypeB)

########################### simulate cell growth ###############################

mod <- inSilicoCellModel(initialNum=initNumCells, runTime=totalTime,
    density=initDensity, outputIncrement=12, cellTypes=allCellTypes,
    cellTypeInitFreq=c(0.5,0.5), syncCycles=FALSE, randSeed=42)

# load("Figure_4_Cell_Model.RData")

############################ helper functions ##################################

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

typeAPwy <- new("Pathway", genes=paste("typeA.", 1:nGenes, sep=""),
    expressionScale=function(model, cell, time)
    {
        ifelse(getCellType(model, time, cell) == 1, 1, 0)
    }
)

typeBPwy <- new("Pathway", genes=paste("typeB.", 1:nGenes, sep=""),
    expressionScale=function(model, cell, time)
    {
        ifelse(getCellType(model, time, cell) == 2, 1, 0)
    }
)

growthPwy <- new("Pathway", genes=paste("growth.", 1:nGenes, sep=""),
    expressionScale=function(model, cell, time)
    {
        min((getCycleLength(model, time, cell) - 16) / 8, 1)
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
params@nCells <- 50
params@sampleFreq <- 12
params@RNAseq <- TRUE
params@singleCell <- TRUE
params@dropoutPresent <- FALSE
params@bcvCommon <- 0.001
params@combineFUN <- mean

calLambda <- 10
calStdDev <- 2

#################################### figure ####################################

whichTimeCourse <- function(times, ica)
{
    whichMax <- 1
    mx <- abs(cor(times, ica$S[,1]))
    for (i in 2:ncol(ica$S))
    {
        if (abs(cor(times, ica$S[,i])) > mx)
        {
            whichMax <- i
            mx <- abs(cor(times, ica$S[,i]))
        }
    }
    return(whichMax)
}

# part A/B

allPwys <- c(typeAPwy, typeBPwy, contactInhibitionPwy)
allPwys <- calibratePathways(allPwys, lambda=calLambda, stddev=calStdDev)
res <- inSilicoGeneExpression(mod, allPwys, params)
sampledTimes <- convertNamesToTimes(colnames(res$expression))
sampledTypes <- convertNamesToCellTypes(colnames(res$expression))

ica <- fastICA(t(res$expression), n.comp=2)
timeCoursePattern <- whichTimeCourse(sampledTimes, ica)

df_AB <- data.frame(time=sampledTimes, type=sampledTypes,
    ic1=ica$S[,timeCoursePattern],
    ic2=ica$S[,setdiff(1:2, timeCoursePattern)]
)

# part A1, ic1 vs time

fig <- ggplot(df_AB, aes(x=time, y=ic1, color=type)) + 
    geom_point() + 
    labs(title="Time Course Component", x="Time", y="ICA component 1") +
    theme_classic() +
    scale_color_manual(values=c("red", "blue"))
ggsave(filename="Fig4A_1.pdf", plot=fig)

# part B1, ic2 vs time

fig <- ggplot(df_AB, aes(x=time, y=ic2, color=type)) + 
    geom_point() + 
    labs(title="Cell Type Component", x="Time", y="ICA component 2") +
    theme_classic() +
    scale_color_manual(values=c("red", "blue"))
ggsave(filename="Fig4B_1.pdf", plot=fig)

# part A2, ic1 vs ic2 colored by time

fig <- ggplot(df_AB, aes(x=ic1, y=ic2, color=time)) + 
    geom_point() + 
    labs(title="Time Course Component",
        x="ICA Component 1", y="ICA Component 2") +
    theme_classic() +
    scale_colour_gradientn(colours=myColorRamp(df_CD$time))
ggsave(filename="Fig4A_2.pdf", plot=fig)

# part B2, ic1 vs ic2 colored by shape

fig <- ggplot(df_AB, aes(x=ic1, y=ic2, color=type)) + 
    geom_point() + 
    labs(title="Cell Type Component",
        x="ICA Component 1", y="ICA Component 2") +
    theme_classic() +
    scale_color_manual(values=c("red", "blue"))
ggsave(filename="Fig4B_2.pdf", plot=fig)

# part C/D

growthPwy@genes <- c(typeAPwy@genes, typeBPwy@genes)
allPwys <- c(growthPwy, typeAPwy, typeBPwy, contactInhibitionPwy)
allPwys <- calibratePathways(allPwys, lambda=calLambda, stddev=calStdDev)
res <- inSilicoGeneExpression(mod, allPwys, params)
sampledTimes <- convertNamesToTimes(colnames(res$expression))
sampledTypes <- convertNamesToCellTypes(colnames(res$expression))

ica <- fastICA(t(res$expression), n.comp=2)
timeCoursePattern <- whichTimeCourse(sampledTimes, ica)

df_CD <- data.frame(time=sampledTimes, type=sampledTypes,
    ic1=ica$S[,timeCoursePattern],
    ic2=ica$S[,setdiff(1:2, timeCoursePattern)]
)

# part C1, ic1 vs time w/ growth overlap

fig <- ggplot(df_CD, aes(x=time, y=ic1, color=type)) + 
    geom_point() + 
    labs(title="Time Course Component, pathway overlap",
        x="Time", y="ICA component 1") +
    theme_classic() + 
    scale_color_manual(values=c("red", "blue"))
ggsave(filename="Fig4C_1.pdf", plot=fig)

# part D1, ic2 vs time w/ growth overlap

fig <- ggplot(df_CD, aes(x=time, y=ic2, color=type)) + 
    geom_point() + 
    labs(title="Cell Type Component, pathway overlap",
        x="Time", y="ICA component 2") +
    theme_classic() + 
    scale_color_manual(values=c("red", "blue"))
ggsave(filename="Fig4D_1.pdf", plot=fig)

# part C2, ic1 vs ic2 colored by time w/ growth overalap

fig <- ggplot(df_CD, aes(x=ic1, y=ic2, color=time)) + 
    geom_point() + 
    labs(title="Time Course Component, pathway overlap",
        x="ICA Component 1", y="ICA Component 2") +
    theme_classic() + 
    scale_colour_gradientn(colours=myColorRamp(df_CD$time))
ggsave(filename="Fig4C_2.pdf", plot=fig)

# part D2, ic1 vs ic2 colored by shape w/ growth overlap

fig <- ggplot(df_CD, aes(x=ic1, y=ic2, color=type)) + 
    geom_point() + 
    labs(title="Cell Type Component, pathway overlap",
        x="ICA Component 1", y="ICA Component 2") +
    theme_classic() + 
    scale_color_manual(values=c("red", "blue"))
ggsave(filename="Fig4D_2.pdf", plot=fig)


# part F, here we show the correlation drop off as the growth pathway
# overlaps with more and more of the cell type pathways

# allOverlaps <- c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1.0)
# cor_vec <- c()
# for (per in allOverlaps)
# {
#     nOverlap <- round(per * length(typeAPwy@genes))
#     growthPwy@genes <- c(typeAPwy@genes[1:nOverlap], typeBPwy@genes[1:nOverlap])
#     allPwys <- c(growthPwy, typeAPwy, typeBPwy, contactInhibitionPwy)
#     allPwys <- calibratePathways(allPwys, lambda=calLambda, stddev=calStdDev)
# 
#     res <- inSilicoGeneExpression(mod, allPwys, params)
#         
#     ica <- fastICA(t(res$expression), n.comp=3)
#     sampledTypes <- convertNamesToCellTypes(colnames(res$expression))
# 
#     vec <- as.numeric(sampledTypes=="type.A")
#     mx <- 0
#     for (i in 1:ncol(ica$S))
#     {
#         mx <- max(mx, abs(cor(ica$S[,i], vec)))
#     }
#     cor_vec <- c(cor_vec, mx)
# }
# 
# pdf("Fig4E.pdf")
# plot(allOverlaps, cor_vec, type="l",
#     main="Identification of cell type as overlap increases",
#     xlab="percent overlap",
#     ylab="correlation")
# dev.off()