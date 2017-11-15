library(ggplot2)
library(gplots)
library(CancerInSilico)
library(scater)
library(methods)
load('Figure_6_cleaned.RData') #ge_pbs, ge_10ug, ge_100ug, pwyActivity

getType <- function(str)
{
  pos <- regexpr('_', str)[1]
  cell <- as.numeric(substr(str, 2, pos-1))
  return(cellType[cell])
}

getPhase <- function(str)
{
  pos <- regexpr('_', str)[1]
  cell <- as.numeric(substr(str, 2, pos-1))
  time <- as.numeric(substr(str, pos + 2, nchar(str)))
  return(cellPhase[cell, time+1])
}

getSCE <- function(counts)
{
    cell_types <- unname(sapply(colnames(counts), getType))
    cell_phases <- unname(sapply(colnames(counts), getPhase))
    phenos <- new('AnnotatedDataFrame', data = data.frame(Cell = colnames(counts),
        Type = cell_types, Phase = cell_phases))
    rownames(phenos) <- colnames(counts)
    features <- new('AnnotatedDataFrame', data = data.frame(Gene = rownames(counts)))
    rownames(features) <- rownames(counts)
    sim <- scater:::newSCESet(countData = counts, phenoData = phenos,
        featureData = features)

    sce <- scater::newSCESet(countData = counts,
        phenoData = new("AnnotatedDataFrame", data = Biobase:::pData(sim)),
        featureData = new("AnnotatedDataFrame", data = Biobase:::fData(sim)))
    
    return(sce)
}

sce <- getSCE(ge$expression)

png(file='fig6.png', width=5600, height=5600, res=300)
scater::plotPCA(sce, exprs_values = 'counts', colour_by = "Type",
    shape_by = "Phase", ncomponents = 3)
warnings()