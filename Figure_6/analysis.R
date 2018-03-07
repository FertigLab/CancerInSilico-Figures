library(ggplot2)
library(gplots)
library(CancerInSilico)
library(scater)
library(methods)
load('Figure_6_cleaned.RData') # pwyActivity, ge, ge_bulk, cellPhase, cellType

set.seed(123)

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

getTime <- function(str)
{
  pos <- regexpr('_', str)[1]
  time <- as.numeric(substr(str, pos + 2, nchar(str)))
}

getSCE <- function(counts)
{
    cell_types <- unname(sapply(colnames(counts), getType))
    cell_phases <- unname(sapply(colnames(counts), getPhase))
    cell_times <- unname(sapply(colnames(counts), getTime))
    I_phase <- which(cell_phases == 'I')
    S_phase <- which(cell_phases == 'S')
    M_phase <- which(cell_phases == 'M')

    sampled <- c(sample(I_phase, floor(length(I_phase) / 3), replace=F),
                 sample(S_phase, floor(length(S_phase) / 3), replace=F),
                 sample(M_phase, floor(length(M_phase) / 3), replace=F))

    sampled <- sample(1:ncol(counts), floor(ncol(counts) / 10), replace=F)
    counts <- counts[,sampled]
    cell_types <- cell_types[sampled]
    cell_phases <- cell_phases[sampled]
    cell_times <- cell_times[sampled]

    phenos <- new('AnnotatedDataFrame', data = data.frame(Cell = colnames(counts),
        Type = cell_types, Phase = cell_phases, Time = cell_times))
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

pdf(file='fig6.pdf', width=11, height=8.5)
scater::plotPCA(sce, exprs_values = 'counts', shape_by = "Type",
  colour_by = "Phase", ncomponents = 2)
warnings()