# load libraries
library(CancerInSilico)
library(foreach)
library(doParallel)

# load cell models
load("Fig4_RawData.RData")

# STC.D PBS/CTX @ 0,24,48,72,96,120 hours
load("../../DataFromKagohara/STCCountsCG.Rda")

# tangExpressionData @ 0,4,8,24,48,72,96,120,144
load("../../DataFromTang/TangExpressionData.Rda")

# for sampling done prior to main functions which set their own seed internally
set.seed(123)

# create pathways for various overlap
data(SamplePathways)
M_none <- pwyMitosis@genes
S_none <- pwySPhase@genes
G_none <- pwyGrowth@genes
CI_none <- pwyContactInhibition@genes

M_half <- c(M_none, sample(G_none, length(G_none) / 2, replace=FALSE))
S_half <- c(S_none, setdiff(G_none, M_half))
G_half <- G_none
CI_half <- CI_none

M_full <- c(M_none, G_none)
S_full <- c(S_none, G_none)
G_full <- G_none
CI_full <- CI_none

# generate gene expression data
simExpression <- function(model, refData, overlap='none', RNAseq=FALSE)
{
    params <- new('GeneExpressionParams')
    params@RNAseq <- RNAseq
    params@singleCell <- FALSE
    params@nCells <- 100
    params@sampleFreq <- 24
    params@perError <- 0.1
    params@randSeed <- 54

    M_genes <- switch(overlap, none=M_none, half=M_half, full=M_full)
    S_genes <- switch(overlap, none=S_none, half=S_half, full=S_full)
    G_genes <- switch(overlap, none=G_none, half=G_half, full=G_full)
    CI_genes <- switch(overlap, none=CI_none, half=CI_half, full=CI_full)

    pwyMitosis@genes <- M_genes[M_genes %in% rownames(refData)]
    pwySPhase@genes <- S_genes[S_genes %in% rownames(refData)]
    pwyGrowth@genes <- G_genes[G_genes %in% rownames(refData)]
    pwyContactInhibition@genes <- CI_genes[CI_genes %in% rownames(refData)]

    pwyContactInhibition <- calibratePathway(pwyContactInhibition, refData)
    pwyGrowth <- calibratePathway(pwyGrowth, refData)
    pwyMitosis <- calibratePathway(pwyMitosis, refData)
    pwySPhase <- calibratePathway(pwySPhase, refData)
    allPwys <- c(pwyGrowth, pwyMitosis, pwySPhase, pwyContactInhibition)
    inSilicoGeneExpression(model, allPwys, params)
}

# generate expression data
cl <- makeCluster(6)
registerDoParallel(cl)
STC_transformed <- floor(2^STC.D)
ge <- foreach(i = 1:48, .packages="CancerInSilico") %dopar%
{
    data(SamplePathways)
    if (i <= 4)
    {
        simExpression(cellModels[[i]], tangExpressionData, "none")
    }
    else if (i <= 8)
    {
        simExpression(cellModels[[i-4]], tangExpressionData, "half")
    }
    else if (i <= 12)
    {
        simExpression(cellModels[[i-8]], tangExpressionData, "full")
    }
    else if (i <= 24)
    {
        simExpression(cellModels[[i-8]], STC_transformed, "none", TRUE)
    }
    else if (i <= 36)
    {
        simExpression(cellModels[[i-20]], STC_transformed, "half", TRUE)
    }
    else
    {
        simExpression(cellModels[[i-32]], STC_transformed, "full", TRUE)
    }
}
stopCluster(cl)

# combine raw data
combineRaw <- function(expList, ndx)
{
    do.call("cbind", lapply(expList[ndx], function(l) l$expression))
}

tang_none_fix <- combineRaw(ge, 1:2)
tang_none_var <- combineRaw(ge, 3:4)
tang_half_fix <- combineRaw(ge, 5:6)
tang_half_var <- combineRaw(ge, 7:8)
tang_full_fix <- combineRaw(ge, 9:10)
tang_full_var <- combineRaw(ge, 11:12)

kagohara_none_fix <- combineRaw(ge, 13:18)
kagohara_none_var <- combineRaw(ge, 19:24)
kagohara_half_fix <- combineRaw(ge, 25:30)
kagohara_half_var <- combineRaw(ge, 31:36)
kagohara_full_fix <- combineRaw(ge, 37:42)
kagohara_full_var <- combineRaw(ge, 43:48)

# needed for smoothing pathway activity
movAvg <- function(data)
{
    window <- 24 / 24 # denom = sample frequency
    avg <- c()    
    for (i in 1:length(data))
    {
        mn <- max(1, i - window)
        mx <- min(length(data), i + window)
        avg <- c(avg, mean(data[mn:mx]))
    }
    return(avg)
}

# only do for tang - overlap does not matter
tangPwyActivity <- data.frame(hour=seq(0,168,24),
    activationGrowth_pbs    =   ge[[1]]$pathways[[1]],
    mitosis_pbs             =   movAvg(ge[[1]]$pathways[[2]]),
    GtoS_pbs                =   movAvg(ge[[1]]$pathways[[3]]),
    contactInhibition_pbs   =   ge[[1]]$pathways[[4]],

    activationGrowth_100ug  =   ge[[2]]$pathways[[1]],
    mitosis_100ug           =   movAvg(ge[[2]]$pathways[[2]]),
    GtoS_100ug              =   movAvg(ge[[2]]$pathways[[3]]),
    contactInhibition_100ug =   ge[[2]]$pathways[[4]]
)

save(tang_none_fix, tang_none_var, tang_half_fix, tang_half_var, tang_full_fix,
    tang_full_var, kagohara_none_fix, kagohara_none_var, kagohara_half_fix,
    kagohara_half_var, kagohara_full_fix, kagohara_full_var,
    tangPwyActivity, M_half, S_half, file='Figure_4_cleaned.RData')


