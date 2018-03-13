# load libraries
library(CancerInSilico)
data(SamplePathways)
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
    params@nCells <- 50
    params@sampleFreq <- 24
    params@perError <- 0.1

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
ge <- foreach(i = 1:24, .packages="CancerInSilico") %dopar%
{
    if (i <= 2)
    {
        j <- ifelse(i == 1, 1, 3)
        simExpression(cellModels[[j]], tangExpressionData, "none")
    }
    else if (i <= 4)
    {
        j <- ifelse(i == 3, 1, 3)
        simExpression(cellModels[[j]], tangExpressionData, "half")
    }
    else if (i <= 6)
    {
        j <- ifelse(i == 3, 1, 3)
        simExpression(cellModels[[j]], tangExpressionData, "full")
    }
    else if (i <= 12)
    {
        simExpression(cellModels[[i-3]], STC_transformed, "none", TRUE)
    }
    else if (i <= 18)
    {
        simExpression(cellModels[[i-9]], STC_transformed, "half", TRUE)
    }
    else
    {
        simExpression(cellModels[[i-15]], STC_transformed, "full", TRUE)
    }
}
stopCluster(cl)

# combine raw data
tang_none <- cbind(ge[[1]]$expression, ge[[2]]$expression)
tang_half <- cbind(ge[[3]]$expression, ge[[4]]$expression)
tang_full <- cbind(ge[[5]]$expression, ge[[6]]$expression)

kagohara_none <- cbind(ge[[7]]$expression, ge[[8]]$expression,
    ge[[9]]$expression, ge[[10]]$expression, ge[[11]]$expression,
    ge[[12]]$expression)
kagohara_half <- cbind(ge[[13]]$expression, ge[[14]]$expression,
    ge[[15]]$expression, ge[[16]]$expression, ge[[17]]$expression,
    ge[[18]]$expression)
kagohara_full <- cbind(ge[[19]]$expression, ge[[20]]$expression,
    ge[[21]]$expression, ge[[22]]$expression, ge[[23]]$expression,
    ge[[24]]$expression)

# needed for smoothing pathwat activity
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

save(tang_none, tang_half, tang_full, kagohara_none, kagohara_half,
    kagohara_full, tangPwyActivity, M_half, S_half,
    file='Figure_4_cleaned.RData')

