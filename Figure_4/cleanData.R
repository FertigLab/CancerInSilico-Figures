# load libraries
library(CancerInSilico)
data(SamplePathways)

# load the fitted simulations from figure 3
load("../Figure_3/Fig3_Fitted_Simulations.RData")

# STC.D PBS/CTX @ 0,24,48,72,96,120 hours
load("../../DataFromKagohara/STCCountsCG.Rda")

# tangExpressionData @ 0,4,8,24,48,72,96,120,144
load(".././DataFromTang/TangExpressionData.Rda")


# for sampling done prior to main functions which set their own seed internally
set.seed(123)

# run cell model
runModel <- function(sim)
{
    type <- new('CellType', name='DEFAULT', minCycle=sim$cycleLength,
        cycleLength=function() sim$cycleLength)
    drug <- new('Drug', name='DEFAULT', timeAdded=24,
        cycleLengthEffect=function(a,b) rnorm(n=1, mean=b*sim$drugEffect, sd=4))
    inSilicoCellModel(initialNum=100, runTime=168, density=sim$initDensity,
        boundary=1, syncCycle=FALSE, randSeed=0, outputIncrement=4,
        recordIncrement=0.25, timeIncrement=0.001, cellTypes=c(type),
        cellTypeInitFreq=c(1), drugs=c(drug), maxDeformation=0.1,
        maxTranslation=0.1, maxRotation=0.3, nG=28, epsilon=10, delta=0.2)
}

# generate gene expression data
params <- new('GeneExpressionParams')
params@RNAseq <- RNAseq
params@singleCell <- FALSE
params@nCells <- 50
params@sampleFreq <- 4
params@perError <- 0.1
simExpression <- function(refData, overlap='none', RNAseq=FALSE)
{
    params@RNAseq <- RNAseq

    if (overlap == 'half')
    {
        N <- length(pwyGrowth@genes)
        M_overlap <- sample(pwyGrowth@genes, N/2, replace=FALSE)
        S_overlap <- setdiff(pwyGrowth@genes, M_overlap)
        pwyMitosis@genes <- c(pwyMitosis@genes, M_overlap)
        pwySPhase@genes <- c(pwySPhase@genes, S_overlap)
    }
    else if (ovelap == 'full')
    {
        overlap <- c(pwyMitosis@genes, pwySPhase@genes, pwyGrowth@genes)
        pwyMitosis@genes <- overlap
        pwySPhase@genes <- overlap
        pwyGrowth@genes <- overlap
    }

    pwyContactInhibition <- calibratePathway(pwyContactInhibition, refData)
    pwyGrowth <- calibratePathway(pwyGrowth, refData)
    pwyMitosis <- calibratePathway(pwyMitosis, refData)
    pwySPhase <- calibratePathway(pwySPhase, refData)
    allPwys <- c(pwyGrowth, pwyMitosis, pwySPhase, pwyContactInhibition)
}

# re-run the model to get the full cell data

tang_pbs_mod <- runModel(tangFit[['pbs']])
tang_10ug_mod <- runModel(tangFit[['10ug']])
tang_100ug_mod <- runModel(tangFit[['100ug']])
kagohara_pbs_scc1_mod <- runModel(kagoharaFit[['pbs_scc1']])
kagohara_pbs_scc6_mod <- runModel(kagoharaFit[['pbs_scc6']])
kagohara_pbs_scc25_mod <- runModel(kagoharaFit[['pbs_scc25']])
kagohara_ctx_scc1_mod <- runModel(kagoharaFit[['pbs_scc1']])
kagohara_ctx_scc6_mod <- runModel(kagoharaFit[['ctx_scc6']])
kagohara_ctx_scc25_mod <- runModel(kagoharaFit[['ctx_scc25']])

# generate expression data

tang_pbs_none <- simExpression(tang_pbs_mod, D.gene, "none")
tang_pbs_full <- simExpression(tang_pbs_mod, D.gene, "full")
tang_100ug_none <- simExpression(tang_100ug_mod, D.gene, "none")
tang_100ug_full <- simExpression(tang_100ug_mod, D.gene, "full")

tang_none <- cbind(tang_pbs_none$expression, tang_100ug_none$expression)
tang_full <- cbind(tang_pbs_full$expression, tang_100ug_full$expression)

STC_transformed <- floor(2^STC.D)
kagohara_pbs_scc1_none <- simExpression(kagohara_pbs_scc1_mod, STC_transformed, "none", TRUE)
kagohara_pbs_scc6_none <- simExpression(kagohara_pbs_scc6_mod, STC_transformed, "none", TRUE)
kagohara_pbs_scc25_none <- simExpression(kagohara_pbs_scc25_mod, STC_transformed, "none", TRUE)
kagohara_ctx_scc1_none <- simExpression(kagohara_ctx_scc1_mod, STC_transformed, "none", TRUE)
kagohara_ctx_scc6_none <- simExpression(kagohara_ctx_scc6_mod, STC_transformed, "none", TRUE)
kagohara_ctx_scc25_none <- simExpression(kagohara_ctx_scc25_mod, STC_transformed, "none", TRUE)

kagohara_none <- cbind(
    kagohara_pbs_scc1_none$expression,
    kagohara_pbs_scc6_none$expression,
    kagohara_pbs_scc25_none$expression,
    kagohara_ctx_scc1_none$expression,
    kagohara_ctx_scc6_none$expression,
    kagohara_ctx_scc25_none$expression
)

kagohara_pbs_scc1_full <- simExpression(kagohara_pbs_scc1_mod, STC_transformed, "full", TRUE)
kagohara_pbs_scc6_full <- simExpression(kagohara_pbs_scc6_mod, STC_transformed, "full", TRUE)
kagohara_pbs_scc25_full <- simExpression(kagohara_pbs_scc25_mod, STC_transformed, "full", TRUE)
kagohara_ctx_scc1_full <- simExpression(kagohara_ctx_scc1_mod, STC_transformed, "full", TRUE)
kagohara_ctx_scc6_full <- simExpression(kagohara_ctx_scc6_mod, STC_transformed, "full", TRUE)
kagohara_ctx_scc25_full <- simExpression(kagohara_ctx_scc25_mod, STC_transformed, "full", TRUE)

kagohara_full <- cbind(
    kagohara_pbs_scc1_full$expression,
    kagohara_pbs_scc6_full$expression,
    kagohara_pbs_scc25_full$expression,
    kagohara_ctx_scc1_full$expression,
    kagohara_ctx_scc6_full$expression,
    kagohara_ctx_scc25_full$expression
)

# subset gene expression data to match real data sample points

# needed for smoothing pathwat activity
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

# only do for tang - overlap does not matter
tangPwyActivity <- data.frame(hour=hours,
    activationGrowth_pbs    =   tang_pbs_none$pathways[[1]],
    mitosis_pbs             =   movAvg(tang_pbs_none$pathways[[2]]),
    GtoS_pbs                =   movAvg(tang_pbs_none$pathways[[3]]),
    contactInhibition_pbs   =   tang_pbs_none$pathways[[4]],

    activationGrowth_100ug  =   tang_100ug_none$pathways[[1]],
    mitosis_100ug           =   movAvg(tang_100ug_none$pathways[[2]]),
    GtoS_100ug              =   movAvg(tang_100ug_none$pathways[[3]]),
    contactInhibition_100ug =   tang_100ug_none$pathways[[4]]
)

save(tang_none, tang_full, kagohara_none, kagohara_full, tangPwyActivity,
    file='Figure_4_cleaned.RData')

