library(CancerInSilico)

## read data, extract neccesary info
allFiles <- list.files(path='/mnt/c/Users/tomsh/Documents/FertigLab/CancerInSilico/fig4_data',
    full.names = TRUE, recursive = TRUE, pattern = "*.RData")

fig4Data <- list()
for (file in allFiles)
{
    load(file)
    fig4Data[[file]] <- output
}

params <- new('GeneExpressionParams')
params@RNAseq <- FALSE
params@singleCell <- FALSE
params@nCells <- 50
params@sampleFreq <- 4
params@perError <- 0.1

hours <- seq(24,168,params@sampleFreq)
colNdx <- hours / params@sampleFreq + 1

# load default pathways
data(SamplePathways)

# add gene overlap to growth pathway
N <- length(pwyGrowth@genes)
M_overlap <- sample(pwyGrowth@genes, N/2, replace=FALSE)
S_overlap <- setdiff(pwyGrowth@genes, M_overlap)
M_overlap <- pwyGrowth@genes
S_overlap <- pwyGrowth@genes
#overlap <- c(pwyMitosis@genes, pwySPhase@genes, pwyGrowth@genes)
#pwyMitosis@genes <- overlap
#pwySPhase@genes <- overlap
# pwyGrowth@genes <- overlap

# set expression levels
pwyContactInhibition <- calibratePathway(pwyContactInhibition,
    referenceGeneExpression)
pwyGrowth <- calibratePathway(pwyGrowth, referenceGeneExpression)
pwyMitosis <- calibratePathway(pwyMitosis, referenceGeneExpression)
pwySPhase <- calibratePathway(pwySPhase, referenceGeneExpression)

allPwys <- c(pwyGrowth, pwyMitosis, pwySPhase, pwyContactInhibition)

ge_pbs <- inSilicoGeneExpression(fig4Data[[1]], allPwys, params)
ge_100ug <- inSilicoGeneExpression(fig4Data[[3]], allPwys, params)

ge_pbs$expression <- ge_pbs$expression[,colNdx]
ge_100ug$expression <- ge_100ug$expression[,colNdx]

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

pwyActivity <- data.frame(hour=hours,
    activationGrowth_pbs    =   ge_pbs$pathways[[1]][colNdx],
    mitosis_pbs             =   movAvg(ge_pbs$pathways[[2]][colNdx]),
    GtoS_pbs                =   movAvg(ge_pbs$pathways[[3]][colNdx]),
    contactInhibition_pbs   =   ge_pbs$pathways[[4]][colNdx],

    activationGrowth_100ug  =   ge_100ug$pathways[[1]][colNdx],
    mitosis_100ug           =   movAvg(ge_100ug$pathways[[2]][colNdx]),
    GtoS_100ug              =   movAvg(ge_100ug$pathways[[3]][colNdx]),
    contactInhibition_100ug =   ge_100ug$pathways[[4]][colNdx]
)

growthGenes <- pwyGrowth@genes
mitosisGenes <- pwyMitosis@genes
SPhaseGenes <- pwySPhase@genes
contactInhibitionGenes <- pwyContactInhibition@genes

save(ge_pbs, ge_100ug, pwyActivity, growthGenes, mitosisGenes,
    SPhaseGenes, contactInhibitionGenes, file='Figure_4_cleaned.RData')

