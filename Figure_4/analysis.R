library(ggplot2)
library(methods)
library(CancerInSilico)
load('../Figure_4/Figure_4_cleaned.RData') #fig4Data
data(SamplePathways)

getMeanScale <- function(mod, pwy)
{
    sapply(0:mod@runTime, function(t)
    {
        mean(sapply(1:getNumberOfCells(mod, t), function(c)
        {
            pwy@expressionScale(mod, c, t)
        }))
    })
}

## Figure 4b - pathway activity for fitted run

pwyActivity <- data.frame(
    activationGrowth=getMeanScale(fig4Data[[1]], pwyGrowth)
    GtoS=getMeanScale(fig4Data[[1]], pwySPhase)
    mitosis=getMeanScale(fig4Data[[1]], pwyMitosis)
    contactInhibition=getMeanScale(fig4Data[[1]], pwyContactInhibition))




