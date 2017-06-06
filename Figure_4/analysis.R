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

pwyActivity <- data.frame(hour=0:168,
    activationGrowth=getMeanScale(fig4Data[[1]], pwyGrowth),
    GtoS=getMeanScale(fig4Data[[1]], pwySPhase),
    mitosis=getMeanScale(fig4Data[[1]], pwyMitosis),
    contactInhibition=1 - getMeanScale(fig4Data[[1]], pwyContactInhibition))

fig <- ggplot(pwyActivity, aes(x=hour)) + geom_line(aes(y=GtoS)) +
    geom_line(aes(y=activationGrowth)) + geom_line(aes(y=mitosis)) + 
    geom_line(aes(y=contactInhibition))
ggsave(filename='fig4b.png', plot=fig)
