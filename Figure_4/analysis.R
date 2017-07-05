library(ggplot2)
library(gplots)
library(methods)
library(CancerInSilico)
load('../Figure_4/Figure_4_cleaned.RData') #fig4Data
data(SamplePathways)

sampFreq <- 24
nCells <- 10
hours <- seq(0,144,sampFreq)
ndx <- 1:length(hours)

## Figure 4d - gene expression data

pwyGrowth <- calibratePathway(pwyGrowth, referenceGeneExpression)
pwyMitosis <- calibratePathway(pwyMitosis, referenceGeneExpression)
pwySPhase <- calibratePathway(pwySPhase, referenceGeneExpression)
pwyContactInhibition <- calibratePathway(pwyContactInhibition, referenceGeneExpression)
allPwys <- c(pwyGrowth, pwyMitosis, pwySPhase, pwyContactInhibition)
ge_pbs <- inSilicoGeneExpression(fig4Data[[1]], allPwys[c(1,2,3,4)], nCells=nCells, sampFreq=sampFreq, perError=0.05)
ge_10ug <- inSilicoGeneExpression(fig4Data[[2]], allPwys[c(1,2,3,4)], nCells=nCells, sampFreq=sampFreq, perError=0.05)
ge_100ug <- inSilicoGeneExpression(fig4Data[[3]], allPwys[c(1,2,3,4)], nCells=nCells, sampFreq=sampFreq, perError=0.05)
ge_pbs$expression <- ge_pbs$expression[,ndx]
ge_10ug$expression <- ge_10ug$expression[,ndx]
ge_100ug$expression <- ge_100ug$expression[,ndx]

#png(file='fig4d.png')
#heatmap.2(ge$expression, col=greenred, scale='row', trace='none',
#    hclust=function(x) hclust(x,method='complete'),
#    distfun=function(x) as.dist((1-cor(t(x)))/2),
#    Colv=F)

#png(file='fig4d_growth.png')
#heatmap.2(ge$expression[pwyGrowth@genes,], col=greenred, scale='row', trace='none',
#    hclust=function(x) hclust(x,method='complete'),
#    distfun=function(x) as.dist((1-cor(t(x)))/2),
#    Colv=F)

#png(file='fig4d_mitosis.png')
#heatmap.2(ge$expression[pwyMitosis@genes,], col=greenred, scale='row', trace='none',
#    hclust=function(x) hclust(x,method='complete'),
#    distfun=function(x) as.dist((1-cor(t(x)))/2),
#    Colv=F)

#png(file='fig4d_sphase.png')
#heatmap.2(ge$expression[pwySPhase@genes,], col=greenred, scale='row', trace='none',
#    hclust=function(x) hclust(x,method='complete'),
#    distfun=function(x) as.dist((1-cor(t(x)))/2),
#    Colv=F)

#png(file='fig4d_contactInhibition.png')
#heatmap.2(ge$expression[pwyContactInhibition@genes,], col=greenred, scale='row', trace='none',
#    hclust=function(x) hclust(x,method='complete'),
#    distfun=function(x) as.dist((1-cor(t(x)))/2),
#    Colv=F)

## Figure 4c - pathway activity for fitted simulations

movAvg <- function(data)
{
    window <- 24 / sampFreq
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
    activationGrowth_pbs    =   ge_pbs$pathways[[1]][ndx],
    mitosis_pbs             =   movAvg(ge_pbs$pathways[[2]][ndx]),
    GtoS_pbs                =   movAvg(ge_pbs$pathways[[3]][ndx]),
    contactInhibition_pbs   =   ge_pbs$pathways[[4]][ndx],

    activationGrowth_10ug   =   ge_10ug$pathways[[1]][ndx],
    mitosis_10ug            =   movAvg(ge_10ug$pathways[[2]][ndx]),
    GtoS_10ug               =   movAvg(ge_10ug$pathways[[3]][ndx]),
    contactInhibition_10ug  =   ge_10ug$pathways[[4]][ndx],

    activationGrowth_100ug  =   ge_100ug$pathways[[1]][ndx],
    mitosis_100ug           =   movAvg(ge_100ug$pathways[[2]][ndx]),
    GtoS_100ug              =   movAvg(ge_100ug$pathways[[3]][ndx]),
    contactInhibition_100ug =   ge_100ug$pathways[[4]][ndx])

fig <- ggplot(pwyActivity, aes(x=hour)) +

#    geom_line(aes(y=GtoS_pbs,               color='pbs')) +
    geom_line(aes(y=activationGrowth_pbs,   color='pbs')) +
#    geom_line(aes(y=mitosis_pbs,            color='pbs')) +
    geom_line(aes(y=contactInhibition_pbs,  color='pbs')) +

#    geom_line(aes(y=GtoS_10ug,              color='10ug')) +
    geom_line(aes(y=activationGrowth_10ug,  color='10ug')) +
#    geom_line(aes(y=mitosis_10ug,           color='10ug')) +
    geom_line(aes(y=contactInhibition_10ug, color='10ug')) +

#    geom_line(aes(y=GtoS_100ug,             color='100ug')) +
    geom_line(aes(y=activationGrowth_100ug, color='100ug')) +
#    geom_line(aes(y=mitosis_100ug,          color='100ug'))
    geom_line(aes(y=contactInhibition_100ug,color='100ug'))

ggsave(filename='fig4c.png', plot=fig)
