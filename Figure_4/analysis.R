library(ggplot2)
library(gplots)
library(methods)
library(CancerInSilico)
load('../Figure_4/Figure_4_cleaned.RData') #fig4Data
data(SamplePathways)

## Figure 4a - gene expression data
pwyGrowth <- calibratePathway(pwyGrowth, referenceGeneExpression)
pwyMitosis <- calibratePathway(pwyMitosis, referenceGeneExpression)
pwySPhase <- calibratePathway(pwySPhase, referenceGeneExpression)
pwyContactInhibition <- calibratePathway(pwyContactInhibition, referenceGeneExpression)
allPwys <- c(pwyGrowth, pwyMitosis, pwySPhase, pwyContactInhibition)
ge <- inSilicoGeneExpression(fig4Data[[1]], allPwys[c(1,2,3,4)], sampFreq=24, perError=0.05)

png(file='fig4a.png')
heatmap.2(ge$expression, col=greenred, scale='row', trace='none',
    hclust=function(x) hclust(x,method='complete'),
    distfun=function(x) as.dist((1-cor(t(x)))/2),
    Colv=F)

png(file='fig4a_growth.png')
heatmap.2(ge$expression[pwyGrowth@genes,], col=greenred, scale='row', trace='none',
    hclust=function(x) hclust(x,method='complete'),
    distfun=function(x) as.dist((1-cor(t(x)))/2),
    Colv=F)

png(file='fig4a_mitosis.png')
heatmap.2(ge$expression[pwyMitosis@genes,], col=greenred, scale='row', trace='none',
    hclust=function(x) hclust(x,method='complete'),
    distfun=function(x) as.dist((1-cor(t(x)))/2),
    Colv=F)

png(file='fig4a_sphase.png')
heatmap.2(ge$expression[pwySPhase@genes,], col=greenred, scale='row', trace='none',
    hclust=function(x) hclust(x,method='complete'),
    distfun=function(x) as.dist((1-cor(t(x)))/2),
    Colv=F)

png(file='fig4a_contactInhibition.png')
heatmap.2(ge$expression[pwyContactInhibition@genes,], col=greenred, scale='row', trace='none',
    hclust=function(x) hclust(x,method='complete'),
    distfun=function(x) as.dist((1-cor(t(x)))/2),
    Colv=F)

#dev.off()

#heatmap.2(D.gene.avg[gp, paste(rep(0:14,2), rep(c('Control','Cetuximab'),each=15),sep=".")],
    #col=greenred,scale='row',
    #trace='none',
    #hclust=function(x) hclust(x,method="complete"),
    #distfun=function(x) as.dist((1-cor(t(x)))/2),
    #Colv=F,
    #labRow="", main=p,
    #ColSideColors = rep(c('black','red'),each=15),
    #labCol = sprintf('%0.1f', rep(c(0,5/60,10/60,30/60, 1,2,4,8,24,48,72,96,120,144,168),2)),
    #xlab='time (hr)')

## Figure 4b - pathway activity for fitted run

pwyActivity <- data.frame(hour=seq(0,168,24),
    activationGrowth=ge$pathways[[1]],
    mitosis=ge$pathways[[2]],
    GtoS=ge$pathways[[3]],
    contactInhibition=ge$pathways[[4]])

fig <- ggplot(pwyActivity, aes(x=hour)) + geom_line(aes(y=GtoS)) +
    geom_line(aes(y=activationGrowth)) + geom_line(aes(y=mitosis)) + 
    geom_line(aes(y=contactInhibition))
ggsave(filename='fig4b.png', plot=fig)

## Figure 4d - comparison of real gene expression data
