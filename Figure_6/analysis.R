library(ggplot2)
library(gplots)
library(CancerInSilico)
library(methods)
load('Figure_6_cleaned.RData') #ge_pbs, ge_10ug, ge_100ug, pwyActivity
data(SamplePathways)

#png(file='fig6d.png')
#heatmap.2(ge_100ug_bulk$expression, col=greenred, scale='row', trace='none',
#    hclust=function(x) hclust(x,method='complete'),
#    distfun=function(x) as.dist((1-cor(t(x)))/2), Colv=F, dendrogram='row')

#png(file='fig6d_growth.png')
#heatmap.2(ge_100ug_bulk$expression[pwyGrowth@genes,], col=greenred, scale='row',
#    trace='none', hclust=function(x) hclust(x,method='complete'),
#    distfun=function(x) as.dist((1-cor(t(x)))/2), Colv=F, dendrogram='row')

#png(file='fig6d_mitosis.png')
#heatmap.2(ge_100ug_bulk$expression[pwyMitosis@genes,], col=greenred, scale='row',
#    trace='none', hclust=function(x) hclust(x,method='complete'),
#    distfun=function(x) as.dist((1-cor(t(x)))/2), Colv=F, dendrogram='row')

#png(file='fig6d_sphase.png')
#heatmap.2(ge_100ug_bulk$expression[pwySPhase@genes,], col=greenred, scale='row',
#    trace='none', hclust=function(x) hclust(x,method='complete'),
#    distfun=function(x) as.dist((1-cor(t(x)))/2), Colv=F, dendrogram='row')

#png(file='fig6d_contactInhibition.png')
#heatmap.2(ge_100ug_bulk$expression[pwyContactInhibition@genes,], col=greenred, scale='row',
#    trace='none', hclust=function(x) hclust(x,method='complete'),
#    distfun=function(x) as.dist((1-cor(t(x)))/2), Colv=F, dendrogram='row')
