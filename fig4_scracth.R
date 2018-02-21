#source("https://bioconductor.org/biocLite.R")
#biocLite('ComplexHeatmap')

library(gplots)
library(ComplexHeatmap)
library(circlize)

load("/mnt/c/Users/tomsh/Documents/FertigLab/CancerInSilico/Figures/Figure_4/Figure_4_cleaned.RData")

times <- c(0,4,8,24,48,72,96,120,144)
timeNdx <- c(1,2,3,7,13,19,25,31,37) # actually starts at 24

ge_mat <- cbind(ge_100ug$expression, ge_pbs$expression)
pathwayGenes <- rownames(ge_pbs$expression)

pdf('sim_G_all.pdf', width=11, height=8.5)
heatmap.2(ge_mat[,c(timeNdx,timeNdx+37)], 
          col=greenred, scale='row',
          trace='none', hclust=function(x) hclust(x,method='complete'),
          distfun=function(x) as.dist((1-cor(t(x)))/2), 
          Colv=F, dendrogram='row',
          ColSideColors = rep(c('red','black'),each=length(timeNdx)),
          RowSideColors = ifelse(pathwayGenes %in% growthGenes, 'black', 'white'),
          labCol = rep(times,2),
          main="sim data - growth genes, growth/M/S overlap")
dev.off()

pdf('sim_M_all.pdf', width=11, height=8.5)
heatmap.2(ge_mat[,c(timeNdx,timeNdx+37)], 
          col=greenred, scale='row',
          trace='none', hclust=function(x) hclust(x,method='complete'),
          distfun=function(x) as.dist((1-cor(t(x)))/2), 
          Colv=F, dendrogram='row',
          ColSideColors = rep(c('red','black'),each=length(timeNdx)),
          RowSideColors = ifelse(pathwayGenes %in% mitosisGenes, 'black', 'white'),
          labCol = rep(times,2),
          main="sim data - mitosis genes, growth/M/S overlap")
dev.off()

pdf('sim_S_all.pdf', width=11, height=8.5)
heatmap.2(ge_mat[,c(timeNdx,timeNdx+37)], 
          col=greenred, scale='row',
          trace='none', hclust=function(x) hclust(x,method='complete'),
          distfun=function(x) as.dist((1-cor(t(x)))/2), 
          Colv=F, dendrogram='row',
          ColSideColors = rep(c('red','black'),each=length(timeNdx)),
          RowSideColors = ifelse(pathwayGenes %in% SPhaseGenes, 'black', 'white'),
          labCol = rep(times,2),
          main="sim data - S phase genes, growth/M/S overlap")
dev.off()

pdf('sim_CI_all.pdf', width=11, height=8.5)
heatmap.2(ge_mat[,c(timeNdx,timeNdx+37)], 
          col=greenred, scale='row',
          trace='none', hclust=function(x) hclust(x,method='complete'),
          distfun=function(x) as.dist((1-cor(t(x)))/2), 
          Colv=F, dendrogram='row',
          ColSideColors = rep(c('red','black'),each=length(timeNdx)),
          RowSideColors = ifelse(pathwayGenes %in% contactInhibitionGenes, 'black', 'white'),
          labCol = rep(times,2),
          main="sim data - contact inhibition genes, growth/M/S overlap")
dev.off()

############################################################################

#     load('../../DataFromTang/Tang_Lumi.Rda')
#     
#     minTime <- 0
#     
#     pData <- lumi.D$pheno
#     
#     load('../../DataFromTang/Tang_GeneAvg.Rda')
#     
#     # annotate data so that replicates are averaged and column names are meaningful
#     D.gene.condition <- matrix(NA_real_,
#                                ncol=length(unique(paste(pData$TimePt,pData$Treatment,sep="."))),
#                                nrow=nrow(D.gene), 
#                                dimnames=list(row.names(D.gene),
#                                              unique(paste(pData$TimePt,pData$Treatment,sep="."))))
#     
#     for (i in 1:nrow(D.gene.condition)) {
#       D.gene.condition[i, ] <-  tapply(D.gene[i,],paste(pData$TimePt,pData$Treatment,sep="."),mean)
#     }
#     
#     realTime = c(0, 5/(60*24), 10/(60*24), 30/(60*24), 1/24, 2/24, 4/24, 8/24, 1:7)
#     
#     timesPlot <- which((24 * realTime) %in% times)
#     
#     
#     pdf('real_G.pdf', width=11, height=8.5)
#     heatmap.2(D.gene.condition[pathwayGenes,paste(rep(timesPlot,2),
#                                                   rep(c('Cetuximab','Control'),each=length(timesPlot)),sep='.')], 
#               col=greenred, scale='row',
#               trace='none', hclust=function(x) hclust(x,method='complete'),
#               distfun=function(x) as.dist((1-cor(t(x)))/2), 
#               Colv=F, dendrogram='row',
#               ColSideColors = rep(c('red','black'),each=length(timesPlot)),
#               labCol = rep(24*realTime[timesPlot],2),
#               RowSideColors = ifelse(pathwayGenes %in% growthGenes, 'black', 'white'),
#               main="real data - growth")
#     dev.off()
#     
#     pdf('real_M.pdf', width=11, height=8.5)
#     heatmap.2(D.gene.condition[pathwayGenes,paste(rep(timesPlot,2),
#                                                   rep(c('Cetuximab','Control'),each=length(timesPlot)),sep='.')], 
#               col=greenred, scale='row',
#               trace='none', hclust=function(x) hclust(x,method='complete'),
#               distfun=function(x) as.dist((1-cor(t(x)))/2), 
#               Colv=F, dendrogram='row',
#               ColSideColors = rep(c('red','black'),each=length(timesPlot)),
#               labCol = rep(24*realTime[timesPlot],2),
#               RowSideColors = ifelse(pathwayGenes %in% mitosisGenes, 'black', 'white'),
#               main="real data - mitosis")
#     dev.off()
#     
#     pdf('real_S.pdf', width=11, height=8.5)
#     heatmap.2(D.gene.condition[pathwayGenes,paste(rep(timesPlot,2),
#                                                   rep(c('Cetuximab','Control'),each=length(timesPlot)),sep='.')], 
#               col=greenred, scale='row',
#               trace='none', hclust=function(x) hclust(x,method='complete'),
#               distfun=function(x) as.dist((1-cor(t(x)))/2), 
#               Colv=F, dendrogram='row',
#               ColSideColors = rep(c('red','black'),each=length(timesPlot)),
#               labCol = rep(24*realTime[timesPlot],2),
#               RowSideColors = ifelse(pathwayGenes %in% SPhaseGenes, 'black', 'white'),
#               main="real data - SPhase")
#     dev.off()
#     
#     pdf('real_CI.pdf', width=11, height=8.5)
#     heatmap.2(D.gene.condition[pathwayGenes,paste(rep(timesPlot,2),
#                 rep(c('Cetuximab','Control'),each=length(timesPlot)),sep='.')], 
#               col=greenred, scale='row',
#               trace='none', hclust=function(x) hclust(x,method='complete'),
#               distfun=function(x) as.dist((1-cor(t(x)))/2), 
#               Colv=F, dendrogram='row',
#               ColSideColors = rep(c('red','black'),each=length(timesPlot)),
#               labCol = rep(24*realTime[timesPlot],2),
#               RowSideColors = ifelse(pathwayGenes %in% contactInhibitionGenes, 'black', 'white'),
#               main="real data - contact inhibition")
#     dev.off()
#     