# load libraries and data
library(ggplot2)
library(gplots)
library(CancerInSilico)
load('Figure_4_cleaned.RData')
# STC.D PBS/CTX @ 0,24,48,72,96,120 hours
load("../../DataFromKagohara/STCCountsCG.Rda")
# tangExpressionData @ 0,4,8,24,48,72,96,120,144
load("../../DataFromTang/TangExpressionData.Rda")
data(SamplePathways)
allGenes <- c(pwyGrowth@genes, pwyMitosis@genes, pwySPhase@genes, 
    pwyContactInhibition@genes)

# helper functions

internalHeatmap <- function(filename, title, data, times, rowColors, colColors)
{
    pdf(filename, width=13, height=10)

    ndx <- apply(data, 1, var) == 0 # remove zero variance rows

    heatmap.2(data[!ndx,], 
        col=greenred, scale='row',
        trace='none', hclust=function(x) hclust(x,method='complete'),
        distfun=function(x) as.dist((1-cor(t(x)))/2), 
        Colv=FALSE, dendrogram='row',
        ColSideColors = colColors,
        RowSideColors = rowColors[!ndx],
        labCol = times,
        main=title)
    dev.off()
}

getGeneNames <- function(pathway, overlap)
{
    g <- pwyGrowth@genes
    ci <- pwyContactInhibition@genes
    m <- switch(overlap, none=pwyMitosis@genes, half=M_half, full=c(pwyMitosis@genes, g))
    s <- switch(overlap, none=pwySPhase@genes, half=S_half, full=c(pwySPhase@genes, g))
    
    switch(pathway, growth=g, mitosis=m, sphase=s, contact_inhibition=ci)
}

plotTangHeatmap <- function(data, simulated, pathway, overlap)
{
    print(paste(simulated, pathway, overlap))
    filename <- paste('tang_', simulated, '_', pathway, '_', overlap, '.pdf',
        sep='')
    plotTitle <- paste('tang ', simulated, ' data - ', pathway, ' genes ',
        switch(overlap, none="", half=", half growth overlap with M/S",
        full=", full growth overlap with M/S"))
    genes <- getGeneNames(pathway, overlap)
    geneColors <- ifelse(rownames(data) %in% genes, 'black', 'white')
    times <- rep(seq(0,144,24),2)
    drugColors <- rep(c('black','red'), each=7)

    internalHeatmap(filename, plotTitle, data, times, geneColors, drugColors)
}

plotKagoharaHeatmap <- function(data, simulated, pathway, overlap, cellLine)
{
    print(paste(simulated, pathway, overlap, cellLine))
    filename <- paste('kagohara_', cellLine, '_', simulated, '_', pathway, '_',
        overlap, '.pdf', sep='')
    plotTitle <- paste('kagohara', cellLine, simulated, 'data -', pathway,
        'genes', switch(overlap, none="", half=", half growth overlap with M/S",
        full=", full growth overlap with M/S"))
    genes <- getGeneNames(pathway, overlap)
    geneColors <- ifelse(rownames(data) %in% genes, 'black', 'white')

    if (cellLine == 'all')
    {
        times <- c(rep(seq(0,120,24),3), rep(seq(24,120,24),3))
        drugColors <- c(rep('black',18), rep('red',15))
    }
    else
    {
        times <- c(seq(0,120,24), seq(24,120,24))
        drugColors <- c(rep('black',6), rep('red',5))
    }

    internalHeatmap(filename, plotTitle, data, times, geneColors, drugColors)
}

plotTangPathways <- function(simulated, overlap)
{
    if (simulated == 'simulated')
    {
        mat <- switch(overlap, none=tang_none, half=tang_half, full=tang_full)
        mat <- mat[,c(1:7,9:15)]
    }
    else
    {
        mat <- tangExpressionData[,c(14,6,8,10,12,2,4,13,5,7,9,11,1,3)]
        mat <- mat[rownames(mat) %in% allGenes,]
    }
    print(colnames(mat))
    plotTangHeatmap(mat, simulated, 'growth', overlap)
    plotTangHeatmap(mat, simulated, 'mitosis', overlap)
    plotTangHeatmap(mat, simulated, 'sphase', overlap)
    plotTangHeatmap(mat, simulated, 'contact_inhibition', overlap)
}

plotKagoharaPathways <- function(simulated, overlap, cellLine)
{
    # drug is added at t=24 in simulation, simulated is 24 hours ahead of real
    if (simulated == 'simulated')
    {
        ndx <- switch(cellLine, c(1:6,17:22,33:38,11:15,27:31,43:47),
            scc1=c(1:6,11:15), scc6=c(17:22,27:31), scc25=c(33:38,43:47))
        mat <- switch(overlap, none=kagohara_none, half=kagohara_half,
            full=kagohara_full)
        mat <- mat[,ndx]
        mat <- log2(mat + 1)
    }
    else
    {
        ndx <- switch(cellLine, all=c(1:6,12:17,23:28,7:11,18:22,29:33),
            scc1=1:11, scc6=12:22, scc25=23:33)
        mat <- STC.D[,ndx]
        mat <- mat[rownames(mat) %in% allGenes,]
    }
    print(colnames(mat))
    plotKagoharaHeatmap(mat, simulated, 'growth', overlap, cellLine)
    plotKagoharaHeatmap(mat, simulated, 'mitosis', overlap, cellLine)
    plotKagoharaHeatmap(mat, simulated, 'sphase', overlap, cellLine)
    plotKagoharaHeatmap(mat, simulated, 'contact_inhibition', overlap, cellLine)
}

plotKagoharaCellLine <- function(cellLine)
{
    plotKagoharaPathways('simulated', 'none', cellLine)
    plotKagoharaPathways('simulated', 'half', cellLine)
    plotKagoharaPathways('simulated', 'full', cellLine)
    plotKagoharaPathways('real', 'none', cellLine)
}

# create heatmaps

plotKagoharaCellLine('all')
plotKagoharaCellLine('scc1')
plotKagoharaCellLine('scc6')
plotKagoharaCellLine('scc25')

plotTangPathways('simulated', 'none')
plotTangPathways('simulated', 'half')
plotTangPathways('simulated', 'full')
plotTangPathways('real', 'none')

## plot pathway activity

fig <- ggplot(tangPwyActivity, aes(x=hour)) +
    geom_line(aes(y=activationGrowth_100ug, color='Growth Rate', linetype='100ug')) +
    geom_line(aes(y=mitosis_100ug,          color='Mitosis', linetype='100ug')) +
    geom_line(aes(y=contactInhibition_100ug,color='Contact Inhibition', linetype='100ug')) +
    geom_line(aes(y=activationGrowth_pbs,   color='Growth Rate', linetype='PBS')) +
    geom_line(aes(y=mitosis_pbs,            color='Mitosis', linetype='PBS')) +
    geom_line(aes(y=contactInhibition_pbs,  color='Contact Inhibition', linetype='PBS')) +
    scale_linetype_manual(values = c('dotted', 'solid')) +
    labs(title = "Pathway Activity", caption = "Figure 4", x = "Hour",
        y = "Percent Pathway Activity", color = "Pathway",
        linetype = "Dosage")
ggsave(filename='fig4_pathway_activity.pdf', plot=fig)


