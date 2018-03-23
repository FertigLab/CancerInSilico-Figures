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

plotTangHeatmap <- function(data, simulated, pathway, overlap, cycle)
{
    print(paste(simulated, pathway, overlap, cycle))
    filename <- paste("tang_", simulated, "_", pathway, "_", overlap, "_",
        cycle, ".pdf", sep="")
    plotTitle <- paste("tang ", simulated, " data - ", pathway, " genes ",
        switch(overlap, none="", half=", half growth overlap with M/S",
        full=", full growth overlap with M/S"),
        ifelse(cycle!="fix", "with cell cycle variance", ""))
    genes <- getGeneNames(pathway, overlap)
    geneColors <- ifelse(rownames(data) %in% genes, "black", "white")
    times <- rep(seq(0,144,24),2)
    drugColors <- rep(c("black","red"), each=7)

    internalHeatmap(filename, plotTitle, data, times, geneColors, drugColors)
}

plotKagoharaHeatmap <- function(data, simulated, pathway, overlap, cellLine, cycle)
{
    print(paste(simulated, pathway, overlap, cellLine, cycle))
    filename <- paste("kagohara_", cellLine, "_", simulated, "_", pathway, "_",
        overlap, "_", cycle, ".pdf", sep="")
    plotTitle <- paste("kagohara", cellLine, simulated, "data -", pathway,
        "genes", switch(overlap, none="", half=", half growth overlap with M/S",
        full=", full growth overlap with M/S"),
        ifelse(cycle!="fix", "with cell cycle variance", ""))
    genes <- getGeneNames(pathway, overlap)
    geneColors <- ifelse(rownames(data) %in% genes, "black", "white")

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

plotTangPathways <- function(simulated, overlap="none", cycle="fix")
{
    if (simulated == 'simulated')
    {
        if (cycle == "fix")
        {
            mat <- switch(overlap, none=tang_none_fix,
                half=tang_half_fix, full=tang_full_fix)
        }
        else
        {
            mat <- switch(overlap, none=tang_none_var,
                half=tang_half_var, full=tang_full_var)
        }
        mat <- mat[,c(2:8,10:16)]
    }
    else
    {
        mat <- tangExpressionData[,c(14,6,8,10,12,2,4,13,5,7,9,11,1,3)]
        mat <- mat[rownames(mat) %in% allGenes,]
    }
    print(colnames(mat))
    plotTangHeatmap(mat, simulated, 'growth', overlap, cycle)
    plotTangHeatmap(mat, simulated, 'mitosis', overlap, cycle)
    plotTangHeatmap(mat, simulated, 'sphase', overlap, cycle)
    plotTangHeatmap(mat, simulated, 'contact_inhibition', overlap, cycle)
}

plotKagoharaPathways <- function(simulated, cellLine, overlap="none", cycle="fix")
{
    if (simulated == 'simulated')
    {
        ndx <- switch(cellLine, c(2:7,16:21,30:35,10:14,24:28,38:42),
            scc1=c(2:7,10:14), scc6=c(16:21,24:28), scc25=c(30:35,38:42))
        if (cycle == "fix")
        {
            mat <- switch(overlap, none=kagohara_none_fix,
                half=kagohara_half_fix, full=kagohara_full_fix)
        }
        else
        {
            mat <- switch(overlap, none=kagohara_none_var,
                half=kagohara_half_var, full=kagohara_full_var)
        }
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
    plotKagoharaHeatmap(mat, simulated, 'growth', overlap, cellLine, cycle)
    plotKagoharaHeatmap(mat, simulated, 'mitosis', overlap, cellLine, cycle)
    plotKagoharaHeatmap(mat, simulated, 'sphase', overlap, cellLine, cycle)
    plotKagoharaHeatmap(mat, simulated, 'contact_inhibition', overlap, cellLine, cycle)
}

plotKagoharaCellLine <- function(cellLine)
{
    plotKagoharaPathways("simulated", cellLine, "none", "fix")
    plotKagoharaPathways("simulated", cellLine, "none", "var")
    plotKagoharaPathways("simulated", cellLine, "half", "fix")
    plotKagoharaPathways("simulated", cellLine, "half", "var")
    plotKagoharaPathways("simulated", cellLine, "full", "fix")
    plotKagoharaPathways("simulated", cellLine, "full", "var")
    plotKagoharaPathways("real", cellLine)
}

# create heatmaps

plotKagoharaCellLine("all")
plotKagoharaCellLine("scc1")
plotKagoharaCellLine("scc6")
plotKagoharaCellLine("scc25")

plotTangPathways("simulated", "none", "fix")
plotTangPathways("simulated", "none", "var")
plotTangPathways("simulated", "half", "fix")
plotTangPathways("simulated", "half", "var")
plotTangPathways("simulated", "full", "fix")
plotTangPathways("simulated", "full", "var")
plotTangPathways("real")

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


