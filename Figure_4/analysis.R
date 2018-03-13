library(ggplot2)
library(gplots)
library(CancerInSilico)
load('Figure_4_cleaned.RData')
# STC.D PBS/CTX @ 0,24,48,72,96,120 hours
load("../../DataFromKagohara/STCCountsCG.Rda")
# tangExpressionData @ 0,4,8,24,48,72,96,120,144
load("../../DataFromTang/TangExpressionData.Rda")
data(SamplePathways)

# wrapper function for heatmap
internalHeatmap <- function(data, times, rowColors, colColors, title)
{
    heatmap.2(data, 
        col=greenred, scale='row',
        trace='none', hclust=function(x) hclust(x,method='complete'),
        distfun=function(x) as.dist((1-cor(t(x)))/2), 
        Colv=FALSE, dendrogram='row',
        ColSideColors = colColors,
        RowSideColors = rowColors,
        labCol = times,
        main=title)
}

getGeneNames <- function(genes, overlap)
{
    g <- pwyGrowth@genes
    m <- pwyMitosis@genes
    s <- pwySPhase@genes
    ci <- pwyContactInhibition@genes
    
    if (overlap == 'half')
    {
        m <- c(m, sample(g, length(g)/2, replace=FALSE))
        s <- c(s, setdiff(g, m))
    }
    else if (overlap == 'full')
    {
        m <- c(m, g)
        s <- c(s, g)
    }
    switch(genes, growth=g, mitosis=m, sphase=s, contact_inhibition=ci)
}

plotHeatmap <- function(data, overlap, genes, simulated)
{
    allGenes <- c(pwyGrowth@genes, pwyMitosis@genes, pwySPhase@genes, 
        pwyContactInhibition@genes)

    if (data == 'tang')
    {
        ndx <- c(1:7,9:15)
        times <- rep(seq(0,144,24),2)
        drugColors <- rep(c('black','red'), each=7)

        if (simulated)
        {
            mat <- switch(overlap, none=tang_none[,ndx], half=tang_half[,ndx],
                full=tang_full[,ndx])
        }
        else
        {
            colOrder <- c(14,6,8,10,12,2,4,13,5,7,9,11,1,3)
            tangExpressionData <- tangExpressionData[,colOrder]
            tangExpressionData <- tangExpressionData[rownames(tangExpressionData) %in% allGenes,]
            mat <- tangExpressionData
        }

    }
    else if (data == 'kagohara')
    {
        ndx <- c(1:6,13:18,25:30,8:12,20:24,32:36)
        mat <- switch(overlap, none=kagohara_none[,ndx],
            half=kagohara_half[,ndx], full=kagohara_full[,ndx])   
    }
    
    geneNames <- getGenes(genes, overlap)

    pdf(paste(data, '_', ifelse(simulated, 'sim', 'real'), '_', genes,
        '_', overlap, '.pdf', sep=''))
    internalHeatmap(mat, times, ifelse(rownames(mat) %in% geneNames,
        'black', 'white'), drugColors, paste(data,
        ifelse(simulated, 'simulated', 'real'), 'data - ', genes, 
        switch(overlap, none="", half=", half growth overlap with M/S",
            full=", full growth overlap with M/S")))
    dev.off()
}


#### TANG DATA ####

tang_ndx <- c(1:7,9:15)

plotHeatmap(tang_none[,tang_ndx], 'tang', 'none', 'growth', TRUE)
plotHeatmap('tang', 'none', 'mitosis', TRUE)
plotHeatmap('tang', 'none', 'sphase', TRUE)
plotHeatmap('tang', 'none', 'contact_inhibition', TRUE)

plotHeatmap('tang', 'half', 'growth', TRUE)
plotHeatmap('tang', 'half', 'mitosis', TRUE)
plotHeatmap('tang', 'half', 'sphase', TRUE)
plotHeatmap('tang', 'half', 'contact_inhibition', TRUE)

plotHeatmap('tang', 'full', 'growth', TRUE)
plotHeatmap('tang', 'full', 'mitosis', TRUE)
plotHeatmap('tang', 'full', 'sphase', TRUE)
plotHeatmap('tang', 'full', 'contact_inhibition', TRUE)

plotHeatmap('tang', 'none', 'growth', FALSE)
plotHeatmap('tang', 'none', 'mitosis', FALSE)
plotHeatmap('tang', 'none', 'sphase', FALSE)
plotHeatmap('tang', 'none', 'contact_inhibition', FALSE)

#### KAGOHARA DATA ####

plotHeatmap('kagohara', 'none', 'growth', TRUE)
plotHeatmap('kagohara', 'none', 'mitosis', TRUE)
plotHeatmap('kagohara', 'none', 'sphase', TRUE)
plotHeatmap('kagohara', 'none', 'contact_inhibition', TRUE)

plotHeatmap('kagohara', 'half', 'growth', TRUE)
plotHeatmap('kagohara', 'half', 'mitosis', TRUE)
plotHeatmap('kagohara', 'half', 'sphase', TRUE)
plotHeatmap('kagohara', 'half', 'contact_inhibition', TRUE)

plotHeatmap('kagohara', 'full', 'growth', TRUE)
plotHeatmap('kagohara', 'full', 'mitosis', TRUE)
plotHeatmap('kagohara', 'full', 'sphase', TRUE)
plotHeatmap('kagohara', 'full', 'contact_inhibition', TRUE)

plotHeatmap('kagohara', 'none', 'growth', FALSE)
plotHeatmap('kagohara', 'none', 'mitosis', FALSE)
plotHeatmap('kagohara', 'none', 'sphase', FALSE)
plotHeatmap('kagohara', 'none', 'contact_inhibition', FALSE)

kagoharaSimNdx <- c(1:6,13:18,25:30,8:12,20:24,32:36)
kagohara_none <- kagohara_none[,kagohara_ndx]
kagohara_half <- kagohara_half[,kagohara_ndx]
kagohara_full <- kagohara_full[,kagohara_ndx]

single_kagohara_times <- c(seq(0,120,24), seq(24,120,24))
all_kagohara_times <- c(rep(seq(0,120,24),3), rep(seq(24,120,24),3))

single_kagohara_colors <- c(rep('black',6), rep('red',5))
all_kagohara_colors <- c(rep('black',18), rep('red',15))

scc1_ndx <- c(1:6,8:12)
scc6_ndx <- c(13:18,20:24)
scc25_ndx <- c(25:30,32:36)


pdf('kagohara_sim_M_none.pdf', width=11, height=8.5)
plotHeatmap(kagohara_none[,kagoharaSimNdx], c(rep(seq(0,120,24),3), rep(seq(24,120,24),3)),
            ifelse(simGenes %in% mitosisGenes, 'black', 'white'),
            c(rep('black',18), rep('red',15)),
            "kagohara simulated data - mitosis genes")
dev.off()

pdf('kagohara_sim_S_none.pdf', width=11, height=8.5)
plotHeatmap(kagohara_none[,kagoharaSimNdx], c(rep(seq(0,120,24),3), rep(seq(24,120,24),3)),
            ifelse(simGenes %in% SPhaseGenes, 'black', 'white'),
            c(rep('black',18), rep('red',15)),
            "kagohara simulated data - S Phase genes")
dev.off()

pdf('kagohara_sim_CI_none.pdf', width=11, height=8.5)
plotHeatmap(kagohara_none[,kagoharaSimNdx], c(rep(seq(0,120,24),3), rep(seq(24,120,24),3)),
            ifelse(simGenes %in% contactInhibitionGenes, 'black', 'white'),
            c(rep('black',18), rep('red',15)),
            "kagohara simulated data - contact inhibition genes")
dev.off()

# 100 overlap

#overlap <- c(pwyGrowth@genes, pwyMitosis@genes, pwySPhase@genes)
#growthGenes <- overlap
mitosisGenes <- c(pwyMitosis@genes, pwyGrowth@genes)
SPhaseGenes <- c(pwySPhase@genes, pwyGrowth@genes)
contactInhibitionGenes <- pwyContactInhibition@genes

pdf('kagohara_sim_G_full.pdf', width=11, height=8.5)
plotHeatmap(kagohara_full[,kagoharaSimNdx], c(rep(seq(0,120,24),3), rep(seq(24,120,24),3)),
            ifelse(simGenes %in% growthGenes, 'black', 'white'),
            c(rep('black',18), rep('red',15)),
            "kagohara simulated data - growth genes, growth/M/S overlap")
dev.off()

pdf('kagohara_sim_G_full_scc25.pdf', width=11, height=8.5)
plotHeatmap(kagohara_full[,scc25_ndx], c(seq(0,120,24), seq(24,120,24)),
            ifelse(simGenes %in% growthGenes, 'black', 'white'),
            c(rep('black',6), rep('red',5)),
            "kagohara simulated data - growth genes, growth/M/S overlap")
dev.off()

pdf('kagohara_sim_M_full.pdf', width=11, height=8.5)
plotHeatmap(kagohara_full[,kagoharaSimNdx], c(rep(seq(0,120,24),3), rep(seq(24,120,24),3)),
            ifelse(simGenes %in% mitosisGenes, 'black', 'white'),
            c(rep('black',18), rep('red',15)),
            "kagohara simulated data - mitosis genes, growth/M/S overlap")
dev.off()

pdf('kagohara_sim_S_full.pdf', width=11, height=8.5)
plotHeatmap(kagohara_full[,kagoharaSimNdx], c(rep(seq(0,120,24),3), rep(seq(24,120,24),3)),
            ifelse(simGenes %in% SPhaseGenes, 'black', 'white'),
            c(rep('black',18), rep('red',15)),
            "kagohara simulated data - S Phase genes, growth/M/S overlap")
dev.off()

pdf('kagohara_sim_CI_full.pdf', width=11, height=8.5)
plotHeatmap(kagohara_full[,kagoharaSimNdx], c(rep(seq(0,120,24),3), rep(seq(24,120,24),3)),
            ifelse(simGenes %in% contactInhibitionGenes, 'black', 'white'),
            c(rep('black',18), rep('red',15)),
            "kagohara simulated data - contact inhibition genes, growth/M/S overlap")
dev.off()

# real data

overlap <- c(pwyGrowth@genes, pwyMitosis@genes, pwySPhase@genes)
growthGenes <- overlap
mitosisGenes <- overlap
SPhaseGenes <- overlap
contactInhibitionGenes <- pwyContactInhibition@genes
allGenes <- c(overlap, pwyContactInhibition@genes)

#colOrder <- c(1:6,12:17,23:28,7:11,18:22,29:33)
#STC.D <- STC.D[,colOrder]
STC.D <- STC.D[rownames(STC.D) %in% allGenes,]
print(colnames(STC.D))
realGenes <- rownames(STC.D)

pdf('kagohara_real_G.pdf', width=11, height=8.5)
plotHeatmap(STC.D, c(rep(seq(0,120,24),3), rep(seq(24,120,24),3)),
            ifelse(realGenes %in% growthGenes, 'black', 'white'),
            c(rep('black',18), rep('red',15)),
            "kagohara data - growth genes")
dev.off()

pdf('kagohara_real_G_scc25.pdf', width=11, height=8.5)
plotHeatmap(STC.D[,23:33], c(seq(0,120,24), seq(24,120,24)),
            ifelse(realGenes %in% growthGenes, 'black', 'white'),
            c(rep('black',6), rep('red',5)),
            "kagohara data - growth genes")
dev.off()

pdf('kagohara_real_M.pdf', width=11, height=8.5)
plotHeatmap(STC.D, c(rep(seq(0,120,24),3), rep(seq(24,120,24),3)),
            ifelse(realGenes %in% mitosisGenes, 'black', 'white'),
            c(rep('black',18), rep('red',15)),
            "kagohara data - mitosis genes")
dev.off()

pdf('kagohara_real_S.pdf', width=11, height=8.5)
plotHeatmap(STC.D, c(rep(seq(0,120,24),3), rep(seq(24,120,24),3)),
            ifelse(realGenes %in% SPhaseGenes, 'black', 'white'),
            c(rep('black',18), rep('red',15)),
            "kagohara data - S Phase genes")
dev.off()

pdf('kagohara_real_CI.pdf', width=11, height=8.5)
plotHeatmap(STC.D, c(rep(seq(0,120,24),3), rep(seq(24,120,24),3)),
            ifelse(realGenes %in% contactInhibitionGenes, 'black', 'white'),
            c(rep('black',18), rep('red',15)),
            "kagohara data - contact inhibition genes")
dev.off()


## plot pathway activity

fig4c <- ggplot(tangPwyActivity, aes(x=hour)) +
    geom_line(aes(y=activationGrowth_100ug, color='Growth Rate', linetype='100ug')) +
    geom_line(aes(y=mitosis_100ug,          color='Mitosis', linetype='100ug')) +
    geom_line(aes(y=contactInhibition_100ug,color='Contact Inhibition', linetype='100ug')) +
    geom_line(aes(y=activationGrowth_pbs,   color='Growth Rate', linetype='PBS')) +
    geom_line(aes(y=mitosis_pbs,            color='Mitosis', linetype='PBS')) +
    geom_line(aes(y=contactInhibition_pbs,  color='Contact Inhibition', linetype='PBS')) +
    labs(title = "Pathway Activity", caption = "Figure 4c", x = "Hour",
        y = "Percent Pathway Activity", color = "Pathway",
        linetype = "Dosage") +
    scale_linetype_manual(values = c('dotted', 'solid'))
ggsave(filename='fig4c.pdf', plot=fig4c)