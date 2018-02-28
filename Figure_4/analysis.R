#library(ggplot2)
#library(gplots)
#library(methods)
load('Figure_4_cleaned.RData')
#data(SamplePathways)


# plot heatmaps

png(file='fig4d.png')
heatmap.2(ge_100ug$expression, col=greenred, scale='row', trace='none',
    hclust=function(x) hclust(x,method='complete'),
    distfun=function(x) as.dist((1-cor(t(x)))/2), Colv=F, dendrogram='row')

png(file='fig4d_growth.png')
heatmap.2(ge_100ug$expression[pwyGrowth@genes,], col=greenred, scale='row',
    trace='none', hclust=function(x) hclust(x,method='complete'),
    distfun=function(x) as.dist((1-cor(t(x)))/2), Colv=F, dendrogram='row')

png(file='fig4d_mitosis.png')
heatmap.2(ge_100ug$expression[pwyMitosis@genes,], col=greenred, scale='row',
    trace='none', hclust=function(x) hclust(x,method='complete'),
    distfun=function(x) as.dist((1-cor(t(x)))/2), Colv=F, dendrogram='row')

png(file='fig4d_sphase.png')
heatmap.2(ge_100ug$expression[pwySPhase@genes,], col=greenred, scale='row',
    trace='none', hclust=function(x) hclust(x,method='complete'),
    distfun=function(x) as.dist((1-cor(t(x)))/2), Colv=F, dendrogram='row')

png(file='fig4d_contactInhibition.png')
heatmap.2(ge_100ug$expression[pwyContactInhibition@genes,], col=greenred, scale='row',
    trace='none', hclust=function(x) hclust(x,method='complete'),
    distfun=function(x) as.dist((1-cor(t(x)))/2), Colv=F, dendrogram='row')

geneNames <- unique(c(pwyContactInhibition@genes, pwyMitosis@genes, pwySPhase@genes,
    pwyGrowth@genes))
png(file='fig4d_realData.png')
heatmap.2(referenceGeneExpression[geneNames,], col=greenred, scale='row',
    trace='none', hclust=function(x) hclust(x,method='complete'),
    distfun=function(x) as.dist((1-cor(t(x)))/2), Colv=F, dendrogram='row')

fig4c <- ggplot(pwyActivity, aes(x=hour)) +
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

fig4c_growth <- ggplot(pwyActivity, aes(x=hour)) +
    geom_line(aes(y=activationGrowth_pbs,   color='pbs')) +
    geom_line(aes(y=activationGrowth_10ug,  color='10ug')) +
    geom_line(aes(y=activationGrowth_100ug, color='100ug')) + 
    labs(title = "Pathway Activity - Growth Rate", caption = "Figure 4c",
        x = "Hour", y = "Percent Pathway Activity", color = "Dosage") 

fig4c_contactInhibition <- ggplot(pwyActivity, aes(x=hour)) +
    geom_line(aes(y=contactInhibition_pbs,  color='pbs')) +
    geom_line(aes(y=contactInhibition_10ug, color='10ug')) +
    geom_line(aes(y=contactInhibition_100ug,color='100ug')) + 
    labs(title = "Pathway Activity - Contact Inhibition", caption = "Figure 4c",
        x = "Hour", y = "Percent Pathway Activity", color = "Dosage") 

fig4c_mitosis <- ggplot(pwyActivity, aes(x=hour)) +
    geom_line(aes(y=mitosis_pbs,            color='pbs')) +
    geom_line(aes(y=mitosis_10ug,           color='10ug')) +
    geom_line(aes(y=mitosis_100ug,          color='100ug')) + 
    labs(title = "Pathway Activity - Mitosis", caption = "Figure 4c",
        x = "Hour", y = "Percent Pathway Activity", color = "Dosage") 

fig4c_sphase <- ggplot(pwyActivity, aes(x=hour)) +
    geom_line(aes(y=GtoS_pbs,               color='pbs')) +
    geom_line(aes(y=GtoS_10ug,              color='10ug')) +
    geom_line(aes(y=GtoS_100ug,             color='100ug')) +
    labs(title = "Pathway Activity - S Phase", caption = "Figure 4c",
        x = "Hour", y = "Percent Pathway Activity", color = "Dosage") 

ggsave(filename='fig4c.png', plot=fig4c)
ggsave(filename='fig4c_growth.png', plot=fig4c_growth)
ggsave(filename='fig4c_contactInhibition.png', plot=fig4c_contactInhibition)
ggsave(filename='fig4c_mitosis.png', plot=fig4c_mitosis)
ggsave(filename='fig4c_sphase.png', plot=fig4c_sphase)
