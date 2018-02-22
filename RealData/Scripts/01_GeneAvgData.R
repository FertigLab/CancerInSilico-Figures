# load in the gene expression data
load('../Data/Tang_Lumi.Rda')

# get the genes associated with each probe in TF
genes2Probes <- tapply(row.names(lumi.D$featureData),
                       lumi.D$featureData$SYMBOL,paste) 
genes2Probes <- genes2Probes[setdiff(names(genes2Probes),"")]

# get gene level averages
D.gene <- matrix(0.,nrow=length(genes2Probes),ncol=nrow(lumi.D$pheno),
                 dimnames=list(names(genes2Probes),
                               row.names(lumi.D$pheno)))
D.gene.1prb <- lumi.D$exprs[unlist(genes2Probes[sapply(genes2Probes,length)==1]),
                            colnames(D.gene)]
row.names(D.gene.1prb) <- names(unlist(genes2Probes[sapply(genes2Probes,length)==1]))

D.gene.nprb <- t(sapply(genes2Probes[sapply(genes2Probes,length)>1],
                        function(x){apply(lumi.D$exprs[x,colnames(D.gene)],
                                          2,mean)}))
D.gene[row.names(D.gene.1prb),colnames(D.gene.1prb)] <- D.gene.1prb
D.gene[row.names(D.gene.nprb),colnames(D.gene.nprb)] <- D.gene.nprb

save(D.gene, file='../Data/Tang_GeneAvg.Rda')
