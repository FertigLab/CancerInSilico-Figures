### R libraries
library('limma')
library('simpleaffy')
library('ClassDiscovery')

### final report with normalized signal for each probe
D <- read.table('../Data/forPartekFinalReport_AVGSignal.txt',sep="\t",stringsAsFactors=F,
                row.names=1)

### header file describing each of the columns in the final report file
hdr <- read.table('../Data/forPartekFinalReport_AVGSignal.txt.fmt',
                  quote="",sep=",",stringsAsFactors=FALSE)

# extracting the column names from the formatted header file
cn <- c(grep('cl',hdr[,1],value=T),hdr[(grep('colgroup',hdr[,1])+1):nrow(hdr),1])

# remove the first column names, as it was set as the row.names of the data
colnames(D) <- cn[2:length(cn)]

### extract objects to create a LumiBatch object!
lumi.D <- list()

# expression from true beads
lumi.D$exprs = t(data.matrix(D[,grep('cl',colnames(D),invert=T,value=T)]))

# data from control probes
lumi.D$controlData <- t(data.matrix(D[,setdiff(colnames(D),
      row.names(lumi.D$exprs))[14:length(setdiff(colnames(D),
                                                 row.names(lumi.D$exprs)))]]))
row.names(lumi.D$controlData) <- 
  sapply(strsplit(row.names(lumi.D$controlData),split='[{}]'),function(x){x[[2]]})

# history information
lumi.D$history <- data.frame(info=as.character(hdr[7:8,1]),stringsAsFactors=F)

### loess normalization
exprsNorm <- normalizeBetweenArrays(lumi.D$exprs,method='cyclicloess')

# store normalized data
lumi.D$exprs <- exprsNorm
lumi.D$history <- rbind(lumi.D$history,'loess normalization with limma')

### information about genes
gAnnot <- read.table('../Data/forPartekFinalReport_IlluminaHT-12_v4.annotation.txt',
                     header=T,sep="\t",quote="",comment.char="",
                     stringsAsFactors=F,row.names=1)
gAnnot <- gAnnot[row.names(lumi.D$exprs),]

lumi.D$featureData <- gAnnot

### information about samples
lumi.D$pheno <- read.table('../Data/SamplesAnnot.csv',header=T,sep=",",
                    stringsAsFactors=F)
row.names(lumi.D$pheno) <- paste(lumi.D$pheno$Array,
                                 lumi.D$pheno$Santrix.Positions,sep="_")

# distinguish treatment conditions from sample names
lumi.D$pheno$Treatment <- 'Control'
lumi.D$pheno$Treatment[grep('M',lumi.D$pheno$Samples.IDS)] <- 'Cetuximab'

# get time and replicate information
lumi.D$pheno[,c('TimePt','Replicate')] <- 
  t(sapply(strsplit(lumi.D$pheno$Samples.IDS,split='[CM]'),
           function(x){c(x[[1]],x[[2]])}))

# set experiment type
lumi.D$pheno$ExpType <- substr(lumi.D$pheno$Samples.IDS,1,nchar(lumi.D$pheno$Samples.IDS)-1)

### ouptut formatted data
save(lumi.D,file='../Data/Tang_Lumi.Rda')
