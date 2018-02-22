
load('fig6_expression.RData')

table(fig6_single_cell_phases)
table(fig6_single_cell_types)

#cells<-sapply(colnames(fig6_two_cell_types_single_cell_RNAseq),function(x) strsplit(x,"_")[[1]][1])
time<-sapply(colnames(fig6_two_cell_types_single_cell_RNAseq),function(x) strsplit(x,"_")[[1]][2])
time<-as.numeric(gsub("t","",time))

phase<-rep("blue",length(fig6_single_cell_phases))
phase[which(fig6_single_cell_phases=="M")]<-"cyan"
phase[which(fig6_single_cell_phases=="S")]<-"orange"

type<-rep("red",length(fig6_single_cell_phases))
type[which(fig6_single_cell_types==2)]<-"purple"

table(time)
table(cells)

library(Rtsne) # Load package
set.seed(42) # Sets seed for reproducibility
tsne_out <- Rtsne(t(fig6_two_cell_types_single_cell_RNAseq), dims = 3, initial_dims = 3, perplexity = 30,
       theta = 0.5, pca = TRUE, max_iter = 1000,verbose = TRUE)

library("viridis")
myColorRamp <- function(values,palette=viridis(255)) {
  if(min(values)<0){
    values<-values+abs(min(values))
  }
  v <- (values - min(values))/diff(range(values))
  #v <- min(values/diff(range(values)),1.0)
  x <- colorRamp(palette)(v)
  print(x)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

COL<-myColorRamp(time)

library(rgl)

plot3d(tsne_out$Y,col=phase)
rgl.postscript("fig6_tsne_coloredbyPHASE.pdf", fmt = "pdf", drawText = TRUE )


plot3d(tsne_out$Y,col=type)
rgl.postscript("fig6_tsne_coloredbyTYPE.pdf", fmt = "pdf", drawText = TRUE )


COL<-myColorRamp(time)
plot3d(tsne_out$Y,col=COL)
rgl.postscript("fig6_tsne_coloredbyTIME.pdf", fmt = "pdf", drawText = TRUE )

#make movie 
rot <- spin3d(axis=c( 0, 0, 1 ))
movie3d(rot, duration= 12, type="gif",dir="~/Downloads/",movie="fig6_tsne_wPCA")