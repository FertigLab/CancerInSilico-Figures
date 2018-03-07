load('Figure_6_cleaned.RData') # pwyActivity, ge, ge_bulk, cellPhase, cellType

table(cellPhase)
table(cellType)

#cells<-sapply(colnames(ge),function(x) strsplit(x,"_")[[1]][1])
ge <- ge$expression
time <- sapply(colnames(ge),function(x) strsplit(x,"_")[[1]][2])
time <- as.numeric(gsub("t","",time))

phase<-rep("blue",length(cellPhase))
phase[which(cellPhase=="M")]<-"cyan"
phase[which(cellPhase=="S")]<-"orange"

type<-rep("red",length(cellPhase))
type[which(cellType==2)]<-"purple"

table(time)
#table(cells)

library(Rtsne) # Load package
set.seed(42) # Sets seed for reproducibility
tsne_out <- Rtsne(t(ge), dims = 3, initial_dims = 3, perplexity = 30,
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
#rot <- spin3d(axis=c( 0, 0, 1 ))
#movie3d(rot, duration= 12, type="gif",movie="fig6_tsne_wPCA")