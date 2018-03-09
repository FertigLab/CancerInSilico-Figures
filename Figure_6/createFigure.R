library(Rtsne)
library(viridis)
library(rgl)
load("fig6.RData") # ge, cellType, cellPhase
set.seed(123)

time <- sapply(colnames(ge), function(x) strsplit(x,"_")[[1]][2])
time <- as.numeric(gsub("t", "", time))

phase <- rep("blue",length(cellPhase))
phase[which(cellPhase == "M")] <- "orange"
phase[which(cellPhase == "S")] <- "cyan"

type <- rep("red", length(cellPhase))
type[which(cellType == 2)] <- "purple"

tsne_out <- Rtsne(t(ge), dims=3, initial_dims=3, perplexity=30,
       theta=0.5, pca=TRUE, max_iter=1000, verbose=TRUE)

myColorRamp <- function(values, palette=viridis(255))
{
    if (min(values) < 0)
    {  
       values <- values + abs(min(values))
    }
    v <- (values - min(values)) / diff(range(values))
    x <- colorRamp(palette)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
}

par3d(windowRect=c(0,0,800,800))

plot3d(tsne_out$Y, col=phase, xlab="", ylab="", zlab="")
rgl.postscript("fig6_tsne_coloredbyPHASE.pdf", fmt="pdf", drawText=TRUE)

plot3d(tsne_out$Y, col=type, xlab="", ylab="", zlab="")
rgl.postscript("fig6_tsne_coloredbyTYPE.pdf", fmt="pdf", drawText=TRUE)

plot3d(tsne_out$Y, col=myColorRamp(time), xlab="", ylab="", zlab="")
rgl.postscript("fig6_tsne_coloredbyTIME.pdf", fmt="pdf", drawText=TRUE)

#make movie 
#rot <- spin3d(axis=c( 0, 0, 1 ))
#movie3d(rot, duration= 12, type="gif",movie="fig6_tsne_wPCA")