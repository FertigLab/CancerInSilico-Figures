library(Rtsne)
library(viridis)
library(rgl)
load("fig6.RData") # ge, cellType, cellPhase
set.seed(123)

times <- unname(sapply(colnames(ge), function(x) strsplit(x,"_")[[1]][2]))
times <- as.numeric(gsub("t", "", times))
if (min(times) != 0 | max(times) != 168) stop('invalid time range')

cells <- unname(sapply(colnames(ge), function(x) strsplit(x,"_")[[1]][1]))
cells <- as.numeric(gsub("c", "", cells))

if (length(times) != length(cells)) stop('dim mismatch')

type <- ifelse(cellType[cells] == 1, "red", "purple")

phase <- sapply(1:length(times), function(i) cellPhase[cells[i], times[i]+1])
phase[phase == "I"] <- "blue"
phase[phase == "S"] <- "orange"
phase[phase == "M"] <- "red"

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

open3d()
par3d(windowRect=c(200,200,1400,1400))

plot3d(tsne_out$Y, col=type, xlab="", ylab="", zlab="")
legend3d("topright", legend=c("Type A", "Type B"), pch=16, cex=1,
    col=c("red", "purple"), inset=c(0.04))
rgl.postscript("fig6_tsne_coloredbyTYPE.pdf", fmt="pdf")

plot3d(tsne_out$Y, col=phase, xlab="", ylab="", zlab="")
legend3d("topright", legend=c("Interphase", "SPhase", "Mitosis"), pch=16, cex=1,
    col=c("blue", "orange", "red"), inset=c(0.02))
rgl.postscript("fig6_tsne_coloredbyPHASE.pdf", fmt="pdf", drawText=TRUE)

plot3d(tsne_out$Y, col=myColorRamp(times), xlab="", ylab="", zlab="")
legend3d("topright", legend=c("0 hours", "168 hours"), pch=16, cex=1,
    col=c(min(myColorRamp(times)), max(myColorRamp(times))), inset=c(0.02))
rgl.postscript("fig6_tsne_coloredbyTIME.pdf", fmt="pdf")

#M <- par3d("userMatrix")
#par3d(userMatrix=rotate3d(M, angle=pi, x=1, y=0, z=0))

save(tsne_out, file="temp.RData")