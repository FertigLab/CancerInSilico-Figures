library('CancerInSilico')
library('ggplot2')
library('RColorBrewer')  
library('reshape2') 
library(methods)
load("Figure_5_cleaned.RData")

singleLow <- fig5Data[sapply(fig5Data,
    function(mod) mod$freqA==1 & mod$sdA < 0.2)]

singleHigh <- fig5Data[sapply(fig5Data,
    function(mod) mod$freqA==1 & mod$sdA > 8)]

doubleFar <- fig5Data[sapply(fig5Data,
    function(mod) mod$freqA!=1 & mod$meanA < 12)]

doubleNear <- fig5Data[sapply(fig5Data,
    function(mod) mod$freqA!=1 & mod$meanA > 20)]

png('fig5a1.png')
plot(NULL, xlim=c(0,168), ylim=c(0,1000))
sapply(singleLow, function(mod) lines(1:length(mod$numCells), mod$numCells))

png('fig5a2.png')
plot(NULL, xlim=c(0,168), ylim=c(0,1000))
sapply(singleHigh, function(mod) lines(1:length(mod$numCells), mod$numCells))

png('fig5a3.png')
plot(NULL, xlim=c(0,168), ylim=c(0,1000))
sapply(doubleFar, function(mod) lines(1:length(mod$numCells), mod$numCells))

png('fig5a4.png')
plot(NULL, xlim=c(0,168), ylim=c(0,1000))
sapply(doubleNear, function(mod) lines(1:length(mod$numCells), mod$numCells))

#fig <- ggplot(subset(fig5data, density == 0.05 & cycleLength %in% c(12,18,24)), aes(x=time)) + 
    #geom_point(aes(y=cellTypeBFreq)) 
#ggsave(filename='fig5a.png', plot=fig)

#mat <- matrix(nrow=length(cellTypeBInitFreq), ncol=length(cellTypeBCycleLength))
#for (data in fig5data)
#{
#	xind <- which(cellTypeBInitFreq == data[1])
#	yind <- which(cellTypeBCycleLength == data[2])
#	mat[xind, yind] <- data[3]
#}
#
#fpm <- mat
#fpm.melted <- melt(fpm)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')
#fig <- ggplot(fpm.melted, aes(x = Var1, y = Var2, fill = value)) +
#        geom_tile() +
#        scale_fill_gradientn(colours = hm.palette(100)) +
#        ylab('Cycle Length of Faster Growing Cell Type') +
#        xlab('Initial Proportion of Slower Growing Cell Type') +
#        ggtitle("With boundary") +
#        theme(plot.title = element_text(hjust = 0.5)) +
#        theme(axis.text=element_text(size=10))

