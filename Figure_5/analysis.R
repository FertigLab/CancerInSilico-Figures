library('CancerInSilico')
library('ggplot2')
library('RColorBrewer')  
library('reshape2') 
library(methods)
load("Figure_5_cleaned.RData")

error.bar <- function(x, y, upper, lower=upper, length=0.1,...)
{
    if (length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

singleTypeRun <- fig5Data[sapply(fig5Data, function(mod) mod$numTypes==1)]
twoTypesRun <- fig5Data[sapply(fig5Data, function(mod) mod$numTypes==2)]
#fourTypesRun <- fig5Data[sapply(fig5Data, function(mod) mod$numTypes==4)]
#sixTypesRun <- fig5Data[sapply(fig5Data, function(mod) mod$numTypes==6)]

# single type

# two types
lowVarianceRuns <- twoTypesRun[sapply(twoTypesRun, function(mod)
    mod$typeFreq[1] > 0.4 & mod$typeFreq[1] < 0.6)]
medVarianceRuns <- twoTypesRun[sapply(twoTypesRun, function(mod)
    mod$typeFreq[1] > 0.2 & mod$typeFreq[1] < 0.8)]
highVarianceRuns <- twoTypesRun[sapply(twoTypesRun, function(mod)
    mod$typeFreq[1] > 0.0 & mod$typeFreq[1] < 1.0)]

# totals
png('fig5_low_total.png')
plot(NULL, xlim=c(0,168), ylim=c(0,1000), xlab="time (hrs)",
    ylab="number of cells", main="Population Size Over Time")
eat <- sapply(lowVarianceRuns, function(mod) lines(1:length(mod$numCells), mod$numCells))

png('fig5_med_total.png')
plot(NULL, xlim=c(0,168), ylim=c(0,1000), xlab="time (hrs)",
    ylab="number of cells", main="Population Size Over Time")
eat <- sapply(medVarianceRuns, function(mod) lines(1:length(mod$numCells), mod$numCells))

png('fig5_high_total.png')
plot(NULL, xlim=c(0,168), ylim=c(0,1000), xlab="time (hrs)",
    ylab="number of cells", main="Population Size Over Time")
eat <- sapply(highVarianceRuns, function(mod) lines(1:length(mod$numCells), mod$numCells))

# init freq
png('fig5_low_freq.png')
freqMat <- unname(sapply(lowVarianceRuns, function(mod) mod$typeFreq))
meanVec <- apply(freqMat, 1, mean)
sdVec <- apply(freqMat, 1, sd)
barx <- barplot(meanVec, ylim=c(0,1.0), axis.lty=1, xlab="Cell Type",
    ylab="Initial Frequency", main="Variance of Initial Cell Type Proportion")
error.bar(barx, meanVec, 1.96 * sdVec / sqrt(length(sdVec)))

png('fig5_med_freq.png')
freqMat <- unname(sapply(medVarianceRuns, function(mod) mod$typeFreq))
meanVec <- apply(freqMat, 1, mean)
sdVec <- apply(freqMat, 1, sd)
barx <- barplot(meanVec, ylim=c(0,1.0), axis.lty=1, xlab="Cell Type",
    ylab="Initial Frequency", main="Variance of Initial Cell Type Proportion")
error.bar(barx, meanVec, 1.96 * sdVec / sqrt(length(sdVec)))

png('fig5_high_freq.png')
freqMat <- unname(sapply(highVarianceRuns, function(mod) mod$typeFreq))
meanVec <- apply(freqMat, 1, mean)
sdVec <- apply(freqMat, 1, sd)
barx <- barplot(meanVec, ylim=c(0,1.0), axis.lty=1, xlab="Cell Type",
    ylab="Initial Frequency", main="Variance of Initial Cell Type Proportion")
error.bar(barx, meanVec, 1.96 * sdVec / sqrt(length(sdVec)))

# distribution
#png('fig5_low_dist.png')
#plot(density(sapply(lowVarianceRuns, function(mod) mod$numCells[48])), 
#    main="distribution of cell population at 48 hours - 2 types",
#    xlab="number of cells", ylab="density")
#
#png('fig5_med_dist.png')
#plot(density(sapply(medVarianceRuns, function(mod) mod$numCells[48])), 
#    main="distribution of cell population at 48 hours - 2 types",
#    xlab="number of cells", ylab="density")
#
#png('fig5_high_dist.png')
#plot(density(sapply(highVarianceRuns, function(mod) mod$numCells[48])), 
#    main="distribution of cell population at 48 hours - 2 types",
#    xlab="number of cells", ylab="density")

#png('fig5_2_total.png')
#plot(NULL, xlim=c(0,168), ylim=c(0,1000), xlab="time (hrs)",
#    ylab="number of cells", main="two cell types")
#eat <- sapply(twoTypesRun, function(mod) lines(1:length(mod$numCells), mod$numCells))
#
#png('fig5_2_freq.png')
#freqMat <- unname(sapply(twoTypesRun, function(mod) mod$typeFreq))
#meanVec <- apply(freqMat, 1, mean)
#sdVec <- apply(freqMat, 1, sd)
#barx <- barplot(meanVec, ylim=c(0,0.75), axis.lty=1, xlab="Cell Type",
#    ylab="Initial Frequency", main="Two cell types - initial proportion")
#error.bar(barx, meanVec, 1.96 * sdVec / sqrt(length(sdVec)))
#
#png('fig5_2_stddev.png')
#plot(density(sapply(twoTypesRun, function(mod) mod$numCells[48])), 
#    main="distribution of cell population at 48 hours - 2 types",
#    xlab="number of cells", ylab="density")

# 4 types
#png('fig5_4_total.png')
#plot(NULL, xlim=c(0,168), ylim=c(0,1000), xlab="time (hrs)",
#    ylab="number of cells", main="four cell types")
#eat <- sapply(fourTypesRun, function(mod) lines(1:length(mod$numCells), mod$numCells))
#
#png('fig5_4_freq.png')
#freqMat <- unname(sapply(fourTypesRun, function(mod) mod$typeFreq))
#meanVec <- apply(freqMat, 1, mean)
#sdVec <- apply(freqMat, 1, sd)
#barx <- barplot(meanVec, ylim=c(0,0.75), axis.lty=1, xlab="Cell Type",
#    ylab="Initial Frequency", main="Four cell types - initial proportion")
#error.bar(barx, meanVec, 1.96 * sdVec / sqrt(length(sdVec)))
#
#png('fig5_4_stddev.png')
#plot(density(sapply(fourTypesRun, function(mod) mod$numCells[48])), 
#    main="distribution of cell population at 48 hours - 4 types",
#    xlab="number of cells", ylab="density")
#
### 6 types
#png('fig5_6_total.png')
#plot(NULL, xlim=c(0,168), ylim=c(0,1000), xlab="time (hrs)",
#    ylab="number of cells", main="six cell types")
#eat <- sapply(sixTypesRun, function(mod) lines(1:length(mod$numCells), mod$numCells))
#
#png('fig5_6_freq.png')
#freqMat <- unname(sapply(sixTypesRun, function(mod) mod$typeFreq))
#meanVec <- apply(freqMat, 1, mean)
#sdVec <- apply(freqMat, 1, sd)
#barx <- barplot(meanVec, ylim=c(0,0.75), axis.lty=1, xlab="Cell Type",
#    ylab="Initial Frequency", main="Six cell types - initial proportion")
#error.bar(barx, meanVec, 1.96 * sdVec / sqrt(length(sdVec)))
#
#png('fig5_6_stddev.png')
#plot(density(sapply(sixTypesRun, function(mod) mod$numCells[48])), 
#    main="distribution of cell population at 48 hours - 6 types",
#    xlab="number of cells", ylab="density")


#png('fig5a2.png')
#plot(NULL, xlim=c(0,168), ylim=c(0,1000))
#eat <- sapply(threeTypesRun, function(mod) lines(1:length(mod$numCells), mod$numCells))
#print(sd(sapply(threeTypesRun, function(mod) mod$numCells[50])))

#png('fig5a3.png')
#plot(NULL, xlim=c(0,168), ylim=c(0,1000))
#eat <- sapply(fiveTypesRun, function(mod) lines(1:length(mod$numCells), mod$numCells))
#print(sd(sapply(fiveTypesRun, function(mod) mod$numCells[50])))

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

