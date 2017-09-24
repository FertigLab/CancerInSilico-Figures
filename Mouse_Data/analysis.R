library(ggplot2)
load('~/FertigLab/CancerInSilico/Figures/Mouse_Data/PDXData_cleaned.RData')
pdxData <- na.omit(pdxData)

## get distribution

pbsDist <- list()
cetDist <- list()
days <- sort(unique(pdxData$day))
for (d in days)
{
    pbsData <- subset(pdxData, day > d - 3 & day < d + 3 & treatment == 'PBS')$tumorVol
    if (length(pbsData) > 30)
        pbsDist[[d]] <- pbsData

    cetData <- subset(pdxData, day > d - 3 & day < d + 3 & treatment == 'CET')$tumorVol
    if (length(cetData) > 30)
        cetDist[[d]] <- cetData
}

## get stats

days <- which(!sapply(pbsDist, is.null))
pbsStats <- data.frame(day=days, treatment=rep('PBS', length(days)), mean=sapply(pbsDist[days], median), sd=sapply(pbsDist[days], sd))

days <- which(!sapply(cetDist, is.null))
cetStats <- data.frame(day=days, treatment=rep('CET', length(days)), mean=sapply(cetDist[days], median), sd=sapply(cetDist[days], sd))

stats <- rbind(pbsStats, cetStats)
ggplot(stats, aes(x=day, color=treatment)) + geom_line(aes(y=mean)) + geom_line(aes(y=sd), linetype=3)

## plots

ggplot(subset(pdxData, model=='409_PDX_F2'), aes(x=day, color=treatment)) + 
    geom_line(aes(y=tumorVol, linetype=as.character(mouseID)))

ggplot(subset(pdxData, treatment=='PBS'), aes(x=day, color=model)) + 
    geom_line(aes(y=tumorVol, linetype=as.character(mouseID)))

ggplot(subset(pdxData, treatment=='CET'), aes(x=day, color=model)) + 
  geom_line(aes(y=tumorVol, linetype=as.character(mouseID)))

ggplot(pdxData, aes(x=day, color=model)) + 
    geom_line(aes(y=tumorVol, linetype=as.character(mouseID)))


