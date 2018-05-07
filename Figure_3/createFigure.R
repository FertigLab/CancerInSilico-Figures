library(ggplot2)
load('Figure_3_cleaned.RData') #fig3Data
load('tangData.RData') #tangData
load('KagoharaGrowthData.RData')

# useful function
l2Norm <- function(sim, real)
{
    sim <- sim * real[1] / sim[1]
    sqrt(sum((sim-real)^2))
}

#### Tang Data (3a)

# get real data
realData <- subset(tangData, dosage %in% c(0, 10, 100))
realData$Fit <- rep(0, nrow(realData))
colnames(realData)[colnames(realData) == "numCells"] <- "cellArea"

# fit pbs data
simNdx <- 24 + seq(0,120,24) + 1 # offset 24 for smoother gene expression
noDrugSims <- fig3Data[sapply(fig3Data, function(d) d$drugEffect==1)]
l2 <- sapply(noDrugSims, function(sim)
{
    simArea <- sim$cellArea[simNdx] 
    realArea <- realData[2:7,]$cellArea
    return(l2Norm(simArea, realArea))
})
pbsFit <- noDrugSims[[which(l2==min(l2))]]

# scale fit data to match initial time point of real data
realData[2:7,]$Fit <- pbsFit$cellArea[simNdx] * realData[2,]$cellArea /
    pbsFit$cellArea[simNdx][1]

print(c(pbsFit$cycleLength, pbsFit$initDensity))
tangFit <- list("pbs"=pbsFit)

# fit ctx data
drugSims <- fig3Data[sapply(fig3Data, function(d)
    d$initDensity==pbsFit$initDensity & d$cycleLength==pbsFit$cycleLength)]
for (i in c(8,15)) # day 1 index for ctx
{
    ndx <- (i+1):(i+6)
    l2 <- sapply(drugSims, function(sim)
    {
        simArea <- sim$cellArea[simNdx] 
        realArea <- realData[ndx,]$cellArea
        return(l2Norm(simArea, realArea))
    })
    drugFit <- drugSims[[which(l2==min(l2))]]
    realData[(i+1):(i+6),]$Fit <- drugFit$cellArea[simNdx] *
        realData[i+1,]$cellArea / drugFit$cellArea[simNdx][1]

    # store fitted simulations
    dose <- c("10ug", "100ug")[round(i/8)]
    tangFit[[dose]] <- drugFit
}

# plot both fits
print(realData)
realData$dosage <- factor(realData$dosage, labels=c("PBS", "10ug", "100ug"))
fig <- ggplot(realData) +
    geom_point(aes(x=day, y=cellArea, shape=factor(dosage), color=factor(dosage))) + 
    geom_line(data=subset(realData, day != 1), aes(x=day, y=Fit, linetype=factor(dosage), color=factor(dosage))) +
    scale_linetype_manual(values=c("solid", "dashed", "dotdash")) +
    labs(title="Tang Data", linetype="Dosage", shape="Dosage", color="Dosage",
        caption="Figure 3a", x="Day", y="fluorescence")
ggsave(filename="fig3a.pdf", plot=fig)

#### Kagohara Data (3bc)

# get real data
realData <- subset(kagoharaData, Experiment==1)[,c(1,2,3,4,10)]
realData$Fit <- rep(0, nrow(realData))
print(realData)

# loop through cell lines, fit each one
simNdx <- 24 + seq(0,120,24) + 1 # offset 24 for smoother gene expression
kagoharaFit <- list()
noDrugSims <- fig3Data[sapply(fig3Data, function(d) d$drugEffect==1)]
for (i in c(1,7,13)) # day 0 index of each cell line for PBS
{
    # calculate indices
    pbsNdx <- i:(i+5)
    ctxNdx <- 18 + pbsNdx

    # fit pbs data
    l2 <- sapply(noDrugSims, function(sim)
    {
        simArea <- sim$cellArea[simNdx] 
        realArea <- realData[pbsNdx,]$Mean
        return(l2Norm(simArea, realArea))
    })
    pbsFit <- noDrugSims[[which(l2==min(l2))]]
    realData[pbsNdx,]$Fit <- pbsFit$cellArea[simNdx] * realData[i,]$Mean /
        pbsFit$cellArea[simNdx][1]
    
    # fit ctx data
    drugSims <- fig3Data[sapply(fig3Data, function(d)
        d$initDensity==pbsFit$initDensity & d$cycleLength==pbsFit$cycleLength)]
    l2 <- sapply(drugSims, function(sim)
    {
        simArea <- sim$cellArea[simNdx] 
        realArea <- realData[ctxNdx,]$Mean
        return(l2Norm(simArea, realArea))
    })
    drugFit <- drugSims[[which(l2==min(l2))]]
    realData[ctxNdx,]$Fit <- drugFit$cellArea[simNdx] * realData[i,]$Mean /
        drugFit$cellArea[simNdx][1]

    # store fitted simulations
    line <- c("1", "6", "25")[(round(i/7) + 1)]
    print(c(pbsFit$cycleLength, pbsFit$initDensity))
    kagoharaFit[[paste("pbs_scc", line, sep="")]] <- pbsFit
    kagoharaFit[[paste("ctx_scc", line, sep="")]] <- drugFit
}

# plot pbs fit
print(realData)
fig <- ggplot(subset(realData, Treatment=='PBS')) +
    geom_point(aes(x=Day, y=Mean, shape=CellLine, color=CellLine)) + 
    geom_line(aes(x=Day, y=Fit, group=CellLine, color=CellLine)) +
    labs(title="Kagohara Data - PBS",
        caption="Figure 3b", x="Day", y="fluorescence")
ggsave(filename="fig3b.pdf", plot=fig)

# plot ctx fit
fig <- ggplot(subset(realData, Treatment=='CTX')) +
    geom_point(aes(x=Day, y=Mean, shape=CellLine, color=CellLine)) + 
    geom_line(aes(x=Day, y=Fit, group=CellLine, color=CellLine)) +
    labs(title="Kagohara Data - CTX",
        caption="Figure 3c", x="Day", y="fluorescence")
ggsave(filename="fig3c.pdf", plot=fig)

# save fitted simulations
save(tangFit, kagoharaFit, file="Fig3_Fitted_Simulations.RData")
