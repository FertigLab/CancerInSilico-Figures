library(ggplot2)
load('Figure_3_cleaned.RData') #fig3Data
load('tangData.RData') #tangData
load('KagoharaGrowthData.RData')

# useful functions
l2norm <- function(a,b)
{
    sqrt(sum((a-b)^2))
}

computeL2 <- function(sim, real, ndx)
{
    nCells <- sim$numCells[ndx]
    nCells <- nCells * real[1] / nCells[1]
    l2norm(nCells, real)
}

#### Tang Data (3a)

# get real data
realData <- subset(tangData, dosage %in% c(0, 10, 100))
realData$Fit <- rep(0, nrow(realData))

# time points of real data translated to indices in simulated data
simNdx <- seq(0,144,24) + 1
fitNdx <- seq(24,144,24) + 1

# fit pbs data
noDrugSims <- fig3Data[sapply(fig3Data, function(d) d$drugEffect==1)]
l2 <- sapply(noDrugSims, computeL2, real=realData[2:7,]$numCells, ndx=fitNdx)
pbsFit <- noDrugSims[[which(l2==min(l2))]]
realData[1:7,]$Fit <- pbsFit$numCells[simNdx] * realData[2,]$numCells /
    pbsFit$numCells[simNdx][2]

tangFit <- list("pbs"=pbsFit)

# fit ctx data
drugSims <- fig3Data[sapply(fig3Data, function(d)
    d$initDensity==pbsFit$initDensity & d$cycleLength==pbsFit$cycleLength)]
for (i in c(8,15)) # day 1 index for ctx
{
    ndx <- (i+1):(i+6)
    l2 <- sapply(drugSims, computeL2, real=realData[ndx,]$numCells, ndx=fitNdx)
    drugFit <- drugSims[[which(l2==min(l2))]]
    realData[i:(i+6),]$Fit <- drugFit$numCells[simNdx] * realData[i+1,]$numCells /
        drugFit$numCells[simNdx][2]

    # store fitted simulations
    dose <- c("10ug", "100ug")[round(i/8)]
    tangFit[[dose]] <- drugFit
}

# plot both fits
realData$dosage <- factor(realData$dosage, labels=c("PBS", "10ug", "100ug"))
fig <- ggplot(realData) +
    geom_point(aes(x=day, y=numCells, shape=factor(dosage), color=factor(dosage))) + 
    geom_line(aes(x=day, y=Fit, linetype=factor(dosage), color=factor(dosage))) +
    scale_linetype_manual(values=c("solid", "dashed", "dotdash")) +
    labs(title="Tang Data", linetype="Dosage", shape="Dosage", color="Dosage",
        caption="Figure 3a", x="Day", y="Number Of Cells")
ggsave(filename="fig3a.pdf", plot=fig)

#### Kagohara Data (3bc)

# get real data
realData <- subset(kagoharaData, Experiment==1)[,c(1,2,3,4,10)]
realData$Fit <- rep(0, nrow(realData))

# time points of real data translated to indices in simulated data
simNdx <- seq(0,120,24) + 1

kagoharaFit <- list()
noDrugSims <- fig3Data[sapply(fig3Data, function(d) d$drugEffect==1)]
for (i in c(1,7,13)) # day 0 index of each cell line for PBS
{
    # calculate indices
    pbsNdx <- i:(i+5)
    ctxNdx <- 18 + pbsNdx

    # fit pbs data
    l2 <- sapply(noDrugSims, computeL2, real=realData[pbsNdx,]$Mean, ndx=simNdx)
    pbsFit <- noDrugSims[[which(l2==min(l2))]]
    realData[pbsNdx,]$Fit <- pbsFit$numCells[simNdx] * realData[i,]$Mean /
        pbsFit$numCells[simNdx][1]
    
    # fit ctx data
    drugSims <- fig3Data[sapply(fig3Data, function(d)
        d$initDensity==pbsFit$initDensity & d$cycleLength==pbsFit$cycleLength)]
    l2 <- sapply(drugSims, computeL2, real=realData[ctxNdx,]$Mean, ndx=simNdx)
    drugFit <- drugSims[[which(l2==min(l2))]]
    realData[ctxNdx,]$Fit <- drugFit$numCells[simNdx] * realData[i,]$Mean /
        drugFit$numCells[simNdx][1]

    # store fitted simulations
    line <- c("1", "6", "25")[(round(i/7) + 1)]
    kagoharaFit[[paste("pbs_scc", line, sep="")]] <- pbsFit
    kagoharaFit[[paste("ctx_scc", line, sep="")]] <- drugFit
}

# plot pbs fit
fig <- ggplot(subset(realData, Treatment=='PBS')) +
    geom_point(aes(x=Day, y=Mean, shape=CellLine, color=CellLine)) + 
    geom_line(aes(x=Day, y=Fit, group=CellLine, color=CellLine)) +
    labs(title="Kagohara Data - PBS",
        caption="Figure 3b", x="Day", y="Number Of Cells")
ggsave(filename="fig3b.pdf", plot=fig)

# plot ctx fit
fig <- ggplot(subset(realData, Treatment=='CTX')) +
    geom_point(aes(x=Day, y=Mean, shape=CellLine, color=CellLine)) + 
    geom_line(aes(x=Day, y=Fit, group=CellLine, color=CellLine)) +
    labs(title="Kagohara Data - CTX",
        caption="Figure 3c", x="Day", y="Number Of Cells")
ggsave(filename="fig3c.pdf", plot=fig)

# save fitted simulations
save(tangFit, kagoharaFit, file="Fig3_Fitted_Simulations.RData")
