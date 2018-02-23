library('ggplot2')
load("Figure_5_cleaned.RData")

nReplicates <- 200
nVars <- 10

getRunID <- function(name, nTypes)
{
    end <- strsplit(name, "output_")[[1]][2]
    id <- as.numeric(strsplit(end, ".RData")[[1]][1])
    if (nTypes == 1)
    {
        return(floor((id - 1) / nReplicates) + 1)
    }
    else
    {
        return(floor((id - (nVars * nReplicates + 1)) / nReplicates) + 1)
    }
}

getReplicateID <- function(name, nTypes)
{
    end <- strsplit(name, "output_")[[1]][2]
    id <- as.numeric(strsplit(end, ".RData")[[1]][1])
    return(((id - 1) %% nReplicates) + 1)
}

processHittingTime <- function(out, name)
{
    hittingTime <- min(which(out$density > 0.60))
    return(c(out$numTypes, getRunID(name, out$numTypes), hittingTime))
}

processRaw <- function(out, name)
{
    N <- length(out$density)
    time <- 1:N - 1
    den <- out$density
    nTypes <- rep(out$numTypes, N)
    runId <- rep(getRunID(name, out$numTypes), N)
    repId <- rep(getReplicateID(name, out$numTypes), N)
    return(cbind(time, den, nTypes, runId, repId))
}

mat <- unname(t(sapply(names(fig5Data), function(name) processHittingTime(fig5Data[[name]], name))))

oneTypeData <- matrix(nrow=nVars, ncol=2)
for (v in 1:nVars)
{
    ndx <- (mat[,1] == 1) & (mat[,2] == v)
    oneTypeData[v,1] <- round(v/3, 2)
    oneTypeData[v,2] <- sd(mat[ndx,3])
}

twoTypeData <- matrix(nrow=nVars, ncol=2)
for (v in 1:nVars)
{
    ndx <- (mat[,1] == 2) & (mat[,2] == v)
    twoTypeData[v,1] <- 2 * 0.05 * v / sqrt(12)
    twoTypeData[v,2] <- sd(mat[ndx,3])
}

scale <- function(data, scale_data)
{
    mn <- min(scale_data)
    mx <- max(scale_data)
    fac <- (mx - mn) / (max(data) - min(data))
    normed <- data - min(data)
    return(normed * fac + mn)
}

rawList <- lapply(names(fig5Data), function(n) processRaw(fig5Data[[n]], n))
N <- length(rawList)
rawMat <- do.call(rbind, rawList[sample(1:N, 800)])
rawDF <- data.frame(time=rawMat[,1], density=rawMat[,2], nTypes=rawMat[,3],
    run_id=rawMat[,4], rep_id=rawMat[,5])

fig <- ggplot(subset(rawDF, nTypes==1 & run_id==1)) +
    geom_line(aes(x=time, y=density, group=rep_id)) +
    labs(title="Variance of Growth Rate within Single Cell Type", 
        caption="Figure 5a, single cell type with small growth rate variance",
        x="time (hrs)", y="Population Density")
ggsave(filename='fig5a.pdf', plot=fig)

fig <- ggplot(subset(rawDF, nTypes==1 & run_id==10)) +
    geom_line(aes(x=time, y=density, group=rep_id)) +
    labs(title="Variance of Growth Rate within Single Cell Type", 
        caption="Figure 5b, single cell type with large growth rate variance",
        x="time (hrs)", y="Population Density")
ggsave(filename='fig5b.pdf', plot=fig)

fig <- ggplot(subset(rawDF, nTypes==2 & run_id==10)) +
    geom_line(aes(x=time, y=density, group=rep_id)) +
    labs(title="Variance of Initial Cell Type Distribution", 
        caption="Figure 5c, two cell types with a large variance in the initial proportions",
        x="time (hrs)", y="Population Density")
ggsave(filename='fig5c.pdf', plot=fig)


df <- data.frame(xa=c(oneTypeData[,1], scale(twoTypeData[,1], oneTypeData[,1])),
    ya=c(oneTypeData[,2], twoTypeData[,2]), nTypes=factor(c(rep(1,nVars), rep(2,nVars))))

fig <- ggplot(df, aes(x=xa, y=ya, linetype=nTypes)) + geom_line() +
    scale_linetype_manual(values=c("dashed", "solid")) +
    labs(title="Heterogeneity within vs between tumors", 
        caption="Figure 5d, variance in growth rate for single cell type, and initial proportion for two cell types",
        x="Growth Rate/Cell Type Variance", y="Population Growth Variance", linetype="N Cell Types")
ggsave(filename='fig5d.pdf', plot=fig)