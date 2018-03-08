library('ggplot2')
load("Figure_5_cleaned.RData")

nReplicates <- 40
nVars <- 10

#### Figures 5abcd

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

processRaw <- function(out, name)
{
    N <- length(out$numCells)
    time <- 1:N - 1
    nCells <- out$numCells
    nTypes <- rep(out$numTypes, N)
    runId <- rep(getRunID(name, out$numTypes), N)
    repId <- rep(getReplicateID(name, out$numTypes), N)
    return(cbind(time, nCells, nTypes, runId, repId))
}

rawList <- lapply(names(fig5Data), function(n) processRaw(fig5Data[[n]], n))
rawMat <- do.call(rbind, rawList)
rawDF <- data.frame(time=rawMat[,1], nCells=rawMat[,2], nTypes=rawMat[,3],
    run_id=rawMat[,4], rep_id=rawMat[,5])

fig <- ggplot(subset(rawDF, nTypes==1 & run_id==1)) +
    geom_line(aes(x=time, y=nCells, group=rep_id)) +
    labs(title="Variance of Growth Rate within Single Cell Type", 
        caption="Figure 5a, single cell type with small growth rate variance",
        x="time (hrs)", y="Population Density")
ggsave(filename='fig5a.pdf', plot=fig)

fig <- ggplot(subset(rawDF, nTypes==1 & run_id==10)) +
    geom_line(aes(x=time, y=nCells, group=rep_id)) +
    labs(title="Variance of Growth Rate within Single Cell Type", 
        caption="Figure 5b, single cell type with large growth rate variance",
        x="time (hrs)", y="Population Density")
ggsave(filename='fig5b.pdf', plot=fig)

fig <- ggplot(subset(rawDF, nTypes==2 & run_id==1)) +
    geom_line(aes(x=time, y=nCells, group=rep_id)) +
    labs(title="Variance of Initial Cell Type Distribution", 
        caption="Figure 5c, two cell types with a small variance in the initial proportions",
        x="time (hrs)", y="Population Density")
ggsave(filename='fig5c.pdf', plot=fig)

fig <- ggplot(subset(rawDF, nTypes==2 & run_id==10)) +
    geom_line(aes(x=time, y=nCells, group=rep_id)) +
    labs(title="Variance of Initial Cell Type Distribution", 
        caption="Figure 5c, two cell types with a large variance in the initial proportions",
        x="time (hrs)", y="Population Density")
ggsave(filename='fig5d.pdf', plot=fig)

#### Figure 5e

mat <- matrix(nrow=(nVars * 2), ncol=3)
finalTime <- max(rawDF$time)
for (n in 1:2)
{
    for (v in 1:nVars)
    {
        temp <- subset(rawDF, time==finalTime & nTypes==n & run_id==v)
        row <- v + (n - 1) * nVars
        mat[row,1] <- n
        mat[row,2] <- v
        mat[row,3] <- var(temp$nCells)
    }
}
df <- data.frame(nTypes=mat[,1], init_var=mat[,2], final_var=mat[,3])

fig <- ggplot(df, aes(x=init_var, y=final_var, linetype=factor(nTypes))) + geom_line() +
    scale_linetype_manual(values=c("dashed", "solid")) +
    labs(title="Heterogeneity within vs between tumors", 
        caption="Figure 5d, variance in growth rate for single cell type, and initial proportion for two cell types",
        x="Growth Rate/Cell Type Variance", y="Population Growth Variance", linetype="N Cell Types")
ggsave(filename='fig5e.pdf', plot=fig)
