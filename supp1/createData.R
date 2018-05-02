library(CancerInSilico)
library(foreach)
library(doParallel)

#load("supp1.RData")

# 1a - variance vs initialNum
nReps <- 12
cl <- makeCluster(6)
registerDoParallel(cl)
data_1a <- sapply(seq(10,100,10), function(n)
{
    finalCount <- foreach(i = 1:nReps, .packages="CancerInSilico") %dopar%
    {
        mod <- inSilicoCellModel(n, 72, 0.1, randSeed=i)
        getNumberOfCells(mod, mod@runTime) / n
    }
    return(finalCount)
})
stopCluster(cl)

# 1b - total vs sync
runTime <- 72
mod_sync <- inSilicoCellModel(initialNum=20, runTime=runTime, density=0.1,
    syncCycles=FALSE)
mod_unsync <- inSilicoCellModel(initialNum=20, runTime=runTime, density=0.1,
    syncCycles=TRUE)
x <- 0:runTime
y_sync <- sapply(x, getNumberOfCells, model=mod_sync)
y_unsync <- sapply(x, getNumberOfCells, model=mod_unsync)
data_1b <- list(x, y_sync, y_unsync)

getDivisionTimes <- function(model, cell)
{
    axis_length <- sapply(0:model@runTime, getAxisLength, model=model, cell=cell)
    axis_length_cleaned <- axis_length[!is.na(axis_length)]
    l_shift <- axis_length_cleaned[1:(length(axis_length_cleaned)-1)]
    r_shift <- axis_length_cleaned[2:length(axis_length_cleaned)]

    if (length(axis_length_cleaned) > 1)
    {
        divisions <- which(r_shift < l_shift)
        return(divisions + sum(is.na(axis_length)))
    }
    else
    {
        return(integer(0))
    }
}

getMeanDivisionTime <- function(model)
{
    times <- sapply(1:getNumberOfCells(model, model@runTime), function(cell_id)
    {
        div <- getDivisionTimes(model, cell_id)
        if (length(div) > 1)
            mean(div[2:length(div)] - div[1:(length(div)-1)])
        else
            NA
    })
    return(mean(times, na.rm=TRUE))
}

# 1c - obs cycle length vs nG
x <- seq(4,48,11)
y <- sapply(x, function(ng)
{
    mod <- inSilicoCellModel(initialNum=50, runTime=72, density=0.01, nG=ng)
    getMeanDivisionTime(mod)
})
data_1c <- list(x,y)

# 1d - obs cycle length vs epsilon
x <- seq(1,50,10)
y <- sapply(x, function(e)
{
    mod <- inSilicoCellModel(initialNum=50, runTime=72, density=0.01, epsilon=e)
    getMeanDivisionTime(mod)
})
data_1d <- list(x,y)

# 1e - obs cycle length & density vs delta
x <- seq(0.1, 0.4, 0.01)
y <- sapply(x, function(d)
{
    mod <- inSilicoCellModel(initialNum=50, runTime=72, density=0.01, delta=d)
    getMeanDivisionTime(mod)
})
data_1e <- list(x,y)

save(data_1a, data_1b, data_1c, data_1d, data_1e, file="supp1.RData")

warnings()