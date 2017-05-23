library(CancerInSilico)
library(methods)

## Figure 1 - General Parameter Sensitivity ##

# a - visual representation of cells (moderately dense)

mod1a <- runCellSimulation(100, 10, 0.4)

# b - multiple growth curves w/ and w/out boundary

mod1b <- list()
mod1b[['noBound']] <- list()
mod1b[['yesBound']] <- list()
cycleLengths <- seq(12,48,6)
initialNum <- 80
runTime <- 168
density <- 0.2

for (cl in cycleLengths)
{
    ct <- new('CellType', name='DEFAULT', minCycle=cl,
        cycleLength=function(x) cl)
    mod1b[['noBound']][[cl]] <- runCellSimulation(initialNum=initialNum,
        runTime=runTime, density=density, cellTypes=c(ct), boundary=0)
    mod1b[['yesBound']][[cl]] <- runCellSimulation(initialNum=initialNum,
        runTime=runTime, density=density, cellTypes=c(ct), boundary=1)
}

# c - density w/ and w/out boundary

mod1c <- list()
mod1c[['noBound']] <- list()
mod1c[['yesBound']] <- list()
density <- seq(0.05, 0.45, 0.05)
initialNum <- 80
runTime <- 168

for (den in density)
{
    mod1c[['noBound']][[den*100]] <- runCellSimulation(initialNum=initialNum,
        runTime=runTime, density=den, boundary=0)
    mod1c[['yesBound']][[den*100]] <- runCellSimulation(initialNum=initialNum,
        runTime=runTime, density=den, boundary=1)
}

# save models

save(mod1a, mod1b, mod1c, file='Figure1_Data.RData')
