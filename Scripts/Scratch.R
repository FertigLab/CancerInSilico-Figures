setwd("~/Biostats_Research/CancerInSilico/Figures/Scripts")

## load data
load('../Data/inheritGrowth.RData')

## run time of each simulation
time <- inherited@params[['runTime']]

## total number of cells over time
inheritedTotal <- sapply(0:time, getNumberOfCells, model=inherited)
randomTotal <- sapply(0:time, getNumberOfCells, model=random)

## plot difference in total
plot(0:time, randomTotal, type = 'l')
lines(0:time, inheritedTotal, col = 'red')

## plot distribution shift  
plot(density(getCycleLengths(inherited, 0)))
lines(density(getCycleLengths(inherited, 168)), col = 'red')

## difference in means
mean(getCycleLengths(inherited, 0)
getCycleLengths(inherited, 168)
