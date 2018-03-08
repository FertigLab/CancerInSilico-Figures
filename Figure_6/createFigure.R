# load libraries

library(CancerInSilico)
library(ggplot2)

# simulate data

typeA <- new('CellType', name='A', minCycle=24-4,
    cycleLength=function() max(24-4, rnorm(1,24,1)))

typeB <- new('CellType', name='B', minCycle=36-4,
    cycleLength=function() max(36-4, rnorm(1,36,1)))

twoTypesModel <- inSilicoCellModel(initialNum=100, runTime=168, density=0.05,
    boundary=1, syncCycle=FALSE, randSeed=123, outputIncrement=4,
    recordIncrement=0.25, timeIncrement=0.001, cellTypes=c(typeA, typeB),
    cellTypeInitFreq=c(0.5,0.5), maxDeformation=0.1,
    maxTranslation=0.1, maxRotation=0.3, nG=28, epsilon=10, delta=0.2)

# save data
save(twoTypesModel, file="Fig6_RawData.RData")
