library('xlsx')

processSheet <- function(sheetName='')
{
    data <- read.xlsx('Cetuximab data Across PDX.xlsx',
        sheetName = sheetName, startRow=2)
    trashCols <- which(colnames(data) == "NA.")
    data <- data[,setdiff(1:ncol(data), trashCols)]
    return(data)
}
pdxModels <- c('409_PDX_F2', '428_PDX_F3', '465_PDX_F3',
    '425_PDX_F4', '425_PDX_F6')

getID <- function(str) as.numeric(substr(str, nchar(str), nchar(str)))

day <- c()
model <- c()
treatment <- c()
mouseID <- c()
tumorVol <- c()

for (mod in pdxModels)
{
    rawData <- processSheet(mod)
    for (c in 2:ncol(rawData))
    {
        N <- length(rawData[,1])
        day <- c(day, rawData[,1])
        model <- c(model, rep(mod, N))
        if (grepl('PBS', colnames(rawData)[c]))
            treatment <- c(treatment, rep('PBS', N))
        else
            treatment <- c(treatment, rep('CET', N))
        mouseID <- c(mouseID, rep(getID(colnames(rawData)[c]), N))
        tumorVol <- c(tumorVol, rawData[,c] / rawData[1,c])
    }
}

pdxData <- data.frame(day=day, model=model, treatment=treatment,
    mouseID=mouseID, tumorVol=tumorVol)
save(pdxData, file='PDXData_cleaned.RData')
