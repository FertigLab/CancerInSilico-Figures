setwd('~/FertigLab/CancerInSilico/Figures/')
load('Figure_4/Figure_4_cleaned.RData')
load('Figure_6/Figure_6_cleaned.RData')

fig4_pbs_microarray <- ge_pbs$expression
fig4_10ug_microarray <- ge_10ug$expression
fig4_100ug_microarray <- ge_100ug$expression

fig6_two_cell_types_single_cell_RNAseq <- ge$expression
fig6_two_cell_types_bulk_RNAseq <- ge_bulk$expression

colnames(ge$expression)[4128]


getType <- function(str)
{
  pos <- regexpr('_', str)[1]
  cell <- as.numeric(substr(str, 2, pos-1))
  return(cellType[cell])
}

getPhase <- function(str)
{
  pos <- regexpr('_', str)[1]
  cell <- as.numeric(substr(str, 2, pos-1))
  time <- as.numeric(substr(str, pos + 2, nchar(str)))
  return(cellPhase[cell, time+1])
}

fig6_single_cell_types <- unname(sapply(colnames(fig6_two_cell_types_single_cell_RNAseq), getType))
fig6_single_cell_phases <- unname(sapply(colnames(fig6_two_cell_types_single_cell_RNAseq), getPhase))

save(fig4_pbs_microarray, fig4_10ug_microarray, fig4_100ug_microarray, file='fig4_expression.RData')
save(fig6_two_cell_types_bulk_RNAseq, fig6_two_cell_types_single_cell_RNAseq, file='fig6_expression.RData')


load('fig4_expression.RData')
load('fig6_expression.RData')
