# library('CancerInSilico')
library('ggplot2')
library('RColorBrewer')  
library('reshape2') 
library(methods)

makeHeatMap <- function(final_proportion_data) {
  
  ordered_init_prop_typeA <- sort(unique(sapply(final_proportion_data, function(x) x[1])))
  ordered_cycle_length_typeB <- sort(unique(sapply(final_proportion_data, function(x) x[2])))
  
  num_rows <- length(ordered_init_prop_typeA)
  num_cols <- length(ordered_cycle_length_typeB)
  
  final_proportion_mat <- matrix(nrow = num_rows, ncol = num_cols, dimnames = list(ordered_init_prop_typeA, ordered_cycle_length_typeB))
  
  for (data in final_proportion_data) {
    
    xind <- which(ordered_init_prop_typeA == data[1])
    yind <- which(ordered_cycle_length_typeB == data[2])
    
    final_proportion_mat[xind, yind] <- data[3]
  }
  
  return(final_proportion_mat)
  
}

load("Figure_5b_cleaned.RData")
fpm <- makeHeatMap(final_proportion_data)
fpm.melted <- melt(fpm)
hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')
fig <- ggplot(fpm.melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colours = hm.palette(100)) +
  ylab('Cycle Length of Faster Growing Cell Type') +
  xlab('Initial Proportion of Slower Growing Cell Type') +
  ggtitle("Without boundary") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=10))
ggsave(filename='fig5b.png', plot=fig)