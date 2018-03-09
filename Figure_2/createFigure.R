library(ggplot2)
load("fig2.RData")

fig <- ggplot(subset(df, initialDensity==0.1), aes(x=time)) +
    geom_line(aes(y=nCells, linetype=boundary, group=interaction(cycleLength, boundary))) +
    scale_linetype_manual(values=c("dashed", "solid")) +    
    scale_y_continuous(trans='log10', breaks=c(100,200,500,1000)) +
    labs(title = "Sensitivity to Growth Rate and Boundary Presence", 
        caption = "Figure 2a, expected cycle lengths are 12,28,44 (hrs) from left to right",
        x = "Time", y = "Number Of Cells (log scale)", linetype = "Boundary") 
ggsave(filename='fig2a.pdf', plot=fig)

fig <- ggplot(subset(df, cycleLength==24), aes(x=time)) +
    geom_line(aes(y=density, linetype=boundary, group=interaction(initialDensity, boundary))) +
    scale_linetype_manual(values=c("dashed", "solid")) +
    labs(title = "Boundary Presence Effects Maximum Population Density", 
        caption = "Figure 2b, varying the initial density and holding growth rates constant",
        x = "Time", y = "Population Density", linetype = "Boundary")
ggsave(filename='fig2b.pdf', plot=fig)
