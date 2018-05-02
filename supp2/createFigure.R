library(CancerInSilico)
data(SamplePathways)

# supp2a = M/S phase activity vs percentage of cells in the cycle
logistic <- function(k, x0, x) 1 / (1 + exp(-k * (x - x0)))
x <- seq(0,0.3,0.01)
y <- logistic(5 / 0.1, 0.05, x)
pdf("supp_2a.pdf", width=13, height=10)
plot(x, y, type="l", xlab="Proportion of Cells in Phase (Mitosis or Interphase)",
    ylab="Pathway Activity")
dev.off()

#'supp2b - growth activity vs cyle length
x <- seq(8,72,1)
y <- exp(-1 * x / 48)
pdf("supp_2b.pdf", width=13, height=10)
plot(x, y, type="l", xlab="Cycle Length (hrs)", ylab="Pathway Activity")
dev.off()