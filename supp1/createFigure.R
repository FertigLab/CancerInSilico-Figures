library(CancerInSilico)
library(foreach)
library(doParallel)

load('supp1.RData')

# 1a - variance vs initialNum
data_1a <- matrix(as.numeric(data_1a), nrow=nrow(data_1a))

pdf("supp_1a.pdf", width=13, height=10)
plot(seq(10,100,10), apply(data_1a, 2, var), type="l", main="supp1a",
    xlab="initial number of cells", ylab="variance of final population")
dev.off()

## 1b - total vs sync
pdf("supp_1b.pdf", width=13, height=10)
x <- data_1b[[1]]
y_sync <- data_1b[[2]]
y_unsync <- data_1b[[3]]
yall <- c(y_sync, y_unsync)
plot(NULL, xlim=c(min(x), max(x)), ylim=c(min(yall), max(yall)),
    main="supplement 1b", xlab="hours", ylab="num cells")
lines(x, y_sync)
lines(x, y_unsync, col="red")
dev.off()

# 1c
pdf("supp_1c.pdf", width=13, height=10)
plot(data_1c[[1]], data_1c[[2]], type="l", main="supp 1c",
    xlab="nG", ylab="effective cycle (should be 24)")
dev.off()

# 1d
pdf("supp_1d.pdf", width=13, height=10)
plot(data_1d[[1]], data_1d[[2]], type="l", main="supp 1d",
    xlab="epsilon", ylab="effective cycle (should be 24)")
dev.off()

# 1e
pdf("supp_1e.pdf", width=13, height=10)
plot(data_1e[[1]], data_1e[[2]], type="l", main="supp 1e",
    xlab="delta", ylab="effective cycle (should be 24)")
dev.off()
