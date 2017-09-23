d1 <- read.table("/Users/jtf/Documents/k-mer_proj/chr01/chr1_dist_shortstep.txt", header = TRUE)
d2 <- read.table("/Users/jtf/Documents/k-mer_proj/chr02/dist02_k5.txt", header = TRUE)
d4 <- read.table("/Users/jtf/Documents/k-mer_proj/chr04/dist04_k5.txt", header = TRUE)
d5 <- read.table("/Users/jtf/Documents/k-mer_proj/chr05/dist05_k5.txt", header = TRUE)
d6 <- read.table("/Users/jtf/Documents/k-mer_proj/chr06/dist06_k5.txt", header = TRUE)
d7 <- read.table("/Users/jtf/Documents/k-mer_proj/chr07/dist07_k5.txt", header = TRUE)
d8 <- read.table("/Users/jtf/Documents/k-mer_proj/chr08/dist08_k5.txt", header = TRUE)
d9 <- read.table("/Users/jtf/Documents/k-mer_proj/chr09/dist09_k5.txt", header = TRUE)
d10 <- read.table("/Users/jtf/Documents/k-mer_proj/chr10/dist10_k5.txt", header = TRUE)
d11 <- read.table("/Users/jtf/Documents/k-mer_proj/chr11/dist11_k5.txt", header = TRUE)
d13 <- read.table("/Users/jtf/Documents/k-mer_proj/chr13/dist13_k5.txt", header = TRUE)
d15 <- read.table("/Users/jtf/Documents/k-mer_proj/chr15/dist15_k5.txt", header = TRUE)
d16 <- read.table("/Users/jtf/Documents/k-mer_proj/chr16/dist16_k5.txt", header = TRUE)
d17 <- read.table("/Users/jtf/Documents/k-mer_proj/chr17/dist17_k5.txt", header = TRUE)
d20 <- read.table("/Users/jtf/Documents/k-mer_proj/chr20/dist20_k5.txt", header = TRUE)
dX <- read.table("/Users/jtf/Documents/k-mer_proj/chrX/distX_k5.txt", header = TRUE)
#d2 <- read.table("/Users/jtf/Documents/k-mer_proj/distoneminus.txt", header = TRUE)
#d.dat <- d[,-4]
#d.vip <- d.dat[d.dat[,3] < 2.61,]venteen
d.vip <- d1[d1[,3] < 2.5,]

plot(1/d1[,3], type="b", cex = 0.3, xlim = c(5900, 6000))

#hist(d[, 3], breaks = 100, main = paste("Histogram of distances between mtdna and chromsome 1 3kb windows"), xlab = "Distance")