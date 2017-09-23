d <- read.table("/Users/jtf/Documents/k-mer_proj/mitoverlap.txt", header = TRUE)
d.dis <- dist(d[,3:ncol(d)])
d.cmds <- cmdscale(d.dis, k = 2)
plot(d.cmds)

d.corr <- cor(d[,3:ncol(d)])
require(lattice)
#levelplot(d.corr)