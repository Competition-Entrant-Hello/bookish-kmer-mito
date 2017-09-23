
#################### by chr1 order
chr1.order <- sort.list(-j13[,2])

par(mar=c(6,6,6,6))
plot(1:512, j13[chr1.order, 2]/100, ylim=c(-0.0015, 0.01), xlab="rank according chr1 5-mer freq", ylab="5-mer freq")
points(1:512, j13[chr1.order, 3]/100, type="b", col="red")

for(i in c(1:512) ){
 ii <- chr1.order[i]
 im1 <- i-1
 text(i, -0.001 + (im1%%4)*0.002 , as.character(j13[ii,1]), srt=90, crt=90, cex=0.7)
}



#################### by MT order
mt.order <- sort.list(-j13[,3])

par(mar=c(6,6,6,6))
plot(1:512, j13[mt.order, 3]/100, ylim=c(-0.0015, 0.01), xlab="rank according MT 5-mer freq", ylab="5-mer freq")
points(1:512, j13[mt.order, 2]/100, type="b", col="red")



for(i in c(1:512) ){
 ii <- mt.order[i]
 im1 <- i-1
 text(i, -0.001 + (im1%%4)*0.002 , as.character(j13[ii,1]), srt=90, crt=90, cex=0.7)
}

