require(lattice)

d <- read.table("/Users/jtf/Documents/k-mer_proj/mitoverlapR.txt", header = TRUE)
data <- d[,3:ncol(d)]

#plot.new()
#plot(toplot[,1])
#for(i in seq (1, length(data), 1) ) {if }
dvars <- read.table("/Users/jtf/Documents/k-mer_proj/mitoverlap_varcount.txt", header = FALSE)
all <- cbind(data, dvars)
corr <- cor(all)
varcor <- corr[,ncol(corr)]
top <- varcor[abs(varcor) > 0.5]


matplot(data, type = 'l', lty = 1, xlab = "Window", ylab = "Percentage in Sequence", main = paste("Each 5mer across 32 overlapping windows of mtdna"))
#plot.new()

for ( i in seq(1,length( data ),1) ){
  if(abs(corr[i,513]) > 0.55 ){
      lines(data[,i],type="l", lwd = 5, col = "black")}}


#matplot(data[data > 1.0], type = 'l', lty = 1, xlab = "Window", ylab = "Percentage in Sequence", main = paste("Greatest kmers across 32 overlapping windows of mtdna"))
#text(data[data > 1.35], labels = names(data))


