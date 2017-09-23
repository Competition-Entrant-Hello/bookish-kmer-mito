
mt5merf.tmp <- read.table("/Users/jtf/Documents/k-mer_proj/mitoverlapR.txt", header=T)
# mtsnpd.tmp  <- read.table("mitoverlap_varcount.txt", header=F)
mtsnpd.tmp  <- scan("/Users/jtf/Documents/k-mer_proj/mitoverlap_varcount.txt")


mt5merf <- mt5merf.tmp[, 3:ncol(mt5merf.tmp)]
mtsnpd <- mtsnpd.tmp/1000

# apply is similar to a loop, but the code is more compact
# 2 means to pick one column at the time (1 is to pick a row)
# once the column is picked, being "x", function defines what to do with that column
# the result is saved to an array (same as your "corr")

# linear correction
 mtcor.5mer.snp <-  apply( mt5merf, 2, function(x){ cor(x, mtsnpd) } )
# non-linear correction (Spearman)
mtncor.5mer.snp <- apply( mt5merf, 2, function(x){ cor(x, mtsnpd, method="spearman") } )

# enough margin on 4 sides
par(mar=c(6,6,6,6))
plot(0,0, pch="", xlim=c(1,32), ylim=c(0, 1.7),
 xlab="window position", ylab="5-mer frequency*100")

j <-0
for(i in 1:(ncol(mt5merf.tmp) - 2)){
 # plot only if corr > 0.5 or corr < -0.5
  if(abs(mtcor.5mer.snp[i]) > 0.5){
 #if(abs(mtncor.5mer.snp[i]) > 0.5){

  # add some noise to separate lines
  xnoise <- rnorm(32, sd=0.1)
  ynoise <- rnorm(32, sd=0.01)

  # here i use color rainbow(512)
  points(1:32+xnoise, mt5merf[,i]+ynoise,  type="l", col=rainbow(512)[i])

  # text is positioned in y direction by j, which increases by 1 after each line 
  text(1, 1.7-j*0.05, i, col=rainbow(512)[i])
  text(5, 1.7-j*0.05, names(mt5merf)[i], col=rainbow(512)[i])
  text(10, 1.7-j*0.05, round(mtcor.5mer.snp[i],3), col=rainbow(512)[i])
  j <- j+1
 }
}

lines(mtsnpd.tmp/100, type = "b", pch = 5)
