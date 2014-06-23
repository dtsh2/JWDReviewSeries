## mean 80, sd = 8
## CI 64 and 96 (+/- 2SD)

# prob selecting sample outside these values, based on these and different sampling
par(mfrow=c(2,2))
# calc z values
z1<-(64-80)/8
z2<-(96-80)/8
pnorm(z2)-pnorm(z1)


yvals10<-rnorm(10,80,8)
 plot(density(yvals10))
 rug(yvals10)
abline(v=64)
abline(v=96)
yvals100<-rnorm(100,80,8)
 plot(density(yvals100))
 rug(yvals100)
abline(v=64)
abline(v=96)
yvals1000<-rnorm(1000,80,8)
 plot(density(yvals1000))
 rug(yvals1000)
abline(v=64)
abline(v=96)
yvals10000<-rnorm(10000,80,8)
  plot(density(yvals10000))
  rug(yvals10000)
abline(v=64)
abline(v=96)

par(mfrow=c(2,3))

yvals10<-rnorm(10,80,8)
plot(density(yvals10),main="",xlab="")
rug(yvals10,side=1,col="red",lwd=1,ticksize=0.2)
abline(v=64)
abline(v=96)
yvals20<-rnorm(20,80,8)
plot(density(yvals20),main="",xlab="")
rug(yvals20,side=1,col="red",lwd=1,ticksize=0.2)
abline(v=64)
abline(v=96)
yvals30<-rnorm(30,80,8)
plot(density(yvals30),main="",xlab="")
rug(yvals30,side=1,col="red",lwd=1,ticksize=0.2)
abline(v=64)
abline(v=96)
yvals40<-rnorm(40,80,8)
plot(density(yvals40),main="",xlab="")
rug(yvals40,side=1,col="red",lwd=1,ticksize=0.2)
abline(v=64)
abline(v=96)
yvals50<-rnorm(50,80,8)
plot(density(yvals50),main="",xlab="")
rug(yvals50,side=1,col="red",lwd=1,ticksize=0.2)
abline(v=64)
abline(v=96)
yvals100<-rnorm(100,80,8)
plot(density(yvals100),main="",xlab="")
rug(yvals100,side=1,col="red",lwd=1,ticksize=0.2)
abline(v=64)
abline(v=96)


plot(yvals100, col= ifelse(yvals100 >= 96, "red", 
                         ifelse(yvals100 <= 64,"blue", "black")),pch=16)
abline(h=64,col="blue")
abline(h=96,col="red")

plot(yvals20, col= ifelse(yvals20 >= 96, "red", 
                           ifelse(yvals20 <= 64,"blue", "black")),pch=16)
abline(h=64,col="blue")
abline(h=96,col="red")

hist(yvals100,breaks=seq(50,120,2))
abline(v=64,col="red")
abline(v=96,col="red")

hist(yvals20,breaks=seq(50,120,2))
abline(v=64,col="red")
abline(v=96,col="red")


plot(yvals1000, col= ifelse(yvals1000 >= 96, "red", 
                          ifelse(yvals1000 <= 64,"blue", "black")),pch=16)
abline(h=64,col="blue")
abline(h=96,col="red")

hist(yvals1000,breaks=seq(50,120,2))
abline(v=64,col="red")
abline(v=96,col="red")

## adjust p values

a <- 5 # mean - use the same one
s <- 2 # sd
n <- 1 # sample size
xbar <- rnorm(10,mean=a,sd=s) # means randomly drawn from the above normal distribution
pvals<-2*(1-pnorm(xbar,mean=a,sd=s/sqrt(20)))
hist(pvals)

padj<-p.adjust(pvals, method = "bonferroni", n = length(pvals))

# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")

#padj<-p.adjust(pvals, method = "BH", n = length(pvals))

plot(pvals)
abline(h=0.05)

plot(padj)
abline(h=0.05)

padj

pv<-do.call( rbind, replicate(10, t.test(rnorm(2,mean=10,sd=s),rnorm(100,mean=10,sd=s)), simplify=FALSE ) )
pvals<-(as.vector(as.numeric(pv[,3])))

plot(pvals,ylim=c(0,max(pvals)))
abline(h=0.05)

padj<-p.adjust(pvals, method = "bonferroni", n = length(pvals))
plot(padj,ylim=c(0,max(padj)))
abline(h=0.05)

