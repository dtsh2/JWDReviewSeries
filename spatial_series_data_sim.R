x<-runif(100)
y<-runif(100)
par(pty="s")
plot(x,y,pch=16)
library(spatstat)
#install.packages("spatstat")
library(spatial)

spat.t<-ppp(x, y, c(0,1), c(0,1))
plot(spat.t)
plot(Kest(spat.t))

plot(envelope(spat.t,Kest))
##
# non-random
x1<-c(runif(50))
x2<-seq(0,1,by=0.1)
x3<-runif(15,min=0.1,max=0.2)
y1<-c(runif(50))
y2<-seq(0,1,by=0.1)
y3<-runif(15,min=0.1,max=0.2)
xall<-c(x1,x2,x3)
yall<-c(y1,y2,y3)
par(pty="s")
par(mfrow=c(2,2))
plot(xall,yall,pch=16,xlim=c(0,1),ylim=c(0,1))
plot(x1,y1,pch=16,xlim=c(0,1),ylim=c(0,1))
plot(x2,y2,pch=16,xlim=c(0,1),ylim=c(0,1))
plot(x3,y3,pch=16,xlim=c(0,1),ylim=c(0,1))

spat.t<-ppp(xall, yall, c(0,1), c(0,1))
plot(spat.t)
plot(Kest(spat.t))

plot(envelope(spat.t,Kest))

## add serological data
pos<-rpois(100,lambda=15)
neg<-100-pos
tot<-neg+pos
seropos<-pos/tot
#plot(seropos)

marks(spat.t) <- seropos
plot((spat.t),fg='grey',bg="grey",main="",
     ylab="Y",'xlab="x')

pos<-seq(from=0.1,to=10,by=0.1)
neg<-100-pos
tot<-neg+pos
seropos<-pos/tot

w<-rlnorm(100,sd=0.5)
Z<-seq(from=0.1,to=10,by=0.1)+w
W<-rev(seq(from=0.1,to=10,by=0.1)+w)
V<-rep(1,100)+w
lin.m<-lm(seropos ~ Z)   
summary(lin.m)
plot(seropos,Z)

lin.m<-lm(seropos ~ W)   
summary(lin.m)
plot(seropos,W)

lin.m<-lm(seropos ~ V)   
summary(lin.m)
plot(seropos,V)

lin.m<-lm(seropos ~ V*W*Z)   
summary(lin.m)

marks(spat.t) <- seropos
plot((spat.t),fg='grey',bg="grey",main="",
     ylab="Y",'xlab="x')

mean(seropos)
diff<-seropos-mean(seropos)

rbPal <- colorRampPalette(c('blue','red'))
dat<-data.frame(x,y,diff,seropos)
dat$Col <- rbPal(10)[as.numeric(cut(diff,breaks = 10))]
plot(dat$x,dat$y,pch = 20,col = dat$Col,cex=2)
text(dat$x, dat$y, round(dat$diff, 2), cex=0.8)

dat$Col <- rbPal(10)[as.numeric(cut(dat$seropos,breaks = 10))]
plot(dat$x,dat$y,pch = 20,col = dat$Col,cex=2)
text(dat$x, dat$y, round(dat$seropos, 2), cex=0.8)

plot(dat$x,dat$y,pch = 20,col = "white",cex=2)
text(dat$x, dat$y, round(dat$seropos, 2), cex=0.8)

plot(dat$x,dat$y,pch = 20,col = "white",cex=2)
text(dat$x, dat$y, round(dat$diff, 2), cex=0.8)

dat<-data.frame(x,y,diff,seropos,Z)
dat$Col <- rbPal(10)[as.numeric(cut(dat$Z,breaks = 10))]
plot(dat$x,dat$y,pch = 20,col = dat$Col,cex=2)
text(dat$x, dat$y, round(dat$Z, 2), cex=0.8)
plot(dat$Z,dat$seropos)

lmdat<-lm(dat$seropos~dat$Z+dat$y+dat$x)
summary(lmdat)

