## clear environment and install packages
rm(list=ls())
library(spatstat)
library(spatial)

## create coordinates for 100 random points
x <- runif(100)
y <- runif(100)

## format plot area and plot above points
par(pty="s")
plot(x, y, pch=16)

## bundle x & y into one point pattern object for spatstat library and plot
spat.t<-ppp(x, y, c(0,1), c(0,1))
plot(spat.t)

## plot Ripley's K graphical description
## (how likely another plot point is within a given radius)
plot(Kest(spat.t, correction="isotropic"), 
     legendargs=list(cex=0.8, y.intersp=0.4))
plot(envelope(spat.t, Kest, correction="isotropic"), 
     legendargs=list(cex=0.8, y.intersp=0.4))

## create some non-random data (100 coordinates)
x1<-c(runif(65)) # random
y1<-c(runif(65))
x2<-seq(0, 1, by=0.07) # trend
y2<-seq(0, 1, by=0.07)
x3<-runif(20, min=0.1, max=0.2) # cluster
y3<-runif(20, min=0.1, max=0.2)
xall<-c(x1, x2, x3) # combine
yall<-c(y1, y2, y3)

## format and plot components and combination
par(pty="s")
par(mfrow=c(2, 2))
plot(x1, y1, pch=16, xlim=c(0, 1), ylim=c(0, 1))
plot(x2, y2, pch=16, xlim=c(0, 1), ylim=c(0, 1))
plot(x3, y3, pch=16, xlim=c(0, 1), ylim=c(0, 1))
plot(xall, yall, pch=16, xlim=c(0, 1), ylim=c(0, 1))

## bundle x & y into one point pattern object for spatstat library and plot
par(mfrow=c(1,1))
spat.t<-ppp(xall, yall, c(0,1), c(0,1))
plot(spat.t)

## plot Ripley's K graphical description
## Here our plot (black line) falls outside range of perfectly random (red line)
plot(Kest(spat.t, correction="isotropic"), 
     legendargs=list(cex=0.8, y.intersp=0.4))
plot(envelope(spat.t, Kest, correction="isotropic"), 
     legendargs=list(cex=0.8, y.intersp=0.4))

## add serological data
# 100 measures, mean seropositive results = 15, random
pos <- rpois(100, lambda=15)
neg <- 100 - pos
tot <- neg + pos
seropos <- pos / tot # proportion seropositive
plot(seropos)

## make seropos data show as the size of the points from our spatial data
marks(spat.t) <- seropos
plot((spat.t),fg='grey',bg="grey",main="",
     ylab="Y",'xlab="x')

## nonrandom seropos data
pos <- seq(from=0.1, to=10, by=0.1)
neg <- 100 - pos
tot <- neg + pos
seropos <- pos / tot

## covariate data, 100 measures
## w = lognormal random, 
## Z = ascending sequence + some random, 
## W = descending sequence + some random, 
## V = same random data as w, larger values
w <- rlnorm(100, sd=0.5) 
Z <- seq(from=0.1, to=10, by=0.1) + w
W <- rev(seq(from=0.1, to=10, by=0.1) + w)
V <- rep(1, 100) + w

## linear models and plots of seropos data vs. above variables
par(mfrow=c(2,2))
lin.m <- lm(seropos ~ Z)   
summary(lin.m)
plot(seropos, Z)

lin.m<-lm(seropos ~ W)   
summary(lin.m)
plot(seropos, W)

lin.m<-lm(seropos ~ V)   
summary(lin.m)
plot(seropos, V)

lin.m<-lm(seropos ~ V * W * Z)   
summary(lin.m)

## create a heatmap based on difference from mean
## red = seroprevalence above mean, blue = seroprevalence below mean
# format plot area and define difference, put coordinates, seropos, 
# and difference in data.frame
par(mfrow=c(1, 2))
diff <- seropos - mean(seropos)
dat <- data.frame(x, y, diff, seropos, Z)
# use diff as plot colors and add labels to points
rbPal <- colorRampPalette(c('blue', 'red'))
dat$Col <- rbPal(10)[as.numeric(cut(diff, breaks = 10))]
plot(dat$x, dat$y, pch = 20, col = dat$Col, cex=2)
text(dat$x, dat$y, round(dat$diff, 2), cex=0.8) 

## take away the colored points, only labels shown
plot(dat$x,dat$y,pch = 20,col = "white",cex=2)
text(dat$x, dat$y, round(dat$diff, 2), cex=0.8)

## create a heatmap based on actual data (only labels should change)
## red = high seropos, blue = low seropos
dat$Col <- rbPal(10)[as.numeric(cut(dat$seropos, breaks = 10))]
plot(dat$x, dat$y, pch = 20, col = dat$Col, cex=2)
text(dat$x, dat$y, round(dat$seropos, 2), cex=0.8)

## take away the colored points, only labels shown
plot(dat$x, dat$y, pch = 20, col = "white", cex=2)
text(dat$x, dat$y, round(dat$seropos, 2), cex=0.8)


## show comparison between using seropos (perfect ascending linear data)
## and using Z (ascending linear with some randomness) as plot point color
# plot heatmap using Z as color for points
par(mfrow=c(2,2))
dat$Col <- rbPal(10)[as.numeric(cut(dat$Z, breaks=10))]
plot(dat$x, dat$y, pch=20, col=dat$Col, cex=2)
text(dat$x, dat$y, round(dat$Z, 2), cex=0.8)
# show relationship between Z and seropos
plot(dat$Z, dat$seropos)
# plot heatmap using seropos as color for points
dat$Col <- rbPal(10)[as.numeric(cut(dat$seropos, breaks = 10))]
plot(dat$x, dat$y, pch = 20, col = dat$Col, cex=2)
text(dat$x, dat$y, round(dat$seropos, 2), cex=0.8)

## show linear model for seropos vs Z, and (x,y) coordinates
lmdat <- lm(dat$seropos ~ dat$Z + dat$y + dat$x)
summary(lmdat)

