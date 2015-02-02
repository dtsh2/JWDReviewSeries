## clear environment and load packages
rm(list=ls())
library(waveslim)
library(dplR)
library(xts)
library(wavelets)
library(wavethresh)
library(gplots)
library(nlme)

## simulate cycles with frequency change and see what happens
##########################
##
## Poisson
##
## format plot areas
par(mfrow=c(1,1))
par(mar=c(5,5,2,7))
## create some random data with a poisson distribution
TIME <- 1:1000
w <- rpois(1000, lambda=20)
plot(TIME, w, pch=20)

## filled contour plot of continuous Morlet wavelet transform of data
test <- morlet(w, x1=seq_along(w), p2=NULL, dj=0.25, siglvl=0.95)
wavelet.plot(test)

## autocorrelation plot
acf(w, lag.max=length(w), plot=T, main="")

## add seasonal pattern to data 
# create and plot seasonal pattern
Seasonal <- 3 * sin(2*pi*TIME/100) + 2 * sin(2*pi*TIME/100) #+ 0.1 * sin(2*pi*4*TIME/14) + 0.1 * cos(2*pi*4*TIME/14)
plot(Seasonal, type="l")

# add seasonal pattern to random data and plot
x <- Seasonal + w
plot(x, pch=20)

# Morlet wavelet transform plot
test <- morlet(x, x1=seq_along(x), p2=NULL, dj=0.25, siglvl=0.95)
wavelet.plot(test)

# autocorrelation plot
acf(x, lag.max = length(x),plot=T,main="")

## add regular peaks to data
# define peak parameters
kappa = 20 # height
s = 25 # width
omega = 10 # frequency
phi = 0 # offset
t = TIME

# create and plot peaks
peaks <- kappa * (1 / sqrt((1/s)*pi) * exp(-((cos(pi*omega*t/1000 - phi))^2) / (1/s)))
plot(peaks, type="l")

# add peaks to above data and plot
x2 <- peaks + x
plot(x2, pch=20)

# Morlet wavelet transform plot
test <- morlet(x2, x1=seq_along(x2), p2=NULL, dj=0.25, siglvl=0.95)
wavelet.plot(test)

# autocorrelation plot
acf(x2, lag.max = length(x2),plot=T,main="")

## add an upward trend to data
# create and plot trend
trend <- 0.05 * TIME
plot(trend, type="l")

# add peaks to above data and plot
x3 <- trend + x2
plot(x3, pch=20)

#Morlet wavelet transform plot
test<-morlet(x3, x1 = seq_along(x3), p2 = NULL, dj = 0.25, siglvl = 0.95)
wavelet.plot(test)

# autocorrelation plot
acf(x3, lag.max = length(x2),plot=T,main="")

## Show all the different components on one graph
plot(TIME, w, ylim=c(0,350), pch=20, axes='n')
lines(TIME, Seasonal+80)
lines(TIME, peaks+120)
lines(TIME, trend+200)
points(TIME, x3+220, pch=20)
axis(side=4, at=c(20, 80, 150, 230, 300), 
     labels=c("Random", "Daily", "Periodic", "Trend", "Combined"), 
     tick=FALSE, las=1)
mtext(side=1, "Time", line=1)
mtext(side=2, "Signal", line=1)

## have above data represent # of seropositives (spos) in a population of 300 over time
tot <- 300
spos <- x3
sneg <- tot - x3
sp <- spos / tot

# plot proportion seropositive (sp) and seroprevalence %
plot(sp, type="l", ylab="Proportion seropositive")
plot(100 * sp, type="l", ylab="Seroprevalence (%)")

###########
### data and confidence intervals

## create matrix of total seropositives, total population, and seroprevalence
res.new <- cbind(round(spos,0), tot, sp)
colnames(res.new)<-c("P","N","Sp")

## create matrix with bounds for confidence intervals based on binomial test for res.new
res.u <- matrix(NA, ncol=2, nrow=1000)
for (ii in 1:1000){
  res.u[ii,1] <- binom.test(res.new[ii], res.new[ii+1000])$conf.int[1]
  res.u[ii,2] <- binom.test(res.new[ii], res.new[ii+1000])$conf.int[2]
}

## combines res.new and res.u and turns them into a data.frame
res.all <- cbind(res.new, res.u)
colnames(res.all) <- c("P", "N", "Sp", "L", "U")
res.all <- as.data.frame(res.all)

## barplots with confidence intervals showing various sampling strategies
# barplot of seroprevalence, 2nd one includes error bars
barplot(res.all$Sp, ylim=c(0,1))
barplot2(res.all$Sp, plot.ci=TRUE, ci.l=res.all$L, ci.u=res.all$U,
         ylim=c(0,1), main="", ci.lty=2, names.arg=c(1:1000))
# barplot of only the first 10 samples
barplot2(res.all$Sp[1:10], plot.ci=TRUE, ci.l=res.all$L[1:10], ci.u=res.all$U[1:10],
         ylim=c(0,1),main="",ci.lty=2, names.arg=c(1:10))
# barplot of samples spaced throughout
time.spacing <- c(1,100,150,300,1000)
barplot2(res.all$Sp[time.spacing], plot.ci=TRUE, 
         ci.l=res.all$L[time.spacing], ci.u=res.all$U[time.spacing],
         ylim=c(0,1), main="",ci.lty=2, names.arg=time.spacing)

## store data corresponding to time.spacing in new vectors 
# (P = seropositive, T = total, N = seronegative, SP = seroprevalence)
res.P <- res.all$P[time.spacing]
res.T <- res.all$N[time.spacing]
res.N <- res.T - res.P
res.SP <- res.P / res.T

## adds linear regression for seropositive proportion to plot
t <- 1:5
lmres <- lm(res.SP ~ t)
abline(lmres, col="red", lty=2)

## 2 methods for getting chi squared test for seropositives, total
prop.test(res.P, res.T) #1
res.data <- cbind(res.P, res.T) #2
chisq.test(res.data)

##############
## autocorrelation function for all data vs. time.spacing data
acf(res.all$Sp, lag.max = length(res.all$Sp), plot=T, main="")
acf(res.SP)

### show how autocorrelation and linear regression change with various sampling strategies
## format graph display for following section
par(mfrow=c(2,4))

## sample every 100 
every.hundred <- seq(from=1, to=1000, by=100)
# autocorrelation
acf(res.all$Sp[every.hundred], lag.max = length(res.all$Sp[every.hundred]), 
    plot=T, main="Autocorrelation every 100")
# make barplot of seropos data counting by 100s
barplot2(res.all$Sp[every.hundred], plot.ci=TRUE, 
         ci.l=res.all$L[every.hundred], ci.u=res.all$U[every.hundred],
         ylim=c(0,1), main="Sample every 100", ci.lty=2, names.arg=every.hundred)
# add linear regression line
lmres<-lm(res.all$Sp[every.hundred] ~ c(1:10))
abline(lmres, col="red", lty=2, lwd=2) 
# write chi^2 and p-value above graph
pt <- prop.test(res.all$P[every.hundred], res.all$N[every.hundred])
text(expression(paste(chi^2, "=")), x=length(every.hundred)/2, y=0.8)
text(round(pt$statistic, 0), x=length(every.hundred)/2, y=0.7)
text(expression(paste("p-value =")), x=length(every.hundred)/2, y=0.6)
text(signif(pt$p.value, 2), x=length(every.hundred)/2, y=0.5)

## sample every 10
every.ten <- seq(from=1, to=1000, by=10)
# autocorrelation
acf(res.all$Sp[every.ten], lag.max = length(res.all$Sp[every.ten]), 
    plot=T, main="Autocorrelation every 10")
# make barplot of seropos data counting by 100s
barplot2(res.all$Sp[every.ten], plot.ci=TRUE, 
         ci.l=res.all$L[every.ten], ci.u=res.all$U[every.ten],
         ylim=c(0,1), main="Sample every 10", ci.lty=2, names.arg=every.ten)
# add linear regression line
lmres<-lm(res.all$Sp[every.ten] ~ c(1:100))
abline(lmres, col="red", lty=2, lwd=2) 
# write chi^2 and p-value above graph
pt <- prop.test(res.all$P[every.ten], res.all$N[every.ten])
text(expression(paste(chi^2, "=")), x=length(every.ten)/2, y=0.8)
text(round(pt$statistic, 0), x=length(every.ten)/2, y=0.7)
text(expression(paste("p-value =")), x=length(every.ten)/2, y=0.6)
text(signif(pt$p.value, 2), x=length(every.ten)/2, y=0.5)

## sample every 10 within the first 100
hundred.by.ten = seq(from=1, to=100, by=10)
# autocorrelation
acf(res.all$Sp[hundred.by.ten], lag.max = length(res.all$Sp[hundred.by.ten]), 
    plot=T, main="Autocorrelation every 10 until 100")
# make barplot of seropos data counting by 100s
barplot2(res.all$Sp[hundred.by.ten], plot.ci=TRUE, 
         ci.l=res.all$L[hundred.by.ten], ci.u=res.all$U[hundred.by.ten],
         ylim=c(0,1), main="Sample every 10 until 100", ci.lty=2, names.arg=hundred.by.ten)
# add linear regression line
lmres<-lm(res.all$Sp[hundred.by.ten] ~ c(1:10))
abline(lmres, col="red", lty=2, lwd=2) 
# write chi^2 and p-value above graph
pt <- prop.test(res.all$P[hundred.by.ten], res.all$N[hundred.by.ten])
text(expression(paste(chi^2, "=")), x=length(hundred.by.ten)/2, y=0.8)
text(round(pt$statistic, 0), x=length(hundred.by.ten)/2, y=0.7)
text(expression(paste("p-value =")), x=length(hundred.by.ten)/2, y=0.6)
text(signif(pt$p.value, 2), x=length(hundred.by.ten)/2, y=0.5)

## sample every 5 within the first 100
hundred.by.five = seq(from=1, to=100, by=5)
# autocorrelation
acf(res.all$Sp[hundred.by.five], lag.max = length(res.all$Sp[hundred.by.five]), 
    plot=T, main="Autocorrelation every 5 until 100")
# make barplot of seropos data counting by 100s
barplot2(res.all$Sp[hundred.by.five], plot.ci=TRUE, 
         ci.l=res.all$L[hundred.by.five], ci.u=res.all$U[hundred.by.five],
         ylim=c(0,1), main="Sample every 5 until 100", ci.lty=2, names.arg=hundred.by.five)
# add linear regression line
lmres<-lm(res.all$Sp[hundred.by.five] ~ c(1:20))
abline(lmres, col="red", lty=2, lwd=2) 
# write chi^2 and p-value above graph
pt <- prop.test(res.all$P[hundred.by.five], res.all$N[hundred.by.five])
text(expression(paste(chi^2, "=")), x=length(hundred.by.five)/2, y=0.8)
text(round(pt$statistic, 0), x=length(hundred.by.five)/2, y=0.7)
text(expression(paste("p-value =")), x=length(hundred.by.five)/2, y=0.6)
text(signif(pt$p.value, 2), x=length(hundred.by.five)/2, y=0.5)

###########
## a few different chi squared methods
# save seropositive and seronegative data for every 5 up to 100
pos <- res.all$P[hundred.by.five]
neg <- res.all$N[hundred.by.five] - pos

# use chi-sq via prop.test
pt <- prop.test(pos, pos+neg)

# use GLM 
t   <- 1:20
mat <- cbind(pos, neg)
mod.glm <- glm(mat ~ as.factor(t), family="binomial")
anova(mod.glm, test="Chisq")

# ideally include auto-correlation, but this is too hard :(
#mod.glm <- glmer(mat ~ as.factor(t), family="binomial")
#anova(mod.glm, test="Chisq")
#glmmPQL(mat ~ as.factor(t), family="binomial", correlation=corAR1())

###########
## Plot autocorrelation and 2 different linear models
## lmres2 accounts for autocorrelation
## NOTE here that we're assuming constant sample size across time.
## with different sample sizes we'd want to weight each datapoint by the (inverse of?) sample size

# sample every 100, take 2 different types of linear model
t <- 1:10
seroprev <- res.all$Sp[every.hundred]
lmres<-lm(seroprev ~ t)
lmres2<-gls(seroprev~t, correlation=corAR1())

# reformat plot area
par(mfrow=c(1,1))

# autocorrelation 
acf(seroprev, lag.max=length(seroprev), plot=T, main="Autocorrelation every 100")

# plot + add linear models
barplot2(seroprev, plot.ci=TRUE, 
         ci.l=res.all$L[every.hundred], ci.u=res.all$U[every.hundred],
         ylim=c(0,1),main="Sample every 100",ci.lty=2, names.arg=every.hundred)
abline(lmres, col="red", lty=2,lwd=2)
abline(lmres2, col="blue")

# add chi^2 and p-value text
pt <- prop.test(res.all$P[every.hundred], res.all$N[every.hundred])
text(expression(paste(chi^2, "=")), x=length(every.hundred)/2, y=0.8)
text(round(pt$statistic, 0), x=length(every.hundred)/2, y=0.7)
text(expression(paste("p-value =")), x=length(every.hundred)/2, y=0.6)
text(signif(pt$p.value, 2), x=length(every.hundred)/2, y=0.5)


