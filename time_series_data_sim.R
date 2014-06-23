rm(list=ls())

library(waveslim)
library(dplR)
library(xts)
library(wavelets)
library(wavethresh)
library(gplots)

## simulate cycles with frequency change and see what happens
##########################
##
## Poisson
##

TIME<-1:1000
w<-rpois(1000,lambda=20)
plot(1:1000,w)
test<-morlet(w, x1 = seq_along(w), p2 = NULL, dj = 0.25, siglvl = 0.95)
wavelet.plot(test)
acf(w, lag.max = length(w),plot=T,main="")

Seasonal<-3*sin(2*pi*TIME/100)+2*sin(2*pi*TIME/100)#+0.1*sin(2*pi*4*TIME/14)+0.1*cos(2*pi*4*TIME/14)
x<-Seasonal+w
plot(Seasonal)
plot(x)#,type="l")
#test<-morlet(x, x1 = seq_along(x), p2 = NULL, dj = 0.25, siglvl = 0.95)
#wavelet.plot(test)

kappa=20; s=25; omega=10;phi=0; t=TIME

peaks<-(kappa*(1/sqrt((1/s)*pi)*exp(-((cos(pi*omega*t/1000-phi))^2)/(1/s))))
plot(peaks,type="l")
x3<-peaks+Seasonal+w
plot(x3)#,type="l")
test<-morlet(x3, x1 = seq_along(x3), p2 = NULL, dj = 0.25, siglvl = 0.95)
wavelet.plot(test)

acf(x3, lag.max = length(x3),plot=T,main="")

trend<-0.05*TIME
x4<-peaks+Seasonal+w+trend
plot(x4)#,type="l")
test<-morlet(x4, x1 = seq_along(x4), p2 = NULL, dj = 0.25, siglvl = 0.95)
wavelet.plot(test)

plot(1:1000,x,
     #type="l",
     ylim=c(0,350),axes="n")
lines(1:1000,Seasonal+80)
lines(1:1000,peaks+120)
lines(1:1000,trend+200)
points(1:1000,x4+220)
mtext(2,text="Signal",line=2,cex=1)
#mtext(2,text="Random  Daily Periodic  Trend Combined",line=1,at=c(25),cex=0.8)

acf(x4, lag.max = length(x4),plot=T,main="")

max(x4)
min(x4)

spos<-x4
sneg<-300-x4
plot(spos/(spos+sneg),type="l",ylab="Proportion seropositive")
plot(100*spos/(spos+sneg),type="l",ylab="Seroprevalence (%)")

###########
## data and confidence intervals

sp<-spos/(spos+sneg)
tot<-spos+sneg
res.new<-rbind(round(spos,0),tot,sp)
rownames(res.new)<-c("P","N","Sp")
res.new<-t(res.new)
class(res.new)

res.u<-matrix(NA,ncol=1000,nrow=2)
for (ii in 1:1000){
  res.u[1,ii]<-binom.test(res.new[ii],res.new[ii+1000])$conf.int[1]
  res.u[2,ii]<-binom.test(res.new[ii],res.new[ii+1000])$conf.int[2]
}

res.u<-t(res.u)
res.all<-cbind(res.new,res.u)
dim(res.all)
colnames(res.all)<-c("P","N","Sp","L","U")

res.all<-as.data.frame(res.all)

## each sample get CIs

barplot(res.all$Sp,ylim=c(0,1))
barplot2(res.all$Sp, plot.ci=TRUE, ci.l=res.all$L, ci.u=res.all$U,
         ylim=c(0,1),main="",ci.lty=2)#,names.arg=c(1:20))

barplot2(res.all$Sp[1:10], plot.ci=TRUE, ci.l=res.all$L[1:10], ci.u=res.all$U[1:10],
         ylim=c(0,1),main="",ci.lty=2)#,names.arg=c(1:20))
####
barplot2(res.all$Sp[c(1,100,150,300,1000)], plot.ci=TRUE, ci.l=res.all$L[c(1,100,150,300,1000)], ci.u=res.all$U[c(1,100,150,300,1000)],
         ylim=c(0,1),main="",ci.lty=2)#,names.arg=c(1:20))

## to here

res.P<-res.all$P[c(1,100,150,300,1000)]
res.T<-res.all$N[c(1,100,150,300,1000)]
res.N<-res.T-res.P
## chi squared test
prop.test(res.P,res.N)
res.data<-cbind(res.P,res.N)
chisq.test(res.data)
res.SP<-res.P/res.T
t<-1:5
lmres<-lm(res.SP~t)
abline(lmres,col="red",
       lty=2)
####

#acf(res.all$Sp, lag.max = length(res.all$Sp),plot=T,main="")
#acf(res.SP)

par(mfrow=c(2,4))
acf(res.all$Sp[seq(from=1,to=1000,by=100)], lag.max = length(res.all$Sp[seq(from=1,to=1000,by=100)]),plot=T,main="")
barplot2(res.all$Sp[seq(from=1,to=1000,by=100)], plot.ci=TRUE, ci.l=res.all$L[seq(from=1,to=1000,by=100)], ci.u=res.all$U[seq(from=1,to=1000,by=100)],
         ylim=c(0,1),main="",ci.lty=2)#,names.arg=c(1:20))
lmres<-lm(res.all$Sp[seq(from=1,to=1000,by=100)]~seq(from=1,to=1000,by=100))
abline(lmres,col="red",
       lty=2,lwd=2)
pt<-prop.test(res.all$P[seq(from=1,to=1000,by=100)],res.all$N[seq(from=1,to=1000,by=100)])
text(expression(paste(chi^2, "=")),x=length(seq(from=1,to=1000,by=100))/2,y=0.8)
text(round(pt$statistic,0),x=length(seq(from=1,to=1000,by=100))/2,y=0.7)
text(expression(paste("p-value =")),x=length(seq(from=1,to=1000,by=100))/2,y=0.6)
text(signif(pt$p.value,2),x=length(seq(from=1,to=1000,by=100))/2,y=0.5)

acf(res.all$Sp[seq(from=1,to=1000,by=10)], lag.max = length(res.all$Sp[seq(from=1,to=1000,by=10)]),plot=T,main="")
barplot2(res.all$Sp[seq(from=1,to=1000,by=10)], plot.ci=TRUE, ci.l=res.all$L[seq(from=1,to=1000,by=10)], ci.u=res.all$U[seq(from=1,to=1000,by=10)],
         ylim=c(0,1),main="",ci.lty=2)#,names.arg=c(1:20))
lmres<-lm(res.all$Sp[seq(from=1,to=1000,by=10)]~seq(from=1,to=1000,by=10))
abline(lmres,col="red",
       lty=2,lwd=2)
pt<-prop.test(res.all$P[seq(from=1,to=1000,by=10)],res.all$N[seq(from=1,to=1000,by=10)])
text(expression(paste(chi^2, "=")),x=length(seq(from=1,to=1000,by=10))/2,y=0.8)
text(round(pt$statistic,0),x=length(seq(from=1,to=1000,by=10))/2,y=0.7)
text(expression(paste("p-value =")),x=length(seq(from=1,to=1000,by=10))/2,y=0.6)
text(signif(pt$p.value,2),x=length(seq(from=1,to=1000,by=10))/2,y=0.5)

acf(res.all$Sp[seq(from=1,to=100,by=10)], lag.max = length(res.all$Sp[seq(from=1,to=100,by=10)]),plot=T,main="")
barplot2(res.all$Sp[seq(from=1,to=100,by=10)], plot.ci=TRUE, ci.l=res.all$L[seq(from=1,to=100,by=10)], ci.u=res.all$U[seq(from=1,to=100,by=10)],
         ylim=c(0,1),main="",ci.lty=2)#,names.arg=c(1:20))
lmres<-lm(res.all$Sp[seq(from=1,to=100,by=10)]~seq(from=1,to=100,by=10))
abline(lmres,col="red",
       lty=2,lwd=2)
pt<-prop.test(res.all$P[seq(from=1,to=100,by=10)],res.all$N[seq(from=1,to=100,by=10)])
text(expression(paste(chi^2, "=")),x=length(seq(from=1,to=100,by=10))/2,y=0.8)
text(round(pt$statistic,0),x=length(seq(from=1,to=100,by=10))/2,y=0.7)
text(expression(paste("p-value =")),x=length(seq(from=1,to=100,by=10))/2,y=0.6)
text(signif(pt$p.value,2),x=length(seq(from=1,to=100,by=10))/2,y=0.5)

acf(res.all$Sp[seq(from=1,to=100,by=5)], lag.max = length(res.all$Sp[seq(from=1,to=100,by=5)]),plot=T,main="")
barplot2(res.all$Sp[seq(from=1,to=100,by=5)], plot.ci=TRUE, ci.l=res.all$L[seq(from=1,to=100,by=5)], ci.u=res.all$U[seq(from=1,to=100,by=5)],
         ylim=c(0,1),main="",ci.lty=2)#,names.arg=c(1:20))
lmres<-lm(res.all$Sp[seq(from=1,to=100,by=5)]~seq(from=1,to=100,by=5))
abline(lmres,col="red",
       lty=2,lwd=2)
pt<-prop.test(res.all$P[seq(from=1,to=100,by=5)],res.all$N[seq(from=1,to=100,by=5)])
text(expression(paste(chi^2, "=")),x=length(seq(from=1,to=100,by=5))/2,y=0.8)
text(round(pt$statistic,0),x=length(seq(from=1,to=100,by=5))/2,y=0.7)
text(expression(paste("p-value =")),x=length(seq(from=1,to=100,by=5))/2,y=0.6)
text(signif(pt$p.value,2),x=length(seq(from=1,to=100,by=5))/2,y=0.5)

