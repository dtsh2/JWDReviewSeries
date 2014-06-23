
######################################################
# captive bat sera anti-LBV analysis ALL POSITIVE ####
######################################################

rm(list=ls())
setwd("~/Cambridge/CSU 2012/captive_serology")
## get data
captive=read.csv("Captive_sero_Jul2011.csv", header=T)

head(captive)

## libraries to be used
library(nlme)
library(lattice)
library(LearnBayes)
library(reshape)
library(epicalc)
library(plyr)
library(lme4)

levels(captive$Age) <- list(B="Neonate", A="SM", JUV="Juv",SI="SI")

summary(captive)
is.factor(captive$Age)
#b A JUV SI
captive$Age = factor(captive$Age,labels=c("Neonate","SM","Juvenile","SIM"))


library(ggplot2)

f=ggplot(captive, aes(Days, LogTitre))
theme_set(theme_gray(base_size = 18))
(f1=f+geom_line(aes(group=ID
                    , color=Age)))
f1+facet_grid(Age~Sex)+ylab(expression(paste
                                       (Log[2],"Titer")))
(f2=f+geom_line(aes(group=ID)))
f2+facet_grid(Age~Sex)+ylab(expression(paste
                                       (Log[2],"Titer")))
# color by ID

(f3=f+geom_line(aes(group=ID, color=ID)))
f3+facet_grid(Age~Sex)+ylab(expression(paste
                                       (Log[2],"Titer")))

# raw titre data

f=ggplot(captive, aes(Days, Titre))
theme_set(theme_gray(base_size = 18))
(f1=f+geom_line(aes(group=ID
                    , color=Age)))
f1+facet_grid(Age~Sex)+ylab(expression(paste
                                       ("Titer")))
(f2=f+geom_line(aes(group=ID)))
f2+facet_grid(Age~Sex)+ylab(expression(paste
                                       ("Titer")))
# color by ID

(f3=f+geom_line(aes(group=ID, color=ID)))
f3+facet_grid(Age~Sex)+ylab(expression(paste
                                       ("Titer")))

## omit NAs
captive1=na.omit(captive)
attach(captive1)
names(captive1)
## check data

xyplot(Titre~Days|ID, data=captive1,type="l",strip=F)
max(captive1$Titre)
captive1[Titre==2048,]
captive1[ID==42,]

xyplot(Titre~Days|Age, data=captive1, type="l", strip=T)#,groups=Age)

# remove ID42

captive42<-captive1[!captive1$ID==42,]
xyplot(Titre~Days|Age, data=captive42, type="l", strip=T)#,groups=Age)
xyplot(Titre~Days|ID, data=captive42,type="l",strip=F)

# plot without 42
f=ggplot(captive42, aes(Days, Titre))
theme_set(theme_gray(base_size = 18))
(f1=f+geom_line(aes(group=ID
                    , color=Age)))
f1+facet_grid(Age~Sex)+ylab(expression(paste
                                       ("Titer")))
(f2=f+geom_line(aes(group=ID)))
f2+facet_grid(Age~Sex)+ylab(expression(paste
                                       ("Titer")))
# color by ID

(f3=f+geom_line(aes(group=ID, color=ID)))
f3+facet_grid(Age~Sex)+ylab(expression(paste
                                       ("Titer")))

## but try to find coef for adults, so try...

test.listA=lmList(LogTitre~Days|ID, subset=Age=="SM", data=captive1)
#plot(intervals(test.listA))
coef.A=coef(test.listA)
coef.A
coef.A=na.omit(coef.A)
meanA=mean(coef.A[,2],)
meanA
sdA=sd(coef.A[,2])
sdA

test.listS=lmList(LogTitre~Days|ID, subset=Age=="SIM", data=captive1)
#plot(intervals(test.listS))
coef.S=coef(test.listS)
coef.S=na.omit(coef.S)
meanS=mean(coef.S[,2])
meanS
sdS=sd(coef.S[,2])
sdS

test.listN=lmList(LogTitre~Days|ID, subset=Age=="Neonate", data=captive1)
#plot(intervals(test.listN))
coef.N=coef(test.listN)
coef.N=na.omit(coef.N)
meanN=mean(coef.N[,2])
meanN
sdN=sd(coef.N[,2])
sdN

test.listJ=lmList(LogTitre~Days|ID, subset=Age=="Juvenile", data=captive1)
#plot(intervals(test.listJ))
coef.J=coef(test.listJ)
coef.J=na.omit(coef.J)
meanJ=mean(coef.J[,2])
meanJ
sdJ=sd(coef.J[,2])
sdA

#########################33
## no individual effects
## but try to find coef for adults, so try...

Ntest.listA=lmList(LogTitre~Days|1, subset=Age=="SM", data=captive1)
#plot(intervals(test.listA))
Ncoef.A=coef(Ntest.listA)
Ncoef.A
Ncoef.A=na.omit(Ncoef.A)
NmeanA=mean(Ncoef.A[,2],)
NmeanA

Ntest.listS=lmList(LogTitre~Days|1, subset=Age=="SIM", data=captive1)
#plot(intervals(test.listS))
Ncoef.S=coef(Ntest.listS)
Ncoef.S=na.omit(Ncoef.S)
NmeanS=mean(Ncoef.S[,2])
NmeanS

Ntest.listN=lmList(LogTitre~Days|1, subset=Age=="Neonate", data=captive1)
#plot(intervals(test.listN))
Ncoef.N=coef(Ntest.listN)
Ncoef.N=na.omit(Ncoef.N)
NmeanN=mean(Ncoef.N[,2])
NmeanN

Ntest.listJ=lmList(LogTitre~Days|1, subset=Age=="Juvenile", data=captive1)
#plot(intervals(test.listJ))
Ncoef.J=coef(Ntest.listJ)
Ncoef.J=na.omit(Ncoef.J)
NmeanJ=mean(Ncoef.J[,2])
NmeanJ

## this appears to work, although all assume 1 ind zero
## compare...
ylab=(expression(paste
                 (Log[2],"Titer")))
ylabc=(expression(paste
                  ("Change (", Log[2],"Titer)")))
#par(mfrow=c(1,2))
#boxplot(coef.N[,1],coef.J[,1],coef.S[,1],coef.A[,1],main="Intercept", names=c("Neo","Juv","SI","SM"))
#boxplot(coef.N[,2],coef.J[,2],coef.S[,2],coef.A[,2],main="Coefficients", names=c("Neo","Juv","SI","SM"))

par(mfrow=c(2,1))
boxplot(coef.N[,1],coef.J[,1],coef.S[,1],coef.A[,1],Ncoef.N[,1],Ncoef.J[,1],Ncoef.S[,1],Ncoef.A[,1],
        main="Initial Titer", names=c("Neo","Juv","SI","SM","Neo","Juv","SI","SM"),
        ylab=ylab,xlab="Age category",col=c(rep("lightgrey",4),rep("red",4)))
boxplot(coef.N[,2],coef.J[,2],coef.S[,2],coef.A[,2],Ncoef.N[,2],Ncoef.J[,2],Ncoef.S[,2],Ncoef.A[,2],
        main="Decay rate", names=c("Neo","Juv","SI","SM","Neo","Juv","SI","SM"),
        ylab=ylabc,xlab="Age category",col="lightgrey")

## it works!!!

## histograms of coeff
par(mfrow=c(2,2))
res.n<-signif(coef.N[,2],2)
hist(res.n, main="Neonates",col="grey",xlab="Decay rate")
abline(v=mean(res.n),lwd=2)
abline(v=Ncoef.N[,2],col="red",lty=2,lwd=2)
res.j<-signif(coef.J[,2],2)
hist(res.j, main="Juveniles",col="grey",xlab="Decay rate")
abline(v=mean(res.j),lwd=2)
abline(v=Ncoef.J[,2],col="red",lty=2,lwd=2)
res.s<-signif(coef.S[,2],2)
hist(res.s, main="Sexually immature",col="grey",xlab="Decay rate")
abline(v=mean(res.s),lwd=2)
abline(v=Ncoef.S[,2],col="red",lty=2,lwd=2)
res.a<-signif(coef.A[,2],2)
hist(res.a, main="Sexually mature",col="grey",xlab="Decay rate")
abline(v=mean(res.a),lwd=2)
abline(v=Ncoef.A[,2],col="red",lty=2,lwd=2)

xx<-0:365
int<-3
coeffs<-cbind(mean(res.n),Ncoef.N[,2],
              mean(res.j),Ncoef.J[,2],
              mean(res.s),Ncoef.S[,2],
              mean(res.a),Ncoef.A[,2])

y<-matrix(NA,ncol=8,nrow=366)

for (i in 1:length(coeffs[1,])) {
  y[,i]<-int+coeffs[1,i]*xx
}

## histograms of coeff
par(mfrow=c(2,4))
res.n<-signif(coef.N[,2],2)
hist(res.n, main="Neonates",col="grey",xlab="Decay rate")
abline(v=mean(res.n),lwd=2)
abline(v=Ncoef.N[,2],col="red",lty=2,lwd=2)
matplot(y[,1:2],type="l",lwd=2,main="Neonates",ylab="Titre",xlab="Days",ylim=c(0,3))

res.j<-signif(coef.J[,2],2)
hist(res.j, main="Juveniles",col="grey",xlab="Decay rate")
abline(v=mean(res.j),lwd=2)
abline(v=Ncoef.J[,2],col="red",lty=2,lwd=2)
matplot(y[,3:4],type="l",lwd=2,main="Juveniles",ylab="Titre",xlab="Days",ylim=c(0,3))

res.s<-signif(coef.S[,2],2)
hist(res.s, main="Sexually immature",col="grey",xlab="Decay rate")
abline(v=mean(res.s),lwd=2)
abline(v=Ncoef.S[,2],col="red",lty=2,lwd=2)
matplot(y[,5:6],type="l",lwd=2,main="Sexually immature",ylab="Titre",xlab="Days",ylim=c(0,3))

res.a<-signif(coef.A[,2],2)
hist(res.a, main="Sexually mature",col="grey",xlab="Decay rate")
abline(v=mean(res.a),lwd=2)
abline(v=Ncoef.A[,2],col="red",lty=2,lwd=2)
matplot(y[,7:8],type="l",lwd=2,main="Sexually mature",ylab="Titre",xlab="Days",ylim=c(0,3))

