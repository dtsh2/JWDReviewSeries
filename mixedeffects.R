
######################################################
# captive bat sera anti-LBV analysis ALL POSITIVE ####
######################################################

rm(list=ls())
setwd("~/Cambridge/CSU 2012/captive_serology")
## get data
captive=read.csv("Captive_sero_Jul2011.csv", header=T)

head(captive)

library(lme4)

levels(captive$Age) <- list(B="Neonate", A="SM", JUV="Juv",SI="SI")

summary(captive)
is.factor(captive$Age)
#b A JUV SI
captive$Age = factor(captive$Age,labels=c("Neonate","SM","Juvenile","SIM"))

## omit NAs
captive1=na.omit(captive)
attach(captive1)
names(captive1)
## check data

per_group_model <- lmList(LogTitre ~ Days|Age, data=captive1)

overall_model <- lm(LogTitre ~ -1 + Days*Age, data=captive1)

per_group_model_ranef <- lmer(LogTitre ~ Days + (Days | ID), data=captive1, subset=Age=="Neonate")
pg_coef <- coef(per_group_model_ranef)$ID
plot(jitter(LogTitre) ~ Days, data=subset(captive1, Age=="Neonate"), col=ID)
for (i in 1:13)
{
  abline(pg_coef[i,1], pg_coef[i,2], col=rownames(pg_coef)[i])
}

################
#
# some stuff to simulate from a clustered population
#
#
sd_ind   <- 0.1
sd_grp   <- 0.5
ovr_mean <- 5
num_ind  <- 10
meas_per_ind <- 10
ind_mean <- rnorm(num_ind, mean=ovr_mean, sd=sd_grp)
ind_val  <- matrix(0,num_ind,meas_per_ind)
for (i in 1:num_ind)
{
  ind_val[i,] <- rnorm(meas_per_ind, mean=ind_mean[i], sd=sd_ind)
}

x <- data.frame(id = rep(1:num_ind,meas_per_ind), grp = rep(rep(c("A","B"), each=num_ind/2), meas_per_ind), val=as.vector(ind_val))
x <- x[order(x$id),]
summary(lmer(val ~ grp + (1|id), data=x))
summary(lm(val ~ grp, data=x))

