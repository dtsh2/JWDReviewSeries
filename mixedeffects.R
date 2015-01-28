
######################################################
# captive bat sera anti-LBV analysis ALL POSITIVE ####
######################################################

## clear environment and load package
rm(list=ls())
library(lme4)

## get data and check first few rows
captive=read.csv("Captive_sero_Jul2011.csv", header=T) 
head(captive) 

## rename codes in Age column
levels(captive$Age) <- list(B="Neonate", A="SM", JUV="Juv",SI="SI")
captive$Age = factor(captive$Age,labels=c("Neonate","SM","Juvenile","SIM"))
#is.factor(captive$Age)

summary(captive)

## omit NAs
captive1=na.omit(captive)

## from now on, work within the captive1 data frame
attach(captive1)
#names(captive1)

## various linear models
# linear model of logTitre against Days, with different fits for each Age group
per_group_model <- lmList(LogTitre ~ Days|Age, 
                          data=captive1)
# linear model of logtitre against Days, Age, and the interaction between Days and Age
overall_model <- lm(LogTitre ~ -1 + Days*Age, 
                    data=captive1)
# linear model of logtitre against days, with random effects from ID (only for Neonates)
# get coefficients from this model
per_group_model_ranef <- lmer(LogTitre ~ Days + (Days | ID), 
                              data=captive1, subset=Age=="Neonate")
pg_coef <- coef(per_group_model_ranef)$ID

## plot logTitre against Days, draw lines defined by above coefficients
plot(jitter(LogTitre) ~ Days, data=subset(captive1, Age=="Neonate"), col=ID)
for (i in 1:13)
{
  abline(pg_coef[i,1], pg_coef[i,2], col=rownames(pg_coef)[i])
}

detach(captive1)
################
### Simulated clustered population

## define parameters
sd_ind   <- 0.1
sd_grp   <- 0.5
ovr_mean <- 5
num_ind  <- 10
meas_per_ind <- 10
ind_mean <- rnorm(num_ind, mean=ovr_mean, sd=sd_grp)

## create matrix with values
ind_val  <- matrix(0,num_ind,meas_per_ind)
for (i in 1:num_ind)
{
  ind_val[i,] <- rnorm(meas_per_ind, mean=ind_mean[i], sd=sd_ind)
}

## turn matrix into a data.frame and sort by id
x <- data.frame(id = rep(1:num_ind,meas_per_ind), grp = rep(rep(c("A","B"), each=num_ind/2), meas_per_ind), val=as.vector(ind_val))
x <- x[order(x$id),]

## various linear models
#linear model of values against group with random effects from id
summary(lmer(val ~ grp + (1|id), data=x))
# linear model of values against group
summary(lm(val ~ grp, data=x))

