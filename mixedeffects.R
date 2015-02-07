
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
pg_coef.Neonate <- coef(per_group_model_ranef)$ID
# linear model of average for whole Neonate group
per_group_model.Neonate <- lm(LogTitre ~ Days, data=captive1, subset=Age=="Neonate")

## plot logTitre against Days, draw linear models for each indvidual (color coded) 
## and the group as a whole (black and dashed)
par(mfrow=c(1, 1))
plot(jitter(LogTitre) ~ Days, data=subset(captive1, Age=="Neonate"), 
     col=ID, main="Neonate")
abline(per_group_model.Neonate, col='black', lwd=3, lty=2)

for (i in 1:13)
{
  abline(pg_coef.Neonate[i,1], pg_coef.Neonate[i,2], col=rownames(pg_coef.Neonate)[i])
}

## Do the same for SM
# individual lms
per_group_model_ranef <- lmer(LogTitre ~ Days + (Days | ID), 
                              data=captive1, subset=Age=="SM")
pg_coef.SM <- coef(per_group_model_ranef)$ID

# average lm
per_group_model.SM <- lm(LogTitre ~ Days, data=captive1, subset=Age=="SM")

# plot 
par(mfrow=c(1, 1))
plot(jitter(LogTitre) ~ Days, data=subset(captive1, Age=="SM"), 
     col=ID, main="SM")
abline(per_group_model.SM, col='black', lwd=3, lty=2)
for (i in 1:57)
{
  abline(pg_coef.SM[i,1], pg_coef.SM[i,2], col=rownames(pg_coef.SM)[i])
}

## Do the same for Juvenile
# individual lms
per_group_model_ranef <- lmer(LogTitre ~ Days + (Days | ID), 
                              data=captive1, subset=Age=="Juvenile")
pg_coef.Juvenile <- coef(per_group_model_ranef)$ID

# average lm
per_group_model.Juvenile <- lm(LogTitre ~ Days, data=captive1, subset=Age=="Juvenile")

# plot 
par(mfrow=c(1, 1))
plot(jitter(LogTitre) ~ Days, data=subset(captive1, Age=="Juvenile"), 
     col=ID, main="Juvenile")
abline(per_group_model.Juvenile, col='black', lwd=3, lty=2)
for (i in 1:8)
{
  abline(pg_coef.Juvenile[i,1], pg_coef.Juvenile[i,2], col=rownames(pg_coef.Juvenile)[i])
}

## Do the same for SIM
# individual lms
per_group_model_ranef <- lmer(LogTitre ~ Days + (Days | ID), 
                              data=captive1, subset=Age=="SIM")
pg_coef.SIM <- coef(per_group_model_ranef)$ID

# average lm
per_group_model.SIM <- lm(LogTitre ~ Days, data=captive1, subset=Age=="SIM")

# plot 
par(mfrow=c(1, 1))
plot(jitter(LogTitre) ~ Days, data=subset(captive1, Age=="SIM"), 
     col=ID, main="SIM")
abline(per_group_model.SIM, col='black', lwd=3, lty=2)
for (i in 1:11)
{
  abline(pg_coef.SIM[i,1], pg_coef.SIM[i,2], col=rownames(pg_coef.SIM)[i])
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
x <- data.frame(id=rep(1:num_ind, meas_per_ind), 
                grp=rep(rep(c("A", "B"), each=num_ind/2), meas_per_ind), 
                val=as.vector(ind_val))
x <- x[order(x$id),]

## various linear models
#linear model of values against group with random effects from id
summary(lmer(val ~ grp + (1|id), data=x))
# linear model of values against group
summary(lm(val ~ grp, data=x))

