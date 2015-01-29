rm(list=ls())

### adjustment of p values

## create matrix of random samples in a normal distribution, repeated tests on samples
a <- 5 # mean - use the same one
s <- 2 # sd
num_tests <- 20
n <- 20 # sample size
xsamples <- matrix(rnorm(num_tests * n, mean=a, sd=s), ncol=num_tests)

## format plots
par(mfrow=c(1,1))
par(oma=c(0.1,1,0.1,0.2))
par(mar=c(2,5,2,2))

## create .tiff with boxplots of data, draw line for population mean
tiff("data.tiff", width=8, height=8, units='in', res=300, compression="lzw")
boxplot(xsamples, xlim=c(0,num_tests+1), ylab="data", cex.lab=1.5, xlab="sample number")
abline(h=a, col="red")

dev.off() 

## calculate sample means and confidence intervals
xbar <- apply(xsamples, 2, mean)
ci   <- apply(xsamples, 2, function(x) { t.test(x)$conf.int })

## create .tiff with means, confidence intervals, p-values, and Bonferri adjusted p-values 
# means & ci's
tiff("mean.tiff", width=8, height=8, units='in', res=300, compression="lzw")
plot(xbar, xlim=c(0,num_tests+1), ylim=range(xsamples), pch=19,
     ylab="mean (95% confidence intervals)", cex.lab=1.5, xlab="sample number")
abline(h=a, col="red")
for (i in 1:num_tests)
  lines(c(i,i), ci[,i])

# calculate p-values and add to plot
pvals <- apply(xsamples, 2, function(x) { t.test(x, mu=a)$p.value })
text(1:num_tests, max(ci)+2, round(pvals, 2), srt=90, cex=1.5)
text(0, max(ci)+2, "p-values", srt=90, cex=1.5)

# calculate Bonferroni adjusted p-values and add to plot 
padj<-p.adjust(pvals, method = "bonferroni", n = length(pvals))
text(1:num_tests, min(ci)-2, round(padj, 2), srt=90, cex=1.5)
text(0, min(ci)-2, "adj. p-values", srt=90, cex=1.5)

dev.off()

## plot p-values and various different adjusted p-values
# create and format new .tiff
tiff("pvals.tiff", width=8, height=8, units='in', res=300, compression="lzw")
par(mfrow=c(6,1))
par(oma=c(1,1,0.1,0.2))
par(mar=c(2,5,2,2))

# plot p-values + dotted line at level of significance
plot(pvals, ylim=c(0,1), pch=16, xlab="", ylab="p-values", cex.lab=1.5)
abline(h=0.05, lty=2)

# plot adjusted p-values by Bonferroni method
plot(padj, ylim=c(0,1), pch=16, xlab="", ylab="Bonferroni", cex.lab=1.5)
abline(h=0.05, lty=2)

# plot adjusted p-values by Holm method
padj<-p.adjust(pvals, method = "holm", n = length(pvals))
plot(padj, ylim=c(0,1), pch=16, xlab="", ylab="Holm", cex.lab=1.5)
abline(h=0.05, lty=2)

# plot adjusted p-values by Hochberg method
padj<-p.adjust(pvals, method = "hochberg", n = length(pvals))
plot(padj, ylim=c(0,1), pch=16, xlab="", ylab="Hochberg", cex.lab=1.5)
abline(h=0.05, lty=2)

# plot adjusted p-values by BH method
padj<-p.adjust(pvals, method = "BH", n = length(pvals)) # same as fdr
plot(padj, ylim=c(0,1), pch=16, xlab="", ylab="BH", cex.lab=1.5)
abline(h=0.05, lty=2)

# plot adjusted p-values by BY method
padj<-p.adjust(pvals, method = "BY", n = length(pvals)) # same as fdr
plot(padj, ylim=c(0,1), pch=16, ylab="BY", cex.lab=1.5, xlab="sample number")
abline(h=0.05, lty=2)

dev.off()

##########
# Bonferroni correction ("bonferroni") 
# in which the p-values are multiplied by the number of comparisons. 
# Less conservative corrections are also included by Holm (1979) ("holm"),
# Hochberg (1988) ("hochberg"), Hommel (1988) ("hommel"), 
# Benjamini & Hochberg (1995) ("BH" or its alias "fdr"), 
# and Benjamini & Yekutieli (2001) ("BY"), respectively.
# Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B 57, 289-300.
#
# Benjamini, Y., and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29, 1165-1188.
#
# Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics 6, 65-70.
#
# Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. Biometrika 75, 383-386.
#
# Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. Biometrika 75, 800-803.
#
# Shaffer, J. P. (1995). Multiple hypothesis testing. Annual Review of Psychology 46, 561-576. (An excellent review of the area.)
#
# Sarkar, S. (1998). Some probability inequalities for ordered MTP2 random variables: a proof of Simes conjecture. Annals of Statistics 26, 494-504.
#
# Sarkar, S., and Chang, C. K. (1997). Simes' method for multiple hypothesis testing with positively dependent test statistics. Journal of the American Statistical Association 92, 1601-1608.
#
# Wright, S. P. (1992). Adjusted P-values for simultaneous inference. Biometrics 48, 1005-1013.
#
# Bonferroni, C. E., Teoria statistica delle classi e calcolo delle probabilit?, Pubblicazioni del R Istituto Superiore di Scienze Economiche e Commerciali di Firenze 1936
# Olive Jean Dunn, Estimation of the Medians for Dependent Variables, The Annals of Mathematical Statistics 1959 http://projecteuclid.org/download/pdf_1/euclid.aoms/1177706374
# Olive Jean Dunn, Multiple Comparisons Among Means, Journal of the American Statistical Association 1961 http://sci2s.ugr.es/keel/pdf/algorithm/articulo/1961-Bonferroni_Dunn-JASA.pdf
