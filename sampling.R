## adjust p values

a <- 5 # mean - use the same one
s <- 2 # sd
num_tests <- 20
n <- 20 # sample size
xsamples <- matrix(rnorm(num_tests * n, mean=a, sd=s), ncol=num_tests)

par(mfrow=c(2,1))
# boxplot of data against population mean
boxplot(xsamples, xlim=c(0,num_tests+1))
abline(h=a, col="red")

# calculate sample means and ci's
xbar <- apply(xsamples, 2, mean)
ci   <- apply(xsamples, 2, function(x) { t.test(x)$conf.int })

# plot means, ci's and p-values
plot(xbar, xlim=c(0,num_tests+1), ylim=range(xsamples), pch=19)
abline(h=a, col="red")
for (i in 1:num_tests)
  lines(c(i,i), ci[,i])
pvals <- apply(xsamples, 2, function(x) { t.test(x, mu=a)$p.value })
text(1:num_tests, max(ci)+0.5, round(pvals, 2))

padj<-p.adjust(pvals, method = "bonferroni", n = length(pvals))
text(1:num_tests, min(ci)-2, round(padj, 2))

# plot p-values and adjusted p-values
plot(pvals)
abline(h=0.05)

plot(padj)
abline(h=0.05)

