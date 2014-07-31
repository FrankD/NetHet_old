#########################################
##This an example of how to use DiffNet##
#########################################


#####clear workspace
rm(list=ls())
set.seed(1)

#####generate data
p <- 100
n <- 60
x1 <- matrix(rnorm(n*p),n,p)
x2 <- x3 <- matrix(rnorm(n*p),n,p)
beta1 <- c(rep(1,2),rep(0,p-2))
beta2 <- rep(0,p)
y1 <- x1%*%as.matrix(beta1)+rnorm(n,sd=1)
y2 <- x2%*%as.matrix(beta2)+rnorm(n,sd=1)
y3 <- x3%*%as.matrix(beta1)+rnorm(n,sd=1)

#####Diffregr (asymptotic p-values)
##alternative case
fit1 <- diffregr_multisplit(y1,y2,x1,x2,b.splits=10)
hist(fit1$pval.onesided,breaks=10)
fit1$aggpval.onesided
fit1$medpval.onesided
##null case
fit2 <- diffregr_multisplit(y1,y3,x1,x2,b.splits=10)
hist(fit2$pval.onesided,breaks=10)
fit2$aggpval.onesided
fit2$medpval.onesided

#####Diffregr (permutation-based p-values; 100 permutations)
##null case
fit3 <- diffregr_multisplit(y1,y3,x1,x2,b.splits=10,n.perm=100)
hist(fit3$pval.onesided,breaks=10)
fit3$aggpval.onesided
fit3$medpval.onesided

