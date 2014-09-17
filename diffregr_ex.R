###############################################################
##This example illustrates the use of Differential Regression##
###############################################################


##set seed
set.seed(1)

##generate data
p <- 100
n <- 60
x1 <- matrix(rnorm(n*p),n,p)
x2 <- matrix(rnorm(n*p),n,p)
x3 <- matrix(rnorm(n*p),n,p)

act1 <- sample
beta1 <- c(rep(1,2),rep(0,p-2))
beta2 <- rep(0,p)
y1 <- x1%*%as.matrix(beta1)+rnorm(n,sd=0.5)
y2 <- x2%*%as.matrix(beta2)+rnorm(n,sd=0.5)
y3 <- x3%*%as.matrix(beta1)+rnorm(n,sd=0.5)

#####Diffregr (asymptotic p-values)
##alternative case
fit1 <- diffregr_multisplit(y1,y2,x1,x2,b.splits=10)
hist(fit1$ms.pval,breaks=10)
fit1$meinshagg.pval
fit1$medagg.pval
##null case
fit2 <- diffregr_multisplit(y1,y3,x1,x2,b.splits=10)
hist(fit2$ms.pval,breaks=10)
fit2$meinshagg.pval
fit2$medagg.pval

#####Diffregr (permutation-based p-values; 100 permutations)
##null case
fit3 <- diffregr_multisplit(y1,y3,x1,x2,b.splits=10,n.perm=100)
hist(fit3$ms.pval,breaks=10)
fit3$meinshagg.pval
fit3$medagg.pval



##Exit R
q('no')
