###############################################################
##This example illustrates the use of Differential Regression##
###############################################################

##set seed
set.seed(1)

##number of predictors / sample size
p <- 200
n <- 80

##predictor matrices
x1 <- matrix(rnorm(n*p),n,p)
x2 <- matrix(rnorm(n*p),n,p)

##active-sets and regression coefficients
act1 <- sample(1:p,5)
act2 <- c(act1[1:3],sample(setdiff(1:p,act1),2))
beta1 <- beta2 <- rep(0,p)
beta1[act1] <- 0.5
beta2[act2] <- 0.5

##response vectors under null-hypothesis
y1 <- x1%*%as.matrix(beta1)+rnorm(n,sd=0.5)
y2 <- x2%*%as.matrix(beta1)+rnorm(n,sd=0.5)

##Diffregr (asymptotic p-values)
fit.null <- diffregr_multisplit(y1,y2,x1,x2,b.splits=10)
plot(fit.null)
summary(fit.null)

##response vectors under alternative-hypothesis
y1 <- x1%*%as.matrix(beta1)+rnorm(n,sd=0.5)
y2 <- x2%*%as.matrix(beta2)+rnorm(n,sd=0.5)

##Diffregr (asymptotic p-values)
fit.alt <- diffregr_multisplit(y1,y2,x1,x2,b.splits=10)
plot(fit.alt)
fit.alt$meinshagg.pval
fit.alt$medagg.pval

##Diffregr (permutation-based p-values; 100 permutations)
fit.alt.perm <- diffregr_multisplit(y1,y2,x1,x2,b.splits=10,n.perm=100)
plot(fit.alt.perm)
summary(fit.alt.perm)


