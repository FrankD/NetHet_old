##Test diffregr

rm(list=ls())

##Model
p <- 100
n <- 100
set.seed(1)
x1 <- matrix(rnorm(n*p),n,p)
x2 <- matrix(rnorm(n*p),n,p)
beta <- c(rep(1,5),rep(0,p-5))
y1 <- x1%*%as.matrix(beta)+rnorm(n,sd=0.5)
y2 <- x2%*%as.matrix(beta)+rnorm(n,sd=0.5)

source('twosample_highdimregr-09082012.r')
set.seed(1)
fit1 <- twosample_regr(y1,y2,x1,x2)
str(fit1)

source('diffregr.r')
set.seed(1)
fit2 <- twosample_regr(y1,y2,x1,x2)
hist(fit1$pval.onesided)
str(fit2)
