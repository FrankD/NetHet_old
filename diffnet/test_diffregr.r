rm(list=ls())
#################
##Test diffregr##
#################

##Model
p <- 60
n <- 60
set.seed(1)
x1 <- matrix(rnorm(n*p),n,p)
x2 <- matrix(rnorm(n*p),n,p)
beta <- c(rep(1,5),rep(0,p-5))
y1 <- x1%*%as.matrix(beta)+rnorm(n,sd=0.5)
y2 <- x2%*%as.matrix(beta)+rnorm(n,sd=0.5)

source('twosample_highdimregr-09082012.r')
set.seed(1)
fit1 <- twosample_regr(y1,y2,x1,x2,b.splits=5)
str(fit1)
fit1$sspval.onesided
## > str(fit1)
## List of 9
##  $ pval.onesided   : num [1:5] 0.442 0.0244 0.0962 0.9142 0.13
##  $ pval.twosided   : num [1:5] 0.884 0.0488 0.1924 0.1715 0.2601
##  $ sspval.onesided : num 0.442
##  $ sspval.twosided : num 0.884
##  $ aggpval.onesided: num 1
##  $ aggpval.twosided: num 1
##  $ LR.last         : num 10.3
##  $ active.last     :List of 3
##   ..$ modJ    : int [1:36] 1 2 3 4 5 11 12 21 22 23 ...
##   ..$ modIpop1: int [1:13] 1 2 3 4 5 11 12 20 23 25 ...
##   ..$ modIpop2: int [1:22] 1 2 3 4 5 12 16 21 22 28 ...
##  $ beta.last       :List of 3
##   ..$ modJ    : num [1:36] 0.967 1.104 1.006 1.067 0.919 ...
##   ..$ modIpop1: num [1:13] 1.239 1.082 1.018 1.154 0.919 ...
##   ..$ modIpop2: num [1:22] 0.965 0.94 1.029 0.804 1.077 ...
## > fit1$sspval.onesided
## [1] 0.4419795

source('diffregr.r')
set.seed(1)
system.time(fit2 <- diffregr_multisplit(y1,y2,x1,x2,b.splits=5))
str(fit2)
fit2$sspval.onesided
## > str(fit2)
## List of 12
##  $ pval.onesided    : num [1:5] 0.442 0.0244 0.0962 0.9142 0.13
##  $ pval.twosided    : num [1:5] 0.884 0.0488 0.1924 0.1715 0.2601
##  $ sspval.onesided  : num 0.442
##  $ sspval.twosided  : num 0.884
##  $ medpval.onesided : num 0.13
##  $ medpval.twosided : num 0.192
##  $ aggpval.onesided : num 1
##  $ aggpval.twosided : num 1
##  $ teststat         : num [1:5] 5.14 47.7 30.27 -5.99 10.26
##  $ weights.nulldistr:List of 5
##   ..$ : num [1:48] 1 1 1 1 0.586 ...
##   ..$ : num [1:34] 1 1 1 1 1 1 1 1 1 1 ...
##   ..$ : num [1:45] 1 1 1 1 1 1 1 1 1 1 ...
##   ..$ : num [1:56] 1 1 1 1 1 ...
##   ..$ : num [1:60] 0.685 0.688 0.694 0.698 0.705 ...
##  $ active.last      :List of 3
##   ..$ modJ    : int [1:36] 1 2 3 4 5 11 12 21 22 23 ...
##   ..$ modIpop1: int [1:13] 1 2 3 4 5 11 12 20 23 25 ...
##   ..$ modIpop2: int [1:22] 1 2 3 4 5 12 16 21 22 28 ...
##  $ beta.last        :List of 3
##   ..$ modJ    : num [1:36] 0.967 1.104 1.006 1.067 0.919 ...
##   ..$ modIpop1: num [1:13] 1.239 1.082 1.018 1.154 0.919 ...
##   ..$ modIpop2: num [1:22] 0.965 0.94 1.029 0.804 1.077 ...
## > fit2$sspval.onesided
## [1] 0.4419795


set.seed(1)
system.time(fit3 <- diffregr_multisplit(y1,y2,x1,x2,b.splits=5,screen.meth='lasso.cvtrunc'))
str(fit3)
fit3$sspval.onesided
## > List of 12
##  $ pval.onesided    : num [1:5] 0.3301 0.9005 0.0825 0.2676 0.1976
##  $ pval.twosided    : num [1:5] 0.66 0.199 0.165 0.535 0.395
##  $ sspval.onesided  : num 0.33
##  $ sspval.twosided  : num 0.66
##  $ medpval.onesided : num 0.268
##  $ medpval.twosided : num 0.395
##  $ aggpval.onesided : num 1
##  $ aggpval.twosided : num 1
##  $ teststat         : num [1:5] 3.107 -0.411 8.189 3.946 5.253
##  $ weights.nulldistr:List of 5
##   ..$ : num [1:15] 1 0.926 0.95 1 1 ...
##   ..$ : num [1:11] 1 1 1 1 1 ...
##   ..$ : num [1:15] 1 0.812 0.917 1 1 ...
##   ..$ : num [1:15] 1 0.666 0.986 1 1 ...
##   ..$ : num [1:15] 1 0.92 0.965 1 1 ...
##  $ active.last      :List of 3
##   ..$ modJ    : int [1:12] 1 2 3 4 5 28 35 37 42 43 ...
##   ..$ modIpop1: int [1:6] 1 2 3 4 5 23
##   ..$ modIpop2: int [1:6] 1 2 3 4 5 22
##  $ beta.last        :List of 3
##   ..$ modJ    : num [1:12] 0.961 1.071 1.043 0.984 0.901 ...
##   ..$ modIpop1: num [1:6] 1.09 1.095 1.041 1.1 0.916 ...
##   ..$ modIpop2: num [1:6] 0.98 1.047 1.117 0.856 0.976 ...
## > fit3$sspval.onesided
## [1] 0.3301193
