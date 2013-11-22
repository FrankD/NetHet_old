#################
##Test diffregr##
#################

#####clear workspace
rm(list=ls())

#####sessionInfo
sessionInfo()
## > sessionInfo()
## R version 2.15.2 (2012-10-26)
## Platform: i386-apple-darwin9.8.0/i386 (32-bit)

## locale:
## [1] C

## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     

#####model&data
p <- 60
n <- 60
set.seed(1)
x1 <- matrix(rnorm(n*p),n,p)
x2 <- matrix(rnorm(n*p),n,p)
beta <- c(rep(1,5),rep(0,p-5))
y1 <- x1%*%as.matrix(beta)+rnorm(n,sd=0.5)
y2 <- x2%*%as.matrix(beta)+rnorm(n,sd=0.5)

#####test diffregr
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
system.time(fit2 <- diffregr_multisplit(y1,y2,x1,x2,b.splits=5,screen.meth='lasso.cvmin',compute.evals='est2.my.ev2',method.compquadform='davies'))
str(fit2)
fit2$sspval.onesided
## > str(fit2)
## List of 12
##  $ pval.onesided    : num [1:5] 0.445081 1 0.000673 0.28719 0.073374
##  $ pval.twosided    : num [1:5] 0.89016 0 0.00135 0.57438 0.14675
##  $ sspval.onesided  : num 0.445
##  $ sspval.twosided  : num 0.89
##  $ medpval.onesided : num 0.287
##  $ medpval.twosided : num 0.147
##  $ aggpval.onesided : num 1
##  $ aggpval.twosided : num 0.0215
##  $ teststat         : num [1:5] 13.8 -52.4 57.4 14.5 37.9
##  $ weights.nulldistr:List of 5
##   ..$ : num [1:39] 1 1 1 1 1 1 1 1 1 1 ...
##   ..$ : num [1:54] 0.638 0.641 0.659 0.674 0.678 ...
##   ..$ : num [1:43] 1 1 1 1 1 1 1 1 1 1 ...
##   ..$ : num [1:31] 1 1 1 1 1 1 1 1 1 1 ...
##   ..$ : num [1:34] 1 1 1 1 1 1 1 1 1 1 ...
##  $ active.last      :List of 3
##   ..$ modJ    : int [1:11] 1 2 3 4 5 16 19 22 25 27 ...
##   ..$ modIpop1: int [1:14] 1 2 3 4 5 7 18 19 27 29 ...
##   ..$ modIpop2: int [1:22] 1 2 3 4 5 6 7 15 16 19 ...
##  $ beta.last        :List of 3
##   ..$ modJ    : num [1:11] 0.955 0.876 1.074 1.03 0.902 ...
##   ..$ modIpop1: num [1:14] 0.96 1.082 1.175 1.034 0.866 ...
##   ..$ modIpop2: num [1:22] 0.852 1.121 1.231 0.786 0.615 ...
## > fit2$sspval.onesided
## [1] 0.4450814

set.seed(1)
system.time(fit3 <- diffregr_multisplit(y1,y2,x1,x2,b.splits=5,screen.meth='lasso.cvtrunc',k.trunc=5,compute.evals='est2.my.ev2',method.compquadform='davies'))
str(fit3)
fit3$sspval.onesided
## > str(fit3)
## List of 12
##  $ pval.onesided    : num [1:5] 0.4998 0.97821 0.00661 0.53039 0.98191
##  $ pval.twosided    : num [1:5] 0.9996 0.0436 0.0132 0.9392 0.0362
##  $ sspval.onesided  : num 0.5
##  $ sspval.twosided  : num 1
##  $ medpval.onesided : num 0.53
##  $ medpval.twosided : num 0.0436
##  $ aggpval.onesided : num 1
##  $ aggpval.twosided : num 0.348
##  $ teststat         : num [1:5] 0.00214 -10.30788 13.41483 -0.32479 -10.80431
##  $ weights.nulldistr: num [1:12, 1:5] 1 1 1 1 1 1 -1 -1 -1 -1 ...
##  $ active.last      :List of 3
##   ..$ modJ    : int [1:11] 1 2 3 4 5 16 19 22 25 27 ...
##   ..$ modIpop1: int [1:5] 1 2 3 4 5
##   ..$ modIpop2: int [1:5] 1 2 3 4 5
##  $ beta.last        :List of 3
##   ..$ modJ    : num [1:11] 0.955 0.876 1.074 1.03 0.902 ...
##   ..$ modIpop1: num [1:5] 1.083 1.095 1.052 1.084 0.886
##   ..$ modIpop2: num [1:5] 0.977 0.774 1.054 0.964 0.993
## > fit3$sspval.onesided
## [1] 0.4997997

source('diffregr.r')
set.seed(1)
system.time(fit4 <- diffregr_multisplit(y1,y2,x1,x2,b.splits=5,screen.meth='lasso.cvmin',compute.evals='est2.my.ev3',method.compquadform='davies'))
fit2$sspval.onesided
fit4$sspval.onesided
fit2$pval.onesided
fit4$pval.onesided
## > fit2$sspval.onesided
## [1] 0.4450814
## > fit4$sspval.onesided
## [1] 0.4390524
## > fit2$pval.onesided
## [1] 0.4450813692 1.0000000000 0.0006729288 0.2871898281 0.0733738510
## > fit4$pval.onesided
## [1] 0.4390523578 1.0000000000 0.0005885578 0.2831460644 0.0719682973

source('diffregr.r')
set.seed(1)
system.time(fit5 <- diffregr_multisplit(y1,y2,x1,x2,b.splits=5,screen.meth='lasso.cvmin',compute.evals='est2.my.ev3',n.perm=250,method.compquadform='davies'))
fit5$pval.onesided# 0.404 0.472 0.756 0.320 0.608
fit5$teststat# 13.84779 83.50377 11.26196 80.92545 11.24971

#####test computations weights
source('diffregr.r')

#dimf>dimg
act1 <- 1:5
act2 <- 1:7
act12 <- c(1:5,6:8)

fit.ev1 <- est2.ww.mat2(y1,y2,x1,x2,beta[act1],beta[act2],beta[act12],act1,act2,act12)
fit.ev1
fit.ev2 <- est2.my.ev2(y1,y2,x1,x2,beta[act1],beta[act2],beta[act12],act1,act2,act12)
fit.ev2
fit.ev3 <- est2.my.ev3(y1,y2,x1,x2,beta[act1],beta[act2],beta[act12],act1,act2,act12)
fit.ev3
## > fit.ev1 <- est2.ww.mat2(y1,y2,x1,x2,beta[act1],beta[act2],beta[act12],act1,act2,act12)
## Warning messages:
## 1: In sqrt(1 - eval.mu) : NaNs produced
## 2: In sqrt(1 - eval.mu) : NaNs produced
## > fit.ev1
## $eval
##  [1]  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00           NaN
##  [6]           NaN           NaN  1.053671e-08  2.356080e-08  5.264434e-01
## [11]  6.445060e-01  9.896041e-01           NaN           NaN           NaN
## [16] -1.053671e-08 -2.356080e-08 -5.264434e-01 -6.445060e-01 -9.896041e-01

## $eval.mu.complex
## [1] 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 0.72285732 0.58461196
## [8] 0.02068375

## > fit.ev2 <- est2.my.ev2(y1,y2,x1,x2,beta[act1],beta[act2],beta[act12],act1,act2,act12)
## > fit.ev2
## $eval
##  [1]  1.0000000  1.0000000  1.0000000  1.0000000  1.0000000  0.0000000
##  [7]  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
## [13]  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.5375301
## [19]  0.6609152  1.0000000 -0.5375301 -0.6609152 -1.0000000

## $ev.aux.complex
## [1] 7.110614e-01 5.631911e-01 1.110223e-16

## > fit.ev3 <- est2.my.ev3(y1,y2,x1,x2,beta[act1],beta[act2],beta[act12],act1,act2,act12)
## > fit.ev3
## $eval
##  [1]  1.0000000  1.0000000  1.0000000  1.0000000  1.0000000  0.0000000
##  [7]  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
## [13]  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.5264434
## [19]  0.6445060  0.9896041 -0.5264434 -0.6445060 -0.9896041

## $ev.aux.complex
## [1] 0.72285732 0.58461196 0.02068375


#dimf<dimg
act1 <- 1:5
act2 <- 1:7
act12 <- c(1:5,6:40)

fit.ev1 <- est2.ww.mat2(y1,y2,x1,x2,beta[act1],beta[act2],beta[act12],act1,act2,act12)
fit.ev1
fit.ev2 <- est2.my.ev2(y1,y2,x1,x2,beta[act1],beta[act2],beta[act12],act1,act2,act12)
fit.ev2
fit.ev3 <- est2.my.ev3(y1,y2,x1,x2,beta[act1],beta[act2],beta[act12],act1,act2,act12)
fit.ev3
## > fit.ev1
## $eval
##  [1] -1.000000e+00 -1.000000e+00 -1.000000e+00 -1.000000e+00 -1.000000e+00
##  [6] -1.000000e+00 -1.000000e+00 -1.000000e+00 -1.000000e+00 -1.000000e+00
## [11] -1.000000e+00 -1.000000e+00 -1.000000e+00 -1.000000e+00 -1.000000e+00
## [16] -1.000000e+00 -1.000000e+00 -1.000000e+00 -1.000000e+00 -1.000000e+00
## [21] -1.000000e+00 -1.000000e+00 -1.000000e+00 -1.000000e+00 -1.000000e+00
## [26] -1.000000e+00 -1.000000e+00 -1.000000e+00           NaN  0.000000e+00
## [31]  0.000000e+00  2.580957e-08  3.494632e-08  4.708960e-01  5.770797e-01
## [36]  7.359270e-01  8.148511e-01  8.389942e-01  9.006553e-01  9.495093e-01
## [41]           NaN  0.000000e+00  0.000000e+00 -2.580957e-08 -3.494632e-08
## [46] -4.708960e-01 -5.770797e-01 -7.359270e-01 -8.148511e-01 -8.389942e-01
## [51] -9.006553e-01 -9.495093e-01

## $eval.mu.complex
##  [1] 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 0.77825693
##  [7] 0.66697901 0.45841141 0.33601763 0.29608874 0.18882002 0.09843215

## > fit.ev2 <- est2.my.ev2(y1,y2,x1,x2,beta[act1],beta[act2],beta[act12],act1,act2,act12)
## > fit.ev2
## $eval
##  [1] -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000
##  [7] -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000
## [13] -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000
## [19] -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000
## [25] -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000
## [31] -1.0000000 -1.0000000 -1.0000000  1.0000000  1.0000000  1.0000000
## [37]  1.0000000  1.0000000  1.0000000  0.0000000  0.0000000  0.0000000
## [43]  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
## [49]  0.0000000  0.0000000  0.0000000  0.4839946  0.5912149 -0.4839946
## [55] -0.5912149

## $ev.aux.complex
## [1] 0.7657493 0.6504649

## > fit.ev3 <- est2.my.ev3(y1,y2,x1,x2,beta[act1],beta[act2],beta[act12],act1,act2,act12)
## Warning message:
## In est2.my.ev3(y1, y2, x1, x2, beta[act1], beta[act2], beta[act12],  :
##   imaginary parts discarded in coercion
## > fit.ev3
## $eval
##  [1]  1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000
##  [7] -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000
## [13] -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000
## [19] -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000
## [25] -1.0000000 -1.0000000 -1.0000000 -1.0000000 -1.0000000  0.0000000
## [31]  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
## [37]  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.4708960
## [43]  0.5770797  0.7359270  0.8148511  0.8389942  0.9006553  0.9495093
## [49] -0.4708960 -0.5770797 -0.7359270 -0.8148511 -0.8389942 -0.9006553
## [55] -0.9495093

## $ev.aux.complex
##  [1]  7.782569e-01+0.000000e+00i  6.669790e-01+0.000000e+00i
##  [3]  4.584114e-01+0.000000e+00i  3.360176e-01+0.000000e+00i
##  [5]  2.960887e-01+0.000000e+00i  1.888200e-01+0.000000e+00i
##  [7]  9.843215e-02+0.000000e+00i -7.393706e-16+0.000000e+00i
##  [9]  6.394172e-16+0.000000e+00i -6.017896e-16+7.853033e-17i
## [11] -6.017896e-16-7.853033e-17i  5.557728e-16+0.000000e+00i
## [13]  4.532550e-16+0.000000e+00i -4.007209e-16+2.006420e-17i
## [15] -4.007209e-16-2.006420e-17i  3.782043e-16+0.000000e+00i
## [17]  2.854656e-16+1.973798e-16i  2.854656e-16-1.973798e-16i
## [19] -2.982079e-16+1.609121e-16i -2.982079e-16-1.609121e-16i
## [21]  2.184219e-16+1.765760e-16i  2.184219e-16-1.765760e-16i
## [23] -1.938774e-16+1.933223e-16i -1.938774e-16-1.933223e-16i
## [25] -2.589748e-16+0.000000e+00i  2.310617e-16+0.000000e+00i
## [27] -1.404250e-16+1.400645e-16i -1.404250e-16-1.400645e-16i
## [29]  1.681037e-16+0.000000e+00i -1.125662e-16+3.689483e-17i
## [31] -1.125662e-16-3.689483e-17i  7.445999e-17+0.000000e+00i
## [33] -7.798125e-18+5.713777e-17i -7.798125e-18-5.713777e-17i
## [35]  1.935940e-17+0.000000e+00i
