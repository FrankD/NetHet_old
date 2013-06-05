##Test diffnet.r (computation weights of w-chi2; diffnet_multisplit)
##Gaussian Graphical Model
##

rm(list=ls())
source('twosample_diffnet-20072012.R')
sessionInfo()

## R version 2.15.2 (2012-10-26)
## Platform: i386-apple-darwin9.8.0/i386 (32-bit)

## locale:
## [1] C

## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     

## other attached packages:
## [1] CompQuadForm_1.4 glasso_1.7       mvtnorm_0.9-9994

## loaded via a namespace (and not attached):
## [1] tools_2.15.2

#################
##Generate Data##
#################

##Model
k <- 10
n <- 200
Sig <- diag(0,k)
for(j in 1:k){
    for (jj in 1:k){
        Sig[j,jj] <- 0.9^{abs(j-jj)}
    }
}
SigInv <- solve(Sig)
SigInv[abs(SigInv)<10^{-6}] <- 0

##Generate Data
set.seed(1315)
xx1 <- rmvnorm(n,mean = rep(0,k), sigma = Sig)
xx2 <- rmvnorm(n,mean = rep(0,k), sigma = Sig)

##Select active-sets
la <- lambdagrid_mult(0.1,20,30)[30:1]
fit.cv1 <- cv.glasso.1se(x=xx1,lambda=la)
fit.cv2 <- cv.glasso.1se(x=xx2,lambda=la)
fit.cv <- cv.glasso.1se(x=rbind(xx1,xx2),lambda=la)

act1 <- which(fit.cv1$wi[upper.tri(Sig,diag=TRUE)]!=0)
act2 <- which(fit.cv2$wi[upper.tri(Sig,diag=TRUE)]!=0)
act <- which(fit.cv$wi[upper.tri(Sig,diag=TRUE)]!=0)

##Compute weights (old function)
source('twosample_diffnet-20072012.R')
p <- (k+1)*k/2
imat <- inf.mat(Sig,1:p)
e1 <- ww.mat(imat,act,act1,act2);round(sort(e1$eval),6)
## Warning message:
## In ww.mat(imat, act, act1, act2) : imaginary parts discarded in coercion
##  [1] -0.913500 -0.903837 -0.718028 -0.692951 -0.406839 -0.401952 -0.287048
##  [8]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
## [15]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
## [22]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
## [29]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
## [36]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
## [43]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
## [50]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.287048
## [57]  0.401952  0.406839  0.692951  0.718028  0.903837  0.913500  1.000000
## [64]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [71]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [78]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [85]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [92]  1.000000  1.000000  1.000000

##Compute weights (new function)
source('diffnet.r')
dyn.load("../code/betamat_diffnet.so")
p <- k+(k+1)*k/2
act.mean <- (p-k+1):p
imat <- inf.mat(Sig,include.mean=TRUE)
e2 <- ww.mat(imat=imat,act=c(act,act.mean),act1=c(act1,act.mean),act2=c(act2,act.mean));round(sort(e2$eval),6)
## Warning message:
## In ww.mat(imat = imat, act = c(act, act.mean), act1 = c(act1, act.mean),  :
##   imaginary parts discarded in coercion
##   [1] -0.913500 -0.903837 -0.718028 -0.692951 -0.406839 -0.401952 -0.287048
##   [8]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [15]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [22]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [29]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [36]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [43]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [50]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [57]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [64]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [71]  0.000000  0.000000  0.000000  0.000000  0.000000  0.287048  0.401952
##  [78]  0.406839  0.692951  0.718028  0.903837  0.913500  1.000000  1.000000
##  [85]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
##  [92]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
##  [99]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [106]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [113]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [120]  1.000000  1.000000  1.000000  1.000000  1.000000

e3 <- est2.my.ev2(Sig,Sig,Sig,act1,act2,act,include.mean=TRUE);round(sort(e3$eval),6)
##   [1] -0.913500 -0.903837 -0.718028 -0.692951 -0.406839 -0.401952 -0.287048
##   [8]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [15]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [22]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [29]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [36]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [43]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [50]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [57]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [64]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [71]  0.000000  0.000000  0.000000  0.000000  0.000000  0.287048  0.401952
##  [78]  0.406839  0.692951  0.718028  0.903837  0.913500  1.000000  1.000000
##  [85]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
##  [92]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
##  [99]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [106]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [113]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [120]  1.000000  1.000000  1.000000  1.000000  1.000000

e4 <- est2.ww.mat2(Sig,Sig,Sig,act1,act2,act,include.mean=TRUE);round(sort(e4$eval),6)

##Test Pval-Aggregation (old function)
set.seed(1)
source('twosample_diffnet-20072012.R')
fit.pval1 <- twosample_diffnet2(xx1,xx2,b.splits=5,frac.split=1/2,screen.lambda='lambda.cv',gamma.min=0.05,compute.evals='est2.my.ev2',diag.invcov=TRUE,lambda=la,folds=10)
fit.pval1$pval.onesided# 0.37426285 0.06497231 0.07142485 0.20118175 0.15120453
fit.pval1$aggpval.onesided# 1
fit.pval1$LR.last# 42.14646

##Test Pval-Aggregation (new function)
set.seed(1)
source('diffnet.r')
dyn.load("../code/betamat_diffnet.so")
fit.pval2 <- diffnet_multisplit(xx1,xx2,b.splits=5,include.mean=FALSE,screen.meth='cv.glasso',algorithm.mleggm='glasso_rho0',lambda=la)
fit.pval2$pval.onesided# 0.37426285 0.06497231 0.07142485 0.20118175 0.15120453
fit.pval2$aggpval.onesided# 1
fit.pval2$teststat#31.54678 53.48060 41.56653 42.62735 42.14646

set.seed(1)
source('diffnet.r')
dyn.load("../code/betamat_diffnet.so")
fit.pval3 <- diffnet_multisplit(xx1,xx2,b.splits=5,include.mean=TRUE,screen.meth='cv.glasso',algorithm.mleggm='glasso_rho0',lambda=la)
fit.pval3$pval.onesided# 0.41430313 0.03970084 0.04280262 0.24598699 0.08501593
fit.pval3$teststat# 38.72533 68.02944 56.17720 51.79453 57.12753


