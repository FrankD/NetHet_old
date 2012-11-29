##Test diffnet.r (computation weights of w-chi2; diffnet_multisplit)
##Gaussian Graphical Model
##

rm(list=ls())
source('twosample_diffnet-20072012.R')
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
##  [1] -0.853664 -0.718908 -0.707107 -0.707107 -0.549387 -0.281072  0.000000
##  [8]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
## [15]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
## [22]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
## [29]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
## [36]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
## [43]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
## [50]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
## [57]  0.000000  0.000000  0.281072  0.549387  0.707107  0.707107  0.718908
## [64]  0.853664  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [71]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [78]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [85]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [92]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000

##Compute weights (new function)
source('diffnet.r')
p <- k+(k+1)*k/2
act.mean <- (p-k+1):p
imat <- inf.mat(Sig,include.mean=TRUE)
e2 <- ww.mat(imat=imat,act=c(act,act.mean),act1=c(act1,act.mean),act2=c(act2,act.mean));round(sort(e2$eval),6)
## Warning message:
## In ww.mat(imat = imat, act = c(act, act.mean), act1 = c(act1, act.mean),  :
##   imaginary parts discarded in coercion
##   [1] -0.853664 -0.718908 -0.707107 -0.707107 -0.549387 -0.281072  0.000000
##   [8]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [15]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [22]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [29]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [36]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [43]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [50]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [57]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [64]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [71]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [78]  0.000000  0.281072  0.549387  0.707107  0.707107  0.718908  0.853664
##  [85]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
##  [92]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
##  [99]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [106]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [113]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [120]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [127]  1.000000

e3 <- est2.my.ev2(Sig,Sig,Sig,act1,act2,act,include.mean=TRUE);round(sort(e3$eval),6)
##   [1] -0.853664 -0.718908 -0.707107 -0.707107 -0.549387 -0.281072  0.000000
##   [8]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [15]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [22]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [29]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [36]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [43]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [50]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [57]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [64]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [71]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
##  [78]  0.000000  0.281072  0.549387  0.707107  0.707107  0.718908  0.853664
##  [85]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
##  [92]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
##  [99]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [106]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [113]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [120]  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000
## [127]  1.000000

e4 <- est2.ww.mat2(Sig,Sig,Sig,act1,act2,act,include.mean=TRUE);round(sort(e4$eval),6)

##Test Pval-Aggregation (old function)
set.seed(1)
source('twosample_diffnet-20072012.R')
fit.pval1 <- twosample_diffnet2(xx1,xx2,b.splits=5,frac.split=1/2,screen.lambda='lambda.cv',gamma.min=0.05,compute.evals='est2.my.ev2',diag.invcov=TRUE,lambda=la,folds=10)
fit.pval1$pval.onesided# 0.8262682 0.8603141 0.9607790 0.2258820 0.2905430
fit.pval1$aggpval.onesided# 1
fit.pval1$LR.last#38.96125

##Test Pval-Aggregation (new function)
set.seed(1)
source('diffnet.r')
fit.pval2 <- diffnet_multisplit(xx1,xx2,b.splits=5,include.mean=FALSE,lambda=la)
fit.pval2$pval.onesided# 0.8262682 0.8603141 0.9607790 0.2258820 0.2905430
fit.pval2$aggpval.onesided# 1
fit.pval2$teststat#26.07755 26.07807 21.08478 50.16631 38.96125

set.seed(1)
source('diffnet.r')
fit.pval3 <- diffnet_multisplit(xx1,xx2,b.splits=5,include.mean=TRUE,lambda=la)
fit.pval3$pval.onesided# 0.9999916 0.9996231 0.9999693 0.9960757 0.9845568
fit.pval3$teststat#2.030517 15.938588  7.216169 28.174944 22.661378


## fit.pval$LR.last#377.3058
## set.seed(1)
## n.sim <- 100
## eval.wwmat <- eval.wwmat2 <- eval.myev2 <- eval.2myev2 <- matrix(NA,length(e3$eval),n.sim)
## n <- 1000
## for (i in 1:n.sim){
##     cat(i,'\n')
##     x.v1 <- rmvnorm(n,mean = rep(0,p), sigma = Sig)
##     x.v2 <- rmvnorm(n,mean = rep(0,p), sigma = Sig)
##     x.v <- rbind(x.v1,x.v2)

##     sig1 <- glasso(crossprod(x.v1)/nrow(x.v1),rho=0,zero=which(fit.cv1$wi==0,arr.in=TRUE))$w
##     sig2 <- glasso(crossprod(x.v2)/nrow(x.v2),rho=0,zero=which(fit.cv2$wi==0,arr.in=TRUE))$w
##     sig <- glasso(crossprod(x.v)/nrow(x.v),rho=0,zero=which(fit.cv$wi==0,arr.in=TRUE))$w

##     source('../code/twosample_diffnet-25062012.r')
##     eval.wwmat[,i] <- sort(est.ww.mat(x.v1,x.v2,sig1,sig2,sig,act1,act2,act)$eval,na.last=TRUE)
##     eval.wwmat2[,i] <- sort(est.ww.mat2(x.v1,x.v2,sig1,sig2,sig,act1,act2,act)$eval,na.last=TRUE)
##     eval.myev2[,i] <- sort(est.my.ev2(x.v1,x.v2,sig1,sig2,sig,act1,act2,act)$eval,na.last=TRUE)
##     source('../code/twosample_diffnet-20072012.r')
##     eval.2myev2[,i] <- sort(est2.my.ev2(sig1,sig2,sig,act1,act2,act)$eval,na.last=TRUE)
## }

## par(mfrow=c(2,2),mar=c(2,3,2,2),mgp=c(2,0.8,0))
## boxplot(t(eval.wwmat),ylab=paste('n=',n));points(1:length(e1$eval),sort(e1$eval),col='red',pch=4)
## #boxplot(t(eval.wwmat2));points(1:length(e1$eval),sort(e1$eval),col='red',pch=4)
## #boxplot(t(eval.myev2));points(1:length(e1$eval),sort(e1$eval),col='red',pch=4)
## boxplot(t(eval.2myev2));points(1:length(e1$eval),sort(e1$eval),col='red',pch=4)


## set.seed(1)
## n <- 500
## lograt <- rep(NA,n)
## for (i in 1:length(lograt)){
##     cat(i,'\n')
##     v.xx1 <- rmvnorm(n,mean = rep(0,p), sigma = Sig)
##     v.xx2 <- rmvnorm(n,mean = rep(0,p), sigma = Sig)
##     v.xx <- rbind(v.xx1,v.xx2)
##     lograt[i] <- logratio(v.xx1,v.xx2,v.xx,Sig,Sig,fit.cv$wi)$twiceLR
## }
## rchi <- rep(NA,n)
## for (i in 1:n){
## rchi[i] <- sum(hate2*rchisq(length(hate2),df=1))
## }
## hist(lograt,breaks=50,prob=TRUE)
## lines(density(rchi),col='blue')
