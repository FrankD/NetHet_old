#########################################
##This an example of how to use DiffNet##
#########################################

rm(list=ls())
set.seed(1)

##genereta data
library(mvtnorm)
k <- 10
n <- 50
Sig1 <- diag(0,k)
for(j in 1:k){
    for (jj in 1:k){
        Sig1[j,jj] <- 0.5^{abs(j-jj)}
    }
}
SigInv1 <- solve(Sig1)
SigInv1[abs(SigInv1)<10^{-6}] <- 0
Sig2 <- SigInv2 <- diag(1,k)
x1 <- rmvnorm(n,mean = rep(0,k), sigma = Sig1)
x2 <- rmvnorm(n,mean = rep(0,k), sigma = Sig1)
x3 <- rmvnorm(n,mean = rep(0,k), sigma = Sig2)

##run diffnet
#Comparison x1,x2
fit.dn1 <- diffnet_multisplit(x1,x2,b.splits=10,
                              screen.meth='screen_bic.glasso',save.mle=TRUE)
par(mfrow=c(2,2))
image(x=1:k,y=1:k,abs(fit.dn1$medwi$modIpop1),xlab='',ylab='',main='median siginv1')
image(x=1:k,y=1:k,abs(fit.dn1$medwi$modIpop2),xlab='',ylab='',main='median siginv2')
hist(fit.dn1$pval.onesided,breaks=10);abline(v=fit.dn1$medpval.onesided,lty=2,col='red')

#Comparison x1,x3
fit.dn2 <- diffnet_multisplit(x1,x3,b.splits=10,
                              screen.meth='screen_bic.glasso',save.mle=TRUE)
par(mfrow=c(2,2))
image(x=1:k,y=1:k,abs(fit.dn2$medwi$modIpop1),xlab='',ylab='',main='median siginv1')
image(x=1:k,y=1:k,abs(fit.dn2$medwi$modIpop2),xlab='',ylab='',main='median siginv3')
hist(fit.dn2$pval.onesided,breaks=10);abline(v=fit.dn2$medpval.onesided,lty=2,col='red')



