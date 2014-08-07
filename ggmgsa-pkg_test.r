#######################################################
##This is an example of how to use DiffNet and GGMGSA##
#######################################################

rm(list=ls())

##genereta data
set.seed(1)
library(mvtnorm)
k <- 20
n <- 75
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

##run DiffNet
set.seed(1)
#Comparison x1,x2
fit.dn1 <- diffnet_multisplit(scale(x1),scale(x2),b.splits=10,
                              screen.meth='screen_bic.glasso',save.mle=TRUE)
par(mfrow=c(2,2))
image(x=1:k,y=1:k,abs(fit.dn1$medwi$modIpop1),xlab='',ylab='',main='median siginv1')
image(x=1:k,y=1:k,abs(fit.dn1$medwi$modIpop2),xlab='',ylab='',main='median siginv2')
hist(fit.dn1$pval.onesided,breaks=10);abline(v=fit.dn1$medpval.onesided,lty=2,col='red')

#Comparison x1,x3
fit.dn2 <- diffnet_multisplit(scale(x1),scale(x3),b.splits=10,
                              screen.meth='screen_bic.glasso',save.mle=TRUE)
par(mfrow=c(2,2))
image(x=1:k,y=1:k,abs(fit.dn2$medwi$modIpop1),xlab='',ylab='',main='median siginv1')
image(x=1:k,y=1:k,abs(fit.dn2$medwi$modIpop2),xlab='',ylab='',main='median siginv3')
hist(fit.dn2$pval.onesided,breaks=10);abline(v=fit.dn2$medpval.onesided,lty=2,col='red')

##run GGM-GSA
gsets <- list(1:5,3:10,11:15,14:20)

set.seed(1)
fit1 <- gsea.diffnet.multisplit(scale(x1),scale(x3),no.splits=5,gsets,1:20)
fit1$pval

p.adjust(apply(fit1$pval,1,aggpval),method='fdr') ##aggregated and fdr-corrected p-values
p.adjust(apply(fit1$pval,1,median),method='fdr')  ##median-aggregated and fdr-corrected p-values

options(cores=1)
set.seed(1)
fit2 <- par.gsea.diffnet.multisplit(scale(x1),scale(x3),no.splits=5,gsets,1:20,
                                    method.p.adjust='fdr',order.adj.agg='adj-agg')
fit2$pval









