
###########################################
##This an example of how to use MixGLasso##
###########################################

rm(list=ls())

##generate data
set.seed(1)
n <- 1000
K <- 3
p <- 10
mix.prob <- rep(1/K,K)
SigInv <- sparse_conc(p,K,s=p,s.common=0,scale.parcor=TRUE)
Sig <- array(NA,dim=c(p,p,K))
for (k in 1:K){
  Sig[,,k] <- solve(SigInv[[k]])
}
Mu <- matrix(0,p,K)
nonzero.mean <- split(sample(1:p),rep(1:K,length=p))
for(k in 1:K){
  Mu[nonzero.mean[[k]],k] <- -2/sqrt(ceiling(p/K))
}
dat <- simMIX(n,K,mix.prob,Mu,Sig)

##run mixglasso
set.seed(1)
fit1 <-  mixglasso_path(dat$X,n.comp=1:6,save.allfits=FALSE,show.trace=TRUE)
fit1$bic#17096.54 16823.18 16529.86 16585.69 16634.75 16698.64
set.seed(1)
fit2 <-  mixglasso_path(dat$X,n.comp=6,save.allfits=FALSE,show.trace=TRUE)
fit2$bic#16698.64
set.seed(1)
fit3 <-  mixglasso_path(dat$X,n.comp=1:6,lambda=0,save.allfits=FALSE,show.trace=TRUE)
set.seed(1)
fit4 <-  mixglasso_path(dat$X,n.comp=1:6,lambda=Inf,save.allfits=FALSE,show.trace=TRUE)
set.seed(1)
fit5 <-  bwprun_mixglasso(dat$X,n.comp=1,n.comp.max=6,selection.crit='bic')

##compare bic
plot(fit1$bic,type='l',lty=1,ylim=range(c(fit1$bic,fit3$bic,fit4$bic)),
     ylab='bic',xlab='no components')
lines(fit3$bic,lty=2)
lines(fit4$bic,lty=3)

