################################################################
################################################################
##This File is for testing mixglasso

rm(list=ls())

##source functions
source('hmmgl_mixgl-22102012-backup-mixglasso-only.r')

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
fit1 <-  mixglasso_par(dat$X,nr.states=1:6,save.allfits=FALSE,show.trace=TRUE)
set.seed(1)
fit2 <-  mixglasso_par(dat$X,nr.states=1:6,lambda=0,save.allfits=FALSE,show.trace=TRUE)
set.seed(1)
fit3 <-  mixglasso_par(dat$X,nr.states=1:6,lambda=Inf,save.allfits=FALSE,show.trace=TRUE)
set.seed(1)
fit4 <-  bwprun_mixglasso(dat$X,nr.states.min=1,nr.states.max=6,selection.crit='bic')

##compare bic
plot(fit1$bic,type='l',lty=1,ylim=range(c(fit1$bic,fit2$bic,fit3$bic,fit4$selcrit)),ylab='bic',xlab='no components')
lines(fit2$bic,lty=2)
lines(fit3$bic,lty=3)
lines(fit4$selcrit,lty=4)
