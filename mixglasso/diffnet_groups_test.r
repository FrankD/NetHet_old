##############################################
##This an example of how to use DiffNetGroup##
##############################################

#################################
##clear work space and set.seed##
#################################
rm(list=ls())
set.seed(1)

##########
## data ##
##########
library(mvtnorm)
k <- 10##number of proteins
n <- 50##number of samples per disease

##network1: Sig1 (autoregressive process AR(1))
Sig1 <- diag(0,k)
for(j in 1:k){
    for (jj in 1:k){
        Sig1[j,jj] <- 0.5^{abs(j-jj)}
    }
}
SigInv1 <- solve(Sig1)
SigInv1[abs(SigInv1)<10^{-6}] <- 0

##network2: identity matrix (fully unconnected network)
Sig2 <- SigInv2 <- diag(1,k)

##generate data
x1 <- rmvnorm(n,mean = rep(0,k), sigma = Sig1)##data disease 1
x2 <- rmvnorm(n,mean = rep(0,k), sigma = Sig1)##data disease 2
x3 <- rmvnorm(n,mean = rep(0,k), sigma = Sig2)##data disease 3
x <- rbind(x1,x2,x3)##data-matrix with all diseases combined
groups <- factor(c(rep(1,50),rep(2,50),rep(3,50)))##vector specifying 3 diseases

########################################
##run diffnet comparing diseases 1,2,3##
########################################
source('diffnet_groups.R')
options(cores=1)
set.seed(1)
fit.dn <- par.diffnet_groups_multisplit(x,groups,no.splits=10,
                                        method.p.adjust='fdr',order.adj.agg='agg-adj')
convert2mat(fit.dn$pvalmed,levels(groups))
## > convert2mat(fit.dn$pvalmed,levels(groups))
##             1          2           3
## 1 0.000000000 0.73025212 0.001276888
## 2 0.730252123 0.00000000 0.084076075
## 3 0.001276888 0.08407607 0.000000000

p.adjust(apply(fit.dn$pval,1,median),method='fdr')
## 0.730252123 0.001276888 0.084076075   

