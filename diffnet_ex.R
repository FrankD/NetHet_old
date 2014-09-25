
############################################################
##This example illustrates the use of Differential Network##
############################################################


##set seed
set.seed(1)

##sample size and number of nodes
n <- 70
p <- 30

##specifiy sparse inverse covariance matrices
gen.net <- generate_2networks(p,graph='random',n.nz=rep(p,2),
                              n.nz.common=ceiling(p*0.8))
invcov1 <- gen.net[[1]]
invcov2 <- gen.net[[2]]
plot_2networks(invcov1,invcov2,label.pos=0,label.cex=0.7)

##get corresponding correlation matrices
cor1 <- cov2cor(solve(invcov1))
cor2 <- cov2cor(solve(invcov2))

##generate data under null-hypothesis
library('mvtnorm')
x1 <- rmvnorm(n,mean = rep(0,p), sigma = cor1)
x2 <- rmvnorm(n,mean = rep(0,p), sigma = cor1)

##run diffnet (under null)
dn.null <- diffnet_multisplit(x1,x2,b.splits=10,verbose=FALSE)

plot(dn.null)#histogram of multi-split p-values
dn.null$medagg.pval#pvalue aggregated using median
dn.null$meinshagg.pval#pvalue aggregated using approach of Meinshausen et al (2009)


##generate data under alternative
x1 <- rmvnorm(n,mean = rep(0,p), sigma = cor1)
x2 <- rmvnorm(n,mean = rep(0,p), sigma = cor2)

##run diffnet (under alternative)
dn.altn <- diffnet_multisplit(x1,x2,b.splits=10,verbose=FALSE)

plot(dn.altn)#histogram of multi-split p-values
dn.altn$medagg.pval#pvalue aggregated using median
dn.altn$meinshagg.pval#pvalue aggregated using approach of Meinshausen et al (2009)


