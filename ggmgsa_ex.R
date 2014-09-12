#######################################################
##This example illustrates the use of GGMGSA         ##
#######################################################

rm(list=ls())

##genereta data
set.seed(1)
p <- 50#network with p nodes
n <- 50
hub.net <- generate.2networks(p,graph='hub',n.hub=5,n.hub.diff=2)#generate hub networks
invcov1 <- hub.net[[1]]
invcov2 <- hub.net[[2]]
plot.2networks(invcov1,invcov2,label.pos=0,label.cex=0.7)

##generate data
x1 <- rmvnorm(n,mean = rep(0,p), sigma = cov2cor(solve(invcov1)))
x2 <- rmvnorm(n,mean = rep(0,p), sigma = cov2cor(solve(invcov2)))

##run DiffNet
fit.dn <- diffnet_multisplit(x1,x2,b.splits=5)
plot(fit.dn)

##identify hubs with 'gene-sets'
gene.names <- paste('G',1:p,sep='')
gsets <- split(gene.names,rep(1:5,each=10))

##run GGMGSA
options(cores=5)
fit.ggmgsa <- ggmgsa_multisplit_par(x1,x2,no.splits=5,gsets,gene.names)
plot(fit.ggmgsa)
summary(fit.ggmgsa)
fit.ggmgsa$medagg.pval
p.adjust(apply(fit.ggmgsa$pval,1,median),method='fdr')






