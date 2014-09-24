#######################################################
##This example illustrates the use of GGMGSA         ##
#######################################################


##genereta data
set.seed(1)
p <- 50#network with p nodes
n <- 50
hub.net <- generate.2networks(p,graph='hub',n.hub=5,n.hub.diff=2)#generate hub networks
invcov1 <- hub.net[[1]]
invcov2 <- hub.net[[2]]
plot_2networks(invcov1,invcov2,label.pos=0,label.cex=0.7)

##generate data
library('mvtnorm')
x1 <- rmvnorm(n,mean = rep(0,p), sigma = cov2cor(solve(invcov1)))
x2 <- rmvnorm(n,mean = rep(0,p), sigma = cov2cor(solve(invcov2)))

##run DiffNet
fit.dn <- diffnet_multisplit(x1,x2,b.splits=10,verbose=FALSE)
fit.dn$medagg.pval

##identify hubs with 'gene-sets'
gene.names <- paste('G',1:p,sep='')
gsets <- split(gene.names,rep(1:5,each=10))

##run GGM-GSA
fit.ggmgsa <- ggmgsa_multisplit(x1,x2,no.splits=10,gsets,gene.names,verbose=FALSE)
plot(fit.ggmgsa)
fit.ggmgsa$medagg.pval
p.adjust(apply(fit.ggmgsa$pval,1,median),method='fdr')






