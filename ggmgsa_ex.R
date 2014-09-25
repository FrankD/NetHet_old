#######################################################
##This example illustrates the use of GGMGSA         ##
#######################################################


##genereta data
set.seed(1)
p <- 9#network with p nodes
n <- 40
hub.net <- generate_2networks(p,graph='hub',n.hub=3,n.hub.diff=1)#generate hub networks
invcov1 <- hub.net[[1]]
invcov2 <- hub.net[[2]]
plot_2networks(invcov1,invcov2,label.pos=0,label.cex=0.7)

##generate data
library('mvtnorm')
x1 <- rmvnorm(n,mean = rep(0,p), sigma = cov2cor(solve(invcov1)))
x2 <- rmvnorm(n,mean = rep(0,p), sigma = cov2cor(solve(invcov2)))

##run DiffNet
# fit.dn <- diffnet_multisplit(x1,x2,b.splits=2,verbose=FALSE)
# fit.dn$medagg.pval

##identify hubs with 'gene-sets'
gene.names <- paste('G',1:p,sep='')
gsets <- split(gene.names,rep(1:3,each=3))

##run GGM-GSA
fit.ggmgsa <- ggmgsa_multisplit(x1,x2,no.splits=2,gsets,gene.names,verbose=FALSE)
summary(fit.ggmgsa)
fit.ggmgsa$medagg.pval#median aggregated p-values
p.adjust(apply(fit.ggmgsa$pval,1,median),method='fdr')#or: first median aggregation,
                                                      #second fdr-correction






