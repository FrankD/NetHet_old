

###########################################
##This an example of how to use MixGLasso##
###########################################

##generate data
set.seed(1)
n <- 1000
n.comp <- 3
p <- 10

# Create different mean vectors
Mu <- matrix(0,p,n.comp)

nonzero.mean <- split(sample(1:p),rep(1:n.comp,length=p))
for(k in 1:n.comp){
  Mu[nonzero.mean[[k]],k] <- -2/sqrt(ceiling(p/n.comp))
}

sim <- sim_mix_networks(n, p, n.comp, Mu=Mu)

##run mixglasso
set.seed(1)
fit1 <-  mixglasso(sim$data,n.comp=1:6)
fit1$bic
set.seed(1)
fit2 <-  mixglasso(sim$data,n.comp=6)
fit2$bic
set.seed(1)
fit3 <-  mixglasso(sim$data,n.comp=1:6,lambda=0)
set.seed(1)
fit4 <-  mixglasso(sim$data,n.comp=1:6,lambda=Inf)
#set.seed(1)
#fit5 <-  bwprun_mixglasso(dat$X,n.comp=1,n.comp.max=6,selection.crit='bic')

##compare bic
library('ggplot2')
plotting.frame <- data.frame(BIC= c(fit1$bic, fit3$bic, fit4$bic), 
														 Num.Comps=rep(1:6, 3), Lambda=rep( 
														 	                        c('Default', 
														 													  'Lambda = 0',
														 														'Lambda = Inf'),
														 													each=6))

p <- ggplot(plotting.frame) + 
	geom_line(aes(x=Num.Comps, y=BIC, colour=Lambda))

print(p)
