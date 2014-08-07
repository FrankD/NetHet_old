
#####################
##Required Packages##
#####################
library(mvtnorm)

##' Simulate from mixture model with multi-variate Gaussian or t-distributed components.
##'
##' 
##' @title Simulate from mixture model.
##' @param n sample size
##' @param n.comp number of mixture components ("comps")
##' @param mix.prob mixing probablities (need to sum to 1)
##' @param Mu matrix of component-specific mean vectors 
##' @param Sig array of component-specific covariance matrices
##' @param dist 'norm' for Gaussian components, 't' for t-distributed components
##' @param df degrees of freedom of the t-distribution (not used for Gaussian distribution), default=2
##' @return  a list consisting of:
##' \item{S}{component assignments}
##' \item{X}{observed data matrix}
##' @author n.stadler
##' @export
simMIX <- function(n,n.comp,mix.prob,Mu,Sig, dist='norm', df=2){
	
  # Only one Mu specified
  if(is.vector(Mu) || dim(Mu)[2] == 1) 
    Mu = matrix(Mu, length(Mu), n.comp)
  
  # Only one Sig specified
  if(length(dim(Sig)) == 2 || dim(Sig)[3] == 1) 
    Sig = array(Sig, c(dim(Sig)[1:2], n.comp))
     
  # Incorrect dimension of Mu or Sig
  if(is.vector(Mu) || dim(Mu)[2] != n.comp ||
     length(dim(Sig)) != 3 || dim(Sig)[3] != n.comp) 
    stop(paste('The last dimension of Mu and Sig needs to be equal to n.comp',
               '(one mean vector and covariance matrix per component).'))
  
  # Incorrect number of mixture probabilities
  if(length(mix.prob) != n.comp) 
    stop('Length of mix.prob needs to be equal to n.comp.')
  
  # Misspecified mixture probabilities
  if(any(mix.prob > 1 | mix.prob < 0) || sum(mix.prob) != 1) 
    stop('Misspecified mixture probabilities.')
  
  K <- n.comp
	p <- dim(Mu)[1]
	x <- matrix(0,nrow=n,ncol=p)
	s <- sample(1:K,n,replace=TRUE,prob=mix.prob)
	
  # Generate from each component
  for (k in 1:K){    
		x[s==k,] <-  if(dist == 'norm') rmvnorm(sum(s==k),mean=Mu[,k],sigma=Sig[,,k])
                 else if(dist == 't') rmvt(sum(s==k),delta=Mu[,k],sigma=Sig[,,k], df=df, type='shifted')
                 else stop(paste('Invalid dist argument:', dist))
	}
	
  return(list(S=s,X=x))
}
