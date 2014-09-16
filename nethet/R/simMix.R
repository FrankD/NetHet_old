
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


##' Generate an inverse covariance matrix with a given sparsity and dimensionality
##'
##' This function generates an inverse covariance matrix, with at most (1-sparsity)*p(p-1)
##' non-zero off-diagonal entries, where the non-zero entries are sampled from a 
##' beta distribution.
##' 
##' @title generateInvCov
##' @param p Dimensionality of the matrix.
##' @param sparsity Determined the proportion of non-zero off-diagonal entries.
##' @return A p by p positive definite inverse covariance matrix.
##' @export
generateInvCov <- function(p=162, sparsity=0.7) {
	num.edges = p*(p-1)/2
	
	edge.prop = 1 - sparsity
	s = round(num.edges*edge.prop) # Number of non-zero entries in each matrix
	
	return(getinvcov(p, s=s))
}

##' Generate an inverse covariance matrix with a given sparsity and dimensionality
##'
##' Internal function
getinvcov<- function(p,s, a.diff=5,b.diff=5,magn.diag=0,emin=0.1){
	#####8!!! act1 are the indices of the upper-diagonal non-zero entries of a pxp matrix 
	
	ind.upper.tri <- which(upper.tri(matrix(NA,p,p)))
	act1 <- sample(ind.upper.tri,size=s,replace=FALSE)
  
	B1 <- matrix(0,p,p)
	B1[act1] <- rbeta(length(act1),a.diff,b.diff)
	diag(B1) <- 0
	B1 <- (t(B1)+B1)
	
	
	####compute (positive definite) concentration matrices
	SigInv <- list()
	
	siginv <- B1
	e.min <- min(eigen(siginv)$values)
	siginv <- siginv+diag(abs(e.min)+emin+magn.diag,p)
	SigInv <- siginv
	
	return(SigInv)
}



##' Generate inverse covariances, means, mixing probabilities, and simulate 
##' data from resulting mixture model.
##'
##' This function generates n.comp mean vectors from a standard Gaussian and
##' n.comp covariance matrices, with at most (1-sparsity)*p(p-1)/2
##' non-zero off-diagonal entries, where the non-zero entries are sampled from a 
##' beta distribution. Then it uses {\link{simMIX}} to simulate from a 
##' mixture model with these means and covariance matrices.
##' 
##' Means Mu and covariance matrices Sig can also be supplied by the user.
##' 
##' @title sim_mix_networks
##' @param n Number of data points to simulate.
##' @param p Dimensionality of the data.
##' @param n.comp Number of components of the mixture model.
##' @param sparsity Determines the proportion of non-zero off-diagonal entries.
##' @param mix.prob Mixture probabilities for the components; defaults to uniform distribution.
##' @param Mu Means for the mixture components, a p by n.comp matrix. If NULL, 
##' sampled from a standard Gaussian.
##' @param Sig Covariances for the mixture components, a p by p by n.comp array. If NULL,
##' generated using {\link{generateInvCov}}.
##' @return A list with components:
##' \item{Mu} Means of the mixture components.
##' \item{Sig} Covariances of the mixture components.
##' \item{data} Simulated data, a n by p matrix.
##' \item{S} Component assignments, a vector of length n.
##' @export
sim_mix_networks <- function(n, p, n.comp, sparsity=0.7, 
														 mix.prob=rep(1/n.comp, n.comp),
														 Mu=NULL, Sig=NULL, ...) {
	
	if(is.null(Mu)) {
	  Mu = sapply(1:n.comp, function(n.comp) rnorm(p,0,3))
	}
	
	if(is.null(Sig)) {
	  Sig = sapply(1:n.comp, 
												function(n.comp) solve(generateInvCov(p, sparsity)), 
												simplify='array')
	}
	
	data = simMIX(n, n.comp, mix.prob, Mu, Sig, ...)
	
	return(list(Mu=Mu, Sig=Sig, data=data$X, comp=data$S))
}
