
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
##' @import mvtnorm
##' @examples
##' n.comp = 4
##' p = 5 # dimensionality
##' Mu = matrix(rep(0, p), p, n.comp)
##' Sigma = array(diag(p), c(p, p, n.comp))
##' mix.prob = rep(0.25, n.comp)
##' 
##' sim_mix(100, n.comp, mix.prob, Mu, Sigma)

sim_mix <- function(n,n.comp,mix.prob,Mu,Sig, dist='norm', df=2){
	
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
##' @title generate_inv_cov
##' @param p Dimensionality of the matrix.
##' @param sparsity Determined the proportion of non-zero off-diagonal entries.
##' @return A p by p positive definite inverse covariance matrix.
##' @export
##' @examples
##' generate_inv_cov(p=162)
generate_inv_cov <- function(p=162, sparsity=0.7) {
	num.edges = p*(p-1)/2
	
	edge.prop = 1 - sparsity
	s = round(num.edges*edge.prop) # Number of non-zero entries in each matrix
	
	return(getinvcov(p, s=s))
}

##' Generate an inverse covariance matrix with a given sparsity and dimensionality
##'
##' @param p Dimensionality
##' @param s Sparsity
##' @param a.diff binomial parameter
##' @param b.diff binomial parameter
##' @param magn.diag Magnitude
##' @param emin e min
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
##' beta distribution. Then it uses \code{\link{sim_mix}} to simulate from a 
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
##' generated using \code{\link{generate_inv_cov}}.
##' @param ... Further arguments passed to \code{\link{sim_mix}}.
##' @return A list with components:
##' \code{Mu} Means of the mixture components.
##' \code{Sig} Covariances of the mixture components.
##' \code{data} Simulated data, a n by p matrix.
##' \code{S} Component assignments, a vector of length n.
##' @export
##' @examples
##' # Generate dataset with 100 samples of dimensionality 30, and 4 components
##' test.data = sim_mix_networks(n=100, p=30, n.comp=4)
sim_mix_networks <- function(n, p, n.comp, sparsity=0.7, 
														 mix.prob=rep(1/n.comp, n.comp),
														 Mu=NULL, Sig=NULL, ...) {
	
	if(is.null(Mu)) {
	  Mu = sapply(1:n.comp, function(n.comp) rnorm(p,0,3))
	}
	
	if(is.null(Sig)) {
	  Sig = sapply(1:n.comp, 
												function(n.comp) solve(generate_inv_cov(p, sparsity)), 
												simplify='array')
	}
	
	data = sim_mix(n, n.comp, mix.prob, Mu, Sig, ...)
	
	return(list(Mu=Mu, Sig=Sig, data=data$X, comp=data$S))
}

##' Generate two sparse inverse covariance matrices with overlap
##'
##' 
##' @title Generate sparse invcov with overlap
##' @param p number of nodes
##' @param graph 'random' or 'hub'
##' @param n.hub number of hubs (only for graph='hub')
##' @param n.hub.diff number of different hubs
##' @param n.nz number of edges per graph (only for graph='random')
##' @param n.nz.common number of edges incommon between graphs (only for graph='random')
##' @param magn.nz.diff default=0.9
##' @param magn.nz.common default=0.9
##' @param magn.diag default=0
##' @param emin default=0.1 (see ?huge.generator)
##' @param verbose If verbose=FALSE then tracing output is disabled.
##' @export
##' @return Two sparse inverse covariance matrices with overlap
##' @examples
##' n <- 70
##' p <- 30
##' 
##' ## Specifiy sparse inverse covariance matrices,
##' ## with number of edges in common equal to ~ 0.8*p
##' gen.net <- generate_2networks(p,graph='random',n.nz=rep(p,2),
##'                               n.nz.common=ceiling(p*0.8))
##' 
##' invcov1 <- gen.net[[1]]
##' invcov2 <- gen.net[[2]]
##' 
##' plot_2networks(invcov1,invcov2,label.pos=0,label.cex=0.7)
generate_2networks<- function(p,graph='random',
															n.nz=rep(p,2),n.nz.common=p,
															n.hub=2,n.hub.diff=1,
															magn.nz.diff=0.8,
															magn.nz.common=0.9,magn.diag=0,emin=0.1,
															verbose=FALSE){
	
	if(graph=='random'){
		
		if(min(n.nz)<n.nz.common){stop('min(n.nz) < n.nz.common')}
		
		ind.upper.tri <- which(upper.tri(matrix(NA,p,p)))
		B.list <- list()
		B1 <- matrix(0,p,p)
		nz1 <- sample(ind.upper.tri,size=n.nz[1],replace=FALSE)
		nz12 <- sample(nz1,size=n.nz.common,replace=FALSE)
		remain.zero <- setdiff(ind.upper.tri,nz1)
		B1[nz12] <- magn.nz.common
		B1[setdiff(nz1,nz12)] <- magn.nz.diff
		B.list[[1]] <- B1+t(B1)
		
		B2 <- matrix(0,p,p)
		nz2 <- c(nz12,sample(remain.zero,size=n.nz[2]-n.nz.common,replace=FALSE))
		B2[nz2] <- magn.nz.common
		B2[setdiff(nz2,nz12)] <- magn.nz.diff
		B.list[[2]] <- B2+t(B2)
		
	}
	
	if(graph=='hub'){
		
		if(n.hub<n.hub.diff){stop("n.hub less than n.hub.diff: choose smaller n.hub.diff")}
		
		##generate hub network (see library(huge); huge.generator)
		theta.hub <- matrix(0,p,p)
		g.large = p%%n.hub#number of large hubs
		g.small = n.hub - g.large#number of small hubs
		n.small = floor(p/n.hub)#size small hub
		if(n.small<=1){
			stop('hub with less than 2 nodes: choose a smaller n.hub')
		}
		n.large = n.small + 1#size large hub
		g.list = c(rep(n.small, g.small), rep(n.large, g.large))
		g.ind = rep(c(1:n.hub), g.list)
		for (i in 1:n.hub) {
			tmp = which(g.ind == i)
			theta.hub[tmp[1], tmp] = 1
			theta.hub[tmp, tmp[1]] = 1
			rm(tmp)
			gc()
		}
		
		B.list <- list()
		if(n.hub.diff==0){
			B.list[[1]] <- B.list[[2]] <- theta.hub*magn.nz.common
		}else{
			B1 <- B2 <- matrix(0,p,p)
			tmp <- which(g.ind%in%1:n.hub.diff)
			
			B1[-tmp,-tmp] <- theta.hub[-tmp,-tmp]*magn.nz.common
			B.list[[1]] <- B1
			
			B2[-tmp,-tmp] <- theta.hub[-tmp,-tmp]*magn.nz.common
			B2[tmp,tmp] <- theta.hub[tmp,tmp]*magn.nz.diff
			B.list[[2]] <- B2
		}
		
	}
	
	####compute (positive definite) concentration matrices
	SigInv <- list()
	for (k in 1:2){
		siginv <- B.list[[k]]
		e.min <- min(eigen(siginv)$values)
		siginv <- siginv+diag(abs(e.min)+emin+magn.diag,p)
		SigInv[[k]] <- siginv
		if(verbose){
			cat('ev.min: ',min(eigen(siginv)$values),'\n')
			cat('condnum: ',max(abs(eigen(siginv)$values))/min(abs(eigen(siginv)$values)),'\n')
		}
	}
	return(list(invcov1=SigInv[[1]],invcov2=SigInv[[2]]))
}