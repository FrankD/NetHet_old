
library(network)

##' Convert inverse covariance to partial correlation
##' 
##'  
##' @param invcov Inverse covariance matrix
##' @export
##' @return The partial correlation matrix.
##' @examples
##' inv.cov = generate_inv_cov(p=25)
##' p.corr = invcov2parcor(inv.cov)
invcov2parcor <- function(invcov){
	return(-invcov*tcrossprod(sqrt(1/diag(invcov))))
}

##' Convert inverse covariance to partial correlation for several inverse
##' covariance matrices collected in an array.
##' 
##' @param invcov.array Array of inverse covariance matrices, of dimension
##' numNodes by numNodes by numComps.
##' @export
##' @return Array of partial correlation matrices of dimension numNodes by 
##' numNodes by numComps
##' @examples
##' invcov.array = sapply(1:5, function(x) generate_inv_cov(p=25), simplify='array')
##' p.corr = invcov2parcor_array(invcov.array)
invcov2parcor_array <- function(invcov.array) {
	return(sapply(1:dim(invcov.array)[3], 
								function(i) invcov2parcor(invcov.array[,,i]), 
						    simplify='array'))
}

##' Export networks as a CSV table.
##'
##' This function takes the output of \code{\link{het_cv_glasso}} or 
##' \code{\link{mixglasso}} and exports it as a text table in CSV format, where each
##' entry in the table records an edge in one group and its partial correlation.
##' 
##' @param net.clustering A network clustering object as returned by 
##' \code{\link{screen_cv.glasso}} or \code{\link{mixglasso}}.
##' @param file Filename to save the network table under.
##' @param node.names Names for the nodes in the network. If NULL, names from 
##' net.clustering will be used.
##' @param group.names Names for the clusters or groups. If NULL, names from 
##' net.clustering will be used (by default these are integets 1:numClusters).
##' @param p.corrs.thresh Threshold applied to the absolute partial correlations. 
##' Edges that are below the threshold in all of the groups are not exported. 
##' Using a negative value will export all possible edges (including those 
##' with zero partial correlation).
##' @param ... Further parameters passed to \link{write.csv}.
##' @return NULL
##' @author Frank Dondelinger
##' @export
##' @examples
##' n = 500
##' p = 10
##' s = 0.9
##' n.comp = 3
##'
##' # Create different mean vectors
##' Mu = matrix(0,p,n.comp)
##'
##' # Define non-zero means in each group (non-overlapping)
##' nonzero.mean = split(sample(1:p),rep(1:n.comp,length=p))
##'
##' # Set non-zero means to fixed value
##' for(k in 1:n.comp){
##' 	Mu[nonzero.mean[[k]],k] = -2/sqrt(ceiling(p/n.comp))
##' }
##'
##' # Generate data
##' sim.result = sim_mix_networks(n, p, n.comp, s, Mu=Mu)
##' mixglasso.result = mixglasso(sim.result$data, n.comp=3)
##' mixglasso.clustering = mixglasso.result$models[[mixglasso.result$bic.opt]]
##' 
##' \dontrun{
##' # Save network in CSV format suitable for Cytoscape import
##' export_network(mixglasso.clustering, file='nethet_network.csv',
##'							 p.corrs.thresh=0.25, quote=FALSE)
##' }
export_network <- function(net.clustering, file='network_table.csv',
													 node.names=rownames(net.clustering$Mu),
													 group.names=sort(unique(net.clustering$comp)),
													 p.corrs.thresh=0.2, ...) {	
	
	if(is.null(node.names)) node.names = 1:dim(net.clustering$Mu)[1]
	
	# Calculate partial correlations
	p.corrs = invcov2parcor_array(net.clustering$SigInv)
	
	edge.table = data.frame()
	
	for(group.i in 1:length(group.names)) {
		group.name = group.names[group.i]
		
		p.corrs.group = p.corrs[,,group.i]
		
		diag(p.corrs.group) = NA # Discard diagonal
		p.corrs.group[lower.tri(p.corrs.group)] = NA # Ignore lower triangular matrix
		
		edges = which(abs(p.corrs.group) > p.corrs.thresh, 
									arr.ind=TRUE, useNames=TRUE)
		
		edge.table.temp = data.frame(Node.1=node.names[edges[,1]], 
																 Node2=node.names[edges[,2]],
																 P.Corr=c(abs(p.corrs.group[edges])), 
																 Edge.Sign=c(sign(p.corrs.group[edges])),
																 Group=rep(group.name, dim(edges)[1]))
		
		edge.table = rbind(edge.table, edge.table.temp)
	}
	
	write.csv(edge.table, file=file, ...)
}

##' Summary function for object of class 'diffnet'
##'
##' 
##' @title Summary function for object of class 'diffnet'
##' @param x object of class 'diffnet'
##' @return aggregated p-values
##' @author nicolas
##' @export
summary.diffnet <- function(x){
	out <- data.frame(medagg.pval=x$medagg.pval,meinshagg.pval=x$meinshagg.pval)
	rownames(out) <- 'aggregated p-values'
	return(out)
}

##' Summary function for object of class 'diffregr'
##'
##' 
##' @title Summary function for object of class 'diffregr'
##' @param x object of class 'diffregr'
##' @return aggregated p-values
##' @author nicolas
##' @export
summary.diffregr <- function(x){
	out <- data.frame(medagg.pval=x$medagg.pval,meinshagg.pval=x$meinshagg.pval)
	rownames(out) <- 'aggregated p-values'
	return(out)
}

##' Summary function for object of class 'ggmgsa'
##'
##' 
##' @title Summary function for object of class 'ggmgsa'
##' @param x object of class 'ggmgsa'
##' @return aggregated p-values
##' @author nicolas
##' @export
summary.ggmgsa <- function(x){
	out <- data.frame(medagg.pval=x$medagg.pval,meinshagg.pval=x$meinshagg.pval)
	rownames(out) <- x$gs.names
	return(out)
}

##' Summary function for object of class 'nethetclustering'
##'
##' 
##' @title Summary function for object of class 'nethetclustering'
##' @param x object of class 'nethetclustering'
##' @return Network statistics (a 'nethetsummary' object)
##' @author frankd
##' @export
summary.nethetclustering <- function(x) {
	
	out = list(mix.prob=x$mix.prob, p=dim(x$Mu)[1], n.comp=dim(x$Mu)[2],
			 loglik=x$loglik, bic=x$bic, mmdl=x$mmdl, compsize=x$compsize)
	
	class(out) = 'nethetsummary'
	return(out)
}

##' Print function for object of class 'nethetsummary'
##'
##' 
##' @title Print function for object of class 'nethetsummmary'
##' @param x object of class 'nethetsummary'
##' @return NULL
##' @author frankd
##' @export
##' 
print.nethetsummary <- function(x) {
	
	cat('Heterogeneous network mixture\n',
			'Number of nodes:', x$p, '\n',
			'Number of components:', x$n.comp, '\n',
			'Log Likelihood:', x$loglik, '\n',
			'BIC Score:', x$bic, 'MMDL Score:', x$mmdl, '\n',
		  'Mixture.probabilities:\n', x$mix.prob, '\n',
			'Component size:\n', x$compsize, '\n')
	
}