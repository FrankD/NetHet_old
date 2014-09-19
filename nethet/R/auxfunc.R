
library(network)

##' Convert inverse covariance to partial correlation
##' 
##'  
##' @param invcov Inverse covariance matrix
##' @export
##' @return The partial correlation matrix.
##' 
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
invcov2parcor_array <- function(invcov.array) {
	return(sapply(1:dim(invcov.array)[3], 
								function(i) invcov2parcor(invcov.array[,,i]), 
						    simplify='array'))
}

##' Export networks as a CSV table.
##'
##' This function takes the output of {\link{screen_cv.glasso} or 
##' {\link{mixglasso} and exports it as a text table in CSV format, where each
##' entry in the table records an edge in one group and its partial correlation.
##' 
##' @param net.clustering A network clustering object as returned by 
##' {\link{screen_cv.glasso} or {\link{mixglasso}.
##' @param file Filename to save the network table under.
##' @param node.names Names for the nodes in the network. If NULL, names from 
##' net.clustering will be used.
##' @param group.names Names for the clusters or groups. If NULL, names from 
##' net.clustering will be used (by default these are integets 1:numClusters).
##' @param p.corr.thresh Threshold applied to the absolute partial correlations. 
##' Edges that are below the threshold in all of the groups are not exported. 
##' Using a negative value will export all possible edges (including those 
##' with zero partial correlation).
##' @param ... Further parameters passed to \link{write.csv}.
##' @return NULL
##' @author Frank Dondelinger
##' @export
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