# Code for plotting networks
library(ggplot2)
library(network)
library(mclust)
library(mvtnorm)


# Suppress check warning
if(getRversion() >= '2.15.1')  
	utils::globalVariables(c('P.Corr', 'Edge', 'Mean', 'Type', 'Node.1', 
													 'Node.2', 'Group', 'Edge.Name', 'Corr'))


#' Build up dataframe for plotting dot plot with ggplot2
#' 
#' Internal function
#' 
#' @param net.clustering Clustering
#' @param cluster.names Cluster names
#' @param node.names Node names
#' @keywords internal
#' 
buildDotPlotDataFrame <- function(net.clustering, cluster.names, node.names) {
	results.frame = data.frame()
	
	for(cluster.i in 1:length(cluster.names)) {
		cluster.name = as.character(cluster.names[cluster.i])
		
		# Calculate partial correlations
		p.corrs = invcov2parcor(net.clustering$SigInv[,,cluster.i])
		
		Mu = net.clustering$Mu[,cluster.i]
		
		mean.edges = c()
		p.corr.edges = c()
		edge.names = c()
		
		# Record pcorr and means for all edges
		for(node.i in 1:(length(node.names) - 1)) {
			for(node.j in (node.i+1):length(node.names)) {
				
				if(abs(p.corrs[node.i, node.j]) >= 0) {
					mean.edges = c(mean.edges, min(Mu[node.i], Mu[node.j]))
					
					p.corr.edges = c(p.corr.edges, p.corrs[node.i, node.j])
					
					edge.name = paste(node.names[node.i], 
														node.names[node.j], sep='-')
					edge.names = c(edge.names, edge.name)
				}
			}
		}
		
		results.frame = rbind(results.frame, 
													data.frame(P.Corr=p.corr.edges,
																		 Mean=mean.edges,
																		 Edge=edge.names,
																		 Type=cluster.name))
	}
	
	return(results.frame)
	
}


#' Create a plot showing the edges with the highest partial correlation in any cluster.
#' 
#' This function takes the output of \code{\link{het_cv_glasso}} or 
#' \code{\link{mixglasso}} and creates a plot of the highest scoring edges along the
#' y axis, where, the edge in each cluster is represented by a circle whose area
#' is proportional to the smallest mean of the two nodes that make up the edge,
#' and the position along the y axis shows the partial correlation of the edge.
#'
#' @param net.clustering A network clustering object as returned by 
#' \code{\link{het_cv_glasso}} or \code{\link{mixglasso}}.
#' @param p.corrs.thresh Cutoff for the partial correlations; only edges with absolute 
#' partial correlation > p.corrs.thresh (in any cluster) will be displayed. 
#' @param hard.limit Additional hard limit on the number of edges to display. If 
#' p.corrs.thresh results in more edges than hard.limit, only hard.limit edges with the
#' highest partial correlation are returned. 
#' @param display If TRUE, print the plot to the current output device.
#' @param node.names Names for the nodes in the network.
#' @param group.names Names for the clusters or groups.
#' @param dot.size.range Graphical parameter for scaling the size of the circles (dots)
#' representing an edge in each cluster.
#' @export
#' @import ggplot2
#' @return Returns a ggplot2 object. If display=TRUE, additionally displays the 
#' plot.
#' @examples
#' n = 500
#' p = 10
#' s = 0.9
#' n.comp = 3
#'
#' # Create different mean vectors
#' Mu = matrix(0,p,n.comp)
#'
#' # Define non-zero means in each group (non-overlapping)
#' nonzero.mean = split(sample(1:p),rep(1:n.comp,length=p))
#'
#' # Set non-zero means to fixed value
#' for(k in 1:n.comp){
#' 	Mu[nonzero.mean[[k]],k] = -2/sqrt(ceiling(p/n.comp))
#' }
#'
#' # Generate data
#' sim.result = sim_mix_networks(n, p, n.comp, s, Mu=Mu)
#' mixglasso.result = mixglasso(sim.result$data, n.comp=3)
#' mixglasso.clustering = mixglasso.result$models[[mixglasso.result$bic.opt]]
#' 
#' dot_plot(mixglasso.clustering, p.corrs.thresh=0.5)
dot_plot <- function(net.clustering, p.corrs.thresh=0.25, hard.limit=50,
										display=TRUE, node.names=rownames(net.clustering$Mu),
										group.names=sort(unique(net.clustering$comp)),
										dot.size.range=c(3,12)) {
	
	if(is.null(node.names)) node.names = 1:dim(net.clustering$Mu)[1]
	
  cluster.assignments = net.clustering$comp
	
  results.frame = buildDotPlotDataFrame(net.clustering, group.names, 
  																			node.names)
  
  # Determine which edges to include
	biggest.pcorr = results.frame[abs(results.frame$P.Corr) > p.corrs.thresh,]  
  
  # If any edges above threshold
  if(dim(biggest.pcorr)[1] > 0) {
  	
  	interesting.edges = data.frame()
  	
  	# Find edges with maximum pcorr
  	for(edge in unique(biggest.pcorr$Edge)) {
  		edge.pcorr = biggest.pcorr[biggest.pcorr$Edge==edge,];
  		edge.max = edge.pcorr[which.max(abs(edge.pcorr$P.Corr)),];
  		interesting.edges = rbind(interesting.edges, edge.max)
  	}
  	  
    # Sorted by max p.corr
    interesting.edges = 
  	  interesting.edges[sort(abs(interesting.edges$P.Corr), index.return=TRUE)$ix,]
  
    # Enforce hard limit
    if(dim(interesting.edges)[1] > hard.limit) {
  	  interesting.edges = interesting.edges[1:hard.limit,]
    }
	
	  results.frame.reduced = results.frame[results.frame$Edge %in% interesting.edges$Edge,]
	  results.frame.reduced$Abs.P.Corr = abs(results.frame.reduced$P.Corr)
	
    # Ensure right order
	  results.frame.reduced$Edge = factor(results.frame.reduced$Edge, levels=interesting.edges$Edge)
  	
	  g <-  ggplot(results.frame.reduced, aes(x=P.Corr, y=Edge, size=Mean, colour=Type)) +
		  geom_point() +
  		scale_size(range=dot.size.range) +
		  theme(axis.title.y=element_blank())
	
	  if(display) print(g)
	
	  return(g)
	  
  } else {
  	warning('No edges with partial correlation above threshold.')
  	return(NULL)
  }
}


#' Create a scatterplot showing correlation between specific nodes in the network
#' for each pre-specified group.
#' 
#' This function takes the output of \code{\link{het_cv_glasso}} or 
#' \code{\link{mixglasso}} and creates a plot showing the correlation between specified
#' node pairs in the network for all groups. The subplots for each node pair are 
#' arranged in a numPairs by numGroups grid. Partial correlations associated 
#' with each node pair are also displayed.
#'
#' @param net.clustering A network clustering object as returned by 
#' \code{\link{het_cv_glasso}} or \code{\link{mixglasso}}.
#' @param data Observed data for the nodes, a numObs by numNodes matrix. Note 
#' that nodes need to be in the same ordering as in node.names.
#' @param node.pairs A matrix of size numPairs by 2, where each row contains a 
#' pair of nodes to display. If node.names is specified, names in node.pairs
#' must correspond to elements of node.names.
#' @param display If TRUE, print the plot to the current output device.
#' @param node.names Names for the nodes in the network. If NULL, names from 
#' net.clustering will be used.
#' @param group.names Names for the clusters or groups. If NULL, names from 
#' net.clustering will be used (by default these are integets 1:numClusters).
#' @param cex Scale factor for text and symbols in plot.
#' @export
#' @import ggplot2
#' @return Returns a ggplot2 object. If display=TRUE, additionally displays the 
#' plot.
#' @examples
#' n = 500
#' p = 10
#' s = 0.9
#' n.comp = 3
#'
#' # Create different mean vectors
#' Mu = matrix(0,p,n.comp)
#'
#' # Define non-zero means in each group (non-overlapping)
#' nonzero.mean = split(sample(1:p),rep(1:n.comp,length=p))
#'
#' # Set non-zero means to fixed value
#' for(k in 1:n.comp){
#' 	Mu[nonzero.mean[[k]],k] = -2/sqrt(ceiling(p/n.comp))
#' }
#'
#' # Generate data
#' sim.result = sim_mix_networks(n, p, n.comp, s, Mu=Mu)
#' mixglasso.result = mixglasso(sim.result$data, n.comp=3)
#' mixglasso.clustering = mixglasso.result$models[[mixglasso.result$bic.opt]]
#' 
#' # Specify edges
#' node.pairs = rbind(c(1,3), c(6,9),c(7,8))
#' 
#' # Create scatter plots of specified edges
#' scatter_plot(mixglasso.clustering, data=sim.result$data,
#'						 node.pairs=node.pairs)
scatter_plot <- function(net.clustering, data, node.pairs, display=TRUE, 
												node.names=rownames(net.clustering$Mu),
												group.names=sort(unique(net.clustering$comp)),
												cex=1) {
	
	if(class(net.clustering) != 'nethetclustering') 
		 stop('net.clustering needs to be an object of class nethetclustering')
	
	if(is.null(node.names)) node.names = 1:dim(net.clustering$Mu)[1]
	
	# Prepare data frame for plotting
	plotting.frame = data.frame()
	corr.frame = data.frame()
	
	num.samples = dim(data)[1]
	
	# Obtain partial correlations from graphical lasso output
	p.corrs = invcov2parcor_array(net.clustering$SigInv)
	
	# Cycle over node pairs/edges
	for(edge.i in 1:dim(node.pairs)[1]) {
		node.1 = node.pairs[edge.i, 1]
		node.2 = node.pairs[edge.i, 2]
		edge.name = paste(node.1, node.2, sep='\n')
					
		# Partial correlation
		p.corrs.values = p.corrs[node.names == node.1, 
												   node.names == node.2,]
		
		names(p.corrs.values) = group.names
		
		# Spearman correlation
		corrs = sapply(group.names, 
									 function(x) cor(data[net.clustering$comp == x,node.names==node.1],
		 														   data[net.clustering$comp == x,node.names==node.2]))
		
		corr.frame = rbind(corr.frame,
											 data.frame(Edge.Name=rep(edge.name, length(group.names)), 
											 					 Group=factor(group.names),
											 					 Corr=round(corrs, digits=2),
											 					 P.Corr=round(p.corrs.values, digits=2)))
		
		temp.frame = data.frame(Edge.Name=rep(edge.name, num.samples), 
														Node.1=data[,node.1], 
														Node.2=data[,node.2], 
														Group=factor(net.clustering$comp))
		
		plotting.frame = rbind(temp.frame, plotting.frame)
	}
	
	p <- ggplot(plotting.frame, aes(Node.1, Node.2)) +
		geom_point(aes(colour=Group), size=2*cex) +		
		geom_abline(linetype=3) + 
		facet_grid(Edge.Name ~ Group) +
		geom_text(data=corr.frame, aes(x=Inf, y=Inf, label=paste('PCorr: ', P.Corr, '\n', 'Corr: ', Corr, sep='')), size=6*cex, family='Times', hjust=1, vjust=1) +
		theme(strip.text.y = element_text(size=cex*15), strip.text.x = element_text(size=cex*20),
					axis.text=element_text(size=cex*18), axis.title=element_text(size=cex*24,face="bold"))
	
	if(display) print(p)
	
	return(p)
}

#' Plot networks
#' 
#' This function takes the output of \code{\link{screen_cv.glasso}} or 
#' \code{\link{mixglasso}} and creates a network plot using the network library.
#' 
#'
#' @param x A network clustering object as returned by 
#' \code{\link{screen_cv.glasso}} or \code{\link{mixglasso}}.
#' @param node.names Names for the nodes in the network. If NULL, names from 
#' net.clustering will be used.
#' @param group.names Names for the clusters or groups. If NULL, names from 
#' net.clustering will be used (by default these are integets 1:numClusters).
#' @param p.corrs.thresh Threshold applied to the absolute partial correlations. 
#' Edges that are below the threshold in all of the groups are not displayed.
#' @param print.pdf If TRUE, save the output as a PDF file.
#' @param pdf.filename If \code{print.pdf} is TRUE, specifies the file name of
#' the output PDF file.
#' @param ... Further arguments
#' @return Returns NULL and prints out the networks (or saves them to pdf if
#' \code{print.pdf} is TRUE. The networks are displayed as a series of nComps+1
#' plots, where in the first plot edge widths are shown according to 
#' the maximum partial correlation of the edge over all groups. The following plots
#' show the edges for each group. Positive partial correlation edges are shown in
#' black, negative ones in blue. If an edge is below the threshold on the absolute
#' partial correlation, it is displayed in gray or light blue respectively.
#' @method plot nethetclustering
#' @export 
#' @import network
plot.nethetclustering <- function(x, 
																	node.names=rownames(net.clustering$Mu),
																	group.names=sort(unique(net.clustering$comp)),
																	p.corrs.thresh=0.2, print.pdf=FALSE, 
																	pdf.filename='networks', ...){
  
	net.clustering = x
	
	if(is.null(node.names)) node.names = 1:dim(net.clustering$Mu)[1]
	
	p.corrs = invcov2parcor_array(net.clustering$SigInv)
	
	k = ncol(net.clustering$Mu)
	
	#combine parcors for all k states by taking max parcor values:
	PCmax = apply(abs(p.corrs), 1:2, max, na.rm=TRUE)
	diag(PCmax) = NA
	
	#plot network:
	if(print.pdf){
		pdf(paste(pdf.filename,'.pdf',sep=''), width=14, height=14, onefile=TRUE)
	}
	
	#create and plot network object for combined states:
	NWbin = PCmax>p.corrs.thresh
	diag(NWbin) = FALSE
	Network = network(NWbin)
	edgeWeight = ((PCmax-p.corrs.thresh)/max(PCmax-p.corrs.thresh,na.rm=TRUE)+0.1)*15
	
	Coord = plot.network(Network, 
											main="combined states",
											#coord=as.matrix(L[,2:3]),
											vertex.col="orange", 
											vertex.cex=1.2,#0.5-min(mean),
											vertex.rot=360/16, 
											vertex.border=grey(0.8),
											displaylabels=TRUE,
											boxed.labels=FALSE,
											label=node.names,
											label.pos=5,
											label.cex=0.7,
											edge.lwd=edgeWeight, 
											usearrows=FALSE, 
											pad=0.02)
	
	edge.sum = 0
	
	#now the same for each of the k states:
	for(i in 1:k) { 
		PC = abs(p.corrs[,,i])
	  NWbin = abs(PC)>0 & PCmax>p.corrs.thresh
	  diag(NWbin) = FALSE
	  Network = network(NWbin) 
	  #edgeWeight<-((PC-p.corrs.thresh)/max(PCmax-p.corrs.thresh,na.rm=TRUE)+0.1)*15
	  #cat('\n', group.names[i], '\n')
	  PC.temp = PC
	  diag(PC.temp) = 0
	  PC.temp[lower.tri(PC.temp)] = 0
	  indices = which(PC.temp >= p.corrs.thresh, arr.ind=TRUE)
	  #cat(paste(node.names[indices[,1]], node.names[indices[,2]], sep='-'), sep='\n')
	  edge.sum = edge.sum + sum(PC.temp >= p.corrs.thresh)
	 
 	  weights.temp = 6*(PC/max(PCmax,na.rm=TRUE))
	  #edgeWeight = (weights.temp + 0.1)*(15/6)
	  edgeWeight = (((1/(1+exp(-weights.temp))) - 0.5)*1.8 + 0.1)*12
	 
	  edge.col = matrix('black', dim(PC)[1], dim(PC)[2])
	  edge.col[p.corrs[,,i]<0] = 'blue'
	  edge.col[p.corrs[,,i] < p.corrs.thresh & p.corrs[,,i] > 0] = 'grey'
	  edge.col[p.corrs[,,i] > -p.corrs.thresh & p.corrs[,,i] < 0] = 'skyblue'
	 
	  plot.network(Network, 
	 						 main=group.names[i],
	 						 coord=Coord,
	 						 vertex.col="orange", 
	 						 vertex.cex=1.2,#0.5+mean[,i]-min(mean),
	 						 vertex.rot=360/16, 
	 						 vertex.border=grey(0.8),
	 						 displaylabels=TRUE,
	 						 boxed.labels=FALSE,
	 						 label=node.names,
	 						 label.pos=5,
	 						 label.cex=0.7,
	 						 #label.cex=0.3*(2^mean[,i]),
	 						 edge.lwd=edgeWeight, 
	 						 usearrows=FALSE, 
	 						 pad=0.02,
	 						 edge.col=edge.col)	 
	 
	  legend('bottomleft', legend=c('Positive PCorr',
	 										 'Negative PCorr',
	 										 'Pos PCorr < Threshold',
	 										 'Neg PCorr > -Threshold'), 
	 			 col=c('black', 'blue', 'grey', 'skyblue'),
	 			 lty=1, cex=0.75,
	 			 lwd=3,
	 			 xpd=TRUE)
	 
	}
	if(print.pdf){
		dev.off() #close PDF
	}
	
}



##' Plot two networks (GGMs)
##'
##' 
##' @title Plot two networks (GGMs)
##' @param invcov1 Inverse covariance matrix of GGM1.
##' @param invcov2 Inverse covariance matrix of GGM2.
##' @param node.label Names of nodes.
##' @param main Vector (two elements) with network names.
##' @param ... Other arguments (see plot.network).
##' @return Figure with two panels (for each network).
##' @author nicolas
##' @export
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
plot_2networks <- function(invcov1,invcov2,
													 node.label=paste('X',1:nrow(invcov1),sep=''),
													 main=c('',''),...){
	par.orig <- par(no.readonly=TRUE)
	par(mfrow=c(1,2),mar=c(1,1,1,1))
	adj1 <- invcov1!=0
	adj2 <- invcov2!=0
	adj <- list(adj1,adj2)
	nonzeros <- sapply(adj,sum)
	
	if(nonzeros[1]!=nonzeros[2]){
		coord <- plot.network(network(adj[[which.max(nonzeros)]]),
													main=main[which.max(nonzeros)],
													displaylabels=TRUE,
													label=node.label,
													usearrows=FALSE,...)
		plot.network(network(adj[[which.min(nonzeros)]]),
								 coord=coord,
								 main=main[which.min(nonzeros)],
								 displaylabels=TRUE,
								 label=node.label,
								 usearrows=FALSE,...)
	}
	
	if(nonzeros[1]==nonzeros[2]){
		coord <- plot.network(network(adj[[1]]),
													main=main[1],
													displaylabels=TRUE,
													label=node.label,
													usearrows=FALSE,...)
		plot.network(network(adj[[2]]),
								 coord=coord,
								 main=main[2],
								 displaylabels=TRUE,
								 label=node.label,
								 usearrows=FALSE,...)
	}
	par(par.orig)
}

##' Plotting function for object of class 'diffnet' 
##'
##' 
##' @title Plotting function for object of class 'diffnet' 
##' @param x object of class 'diffnet'
##' @param ... Further arguments.
##' @return Histogram over multi-split p-values.
##' @author nicolas
##' @method plot diffnet
##' @export 
plot.diffnet <- function(x,...){
	#if(is.null(x$medwi)){
	hh <- hist(x$ms.pval,
						 main='histogram single-split p-values',xlab='p-values',ylab='frequency',...)
	abline(v=x$medagg.pval,lty=2,col='red')
	abline(v=x$meinshagg.pval,lty=2,col='green')
	legend(x=min(hh$mids),y=max(hh$counts),lty=c(2,2),col=c('red','green'),legend=c('median aggregated','meinshausen aggregated'))
	#}else{
	#    medwi <- x$medwi
	#    k <- ncol(x$medwi[[1]])
	#    par(mfrow=c(2,2))
	#    image(x=1:k,y=1:k,abs(medwi$modIpop1),xlab='',ylab='',main='median invcov1')
	#    image(x=1:k,y=1:k,abs(medwi$modIpop2),xlab='',ylab='',main='median invcov2')
	#    hist(x$ms.pval,breaks=10,
	#         main='histogram single-split p-values',xlab='p-values',ylab='frequency')
	#    abline(v=x$medagg.pval,lty=2,col='red')
	#}
}

##' Plotting function for object of class 'diffregr' 
##'
##' 
##' @title Plotting function for object of class 'diffregr' 
##' @param x object of class 'diffregr'
##' @param ... Further arguments.
##' @return Histogram over multi-split p-values.
##' @author nicolas
##' @method plot diffregr
##' @export
plot.diffregr <- function(x,...){
	hh <- hist(x$ms.pval,
						 main='histogram single-split p-values',xlab='p-values',ylab='frequency',...)
	abline(v=x$medagg.pval,lty=2,col='red')
	abline(v=x$meinshagg.pval,lty=2,col='green')
	legend(x=min(hh$mids),y=max(hh$counts),lty=c(2,2),col=c('red','green'),legend=c('median aggregated','meinshausen aggregated'))
	
}

##' Plotting function for object of class 'ggmgsa' 
##'
##' 
##' @title Plotting function for object of class 'ggmgmsa' 
##' @param x object of class 'ggmgsa'
##' @param ... Further arguments.
##' @return Boxplot of single-split p-values.
##' @author nicolas
##' @method plot ggmgsa
##' @export 
plot.ggmgsa <- function(x,...){
	boxplot(t(x$pval),names=x$gs.names,xlab='gene-sets',ylab='single-split p-values (uncorrected)',...)
}