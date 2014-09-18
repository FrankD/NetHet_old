# Code for plotting networks
require(ggplot2)
require(network)

#' Build up dataframe for plotting dot plot with ggplot2
#' 
#'
#' Internal function
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
#' This function takes the output of {\link{screen_cv.glasso} or 
#' {\link{mixglasso} and creates a plot of the highest scoring edges along the
#' y axis, where, the edge in each cluster is represented by a circle whose area
#' is proportional to the smallest mean of the two nodes that make up the edge,
#' and the position along the y axis shows the partial correlation of the edge.
#'
#' @param net.clustering A network clustering object as returned by {\link{screen_cv.glasso} or 
#' {\link{mixglasso}.
#' @param pcor.cutoff Cutoff for the partial correlations; only edges with absolute 
#' partial correlation > pcor.cutoff (in any cluster) will be displayed. 
#' @param hard.limit Additional hard limit on the number of edges to display. If 
#' pcor.cutoff results in more edges than hard.limit, only hard.limit edges with the
#' highest partial correlation are returned. 
#' @param display If TRUE, print the plot to the current output device.
#' @param node.names Names for the nodes in the network.
#' @param group.names Names for the clusters or groups.
#' @param dot.size.range Graphical parameter for scaling the size of the circles (dots)
#' representing an edge in each cluster.
#' @export
#' @return Returns a ggplot2 object. If display=TRUE, additionally displays the 
#' plot.
dotPlot <- function(net.clustering, pcor.cutoff=0.25, hard.limit=50,
										display=TRUE, node.names=rownames(net.clustering$Mu),
										group.names=sort(unique(net.clustering$comp)),
										dot.size.range=c(3,12)) {
	
	if(is.null(node.names)) node.names = 1:dim(net.clustering$Mu)[1]
	
  cluster.assignments = net.clustering$comp
	
  results.frame = buildDotPlotDataFrame(net.clustering, group.names, 
  																			node.names)
  
  # Determine which edges to include
	biggest.pcorr = results.frame[abs(results.frame$P.Corr) > pcor.cutoff,]  
  
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
#' This function takes the output of {\link{screen_cv.glasso} or 
#' {\link{mixglasso} and creates a plot showing the correlation between specified
#' node pairs in the network for all groups. The subplots for each node pair are 
#' arranged in a numPairs by numGroups grid. Partial correlations associated 
#' with each node pair are also displayed.
#'
#' @param net.clustering A network clustering object as returned by {\link{screen_cv.glasso} or 
#' {\link{mixglasso}.
#' @param data Observed data for the nodes, a numObs by numNodes matrix. Note 
#' that nodes need to be in the same ordering as in node.names.
#' @param nodes.pairs A matrix of size numPairs by 2, where each row contains a 
#' pair of nodes to display. If node.names is specified, names in node.pairs
#' must correspond to elements of node.names.
#' @param display If TRUE, print the plot to the current output device.
#' @param node.names Names for the nodes in the network. If NULL, names from 
#' net.clustering will be used.
#' @param group.names Names for the clusters or groups. If NULL, names from 
#' net.clustering will be used (by default these are integets 1:numClusters).
#' @export
#' @return Returns a ggplot2 object. If display=TRUE, additionally displays the 
#' plot.
#' 
scatterPlot <- function(net.clustering, data, node.pairs, display=TRUE, 
												node.names=rownames(net.clustering$Mu),
												group.names=net.clustering$group.names) {
	
	if(class(net.clustering) != 'nethetclustering') 
		 stop('net.clustering needs to be an object of class nethetclustering')
	
	if(is.null(node.names)) node.names = 1:dim(net.clustering$Mu)[1]
	
	# Prepare data frame for plotting
	plotting.frame = data.frame()
	corr.frame = data.frame()
	
	num.samples = dim(data)[1]
	
	# Obtain partial correlations from graphical lasso output
	p.corrs = sapply(1:length(group.names), 
				 function(i) invcov2parcor(net.clustering$SigInv[,,i]), simplify='array')
	
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
		geom_point(aes(colour=Group)) +		
		geom_abline(linetype=3) + 
		facet_grid(Edge.Name ~ Group) +
		geom_text(data=corr.frame, aes(x=Inf, y=Inf, label=paste('PCorr: ', P.Corr, '\n', 'Corr: ', Corr, sep='')), size=6, family='Times', hjust=1, vjust=1) +
		theme(strip.text.y = element_text(size=15), strip.text.x = element_text(size=20),
					axis.text=element_text(size=18), axis.title=element_text(size=24,face="bold"))
	
	if(display) print(p)
	
	return(p)
}

#' Plot networks
#' 
#' This function takes the output of {\link{screen_cv.glasso} or 
#' {\link{mixglasso} and creates a network plot using the network library.
#' 
#'
#' @param net.clustering A network clustering object as returned by {\link{screen_cv.glasso} or 
#' {\link{mixglasso}.
#' @param node.names Names for the nodes in the network. If NULL, names from 
#' net.clustering will be used.
#' @param group.names Names for the clusters or groups. If NULL, names from 
#' net.clustering will be used (by default these are integets 1:numClusters).
#' @param p.corr.thresh Threshold applied to the absolute partial correlations. 
#' Edges that are below the threshold in all of the groups are not displayed.
#' @param print.pdf If TRUE, save the output as a PDF file.
#' @param pdf.filename If \code{print.pdf} is TRUE, specifies the file name of
#' the output PDF file.
#' @export
#' @return Returns NULL and prints out the networks (or saves them to pdf if
#' \code{print.pdf} is TRUE. The networks are displayed as a series of nComps+1
#' plots, where in the first plot edge widths are shown according to 
#' the maximum partial correlation of the edge over all groups. The following plots
#' show the edges for each group. Positive partial correlation edges are shown in
#' black, negative ones in blue. If an edge is below the threshold on the absolute
#' partial correlation, it is displayed in gray or light blue respectively.
#' 
plot.nethetclustering <- function(net.clustering, 
																	node.names=rownames(net.clustering$Mu),
																	group.names=sort(unique(net.clustering$comp)),
																	p.corrs.thresh=0.2, print.pdf=FALSE, 
																	pdf.filename='networks'){
  
	if(is.null(node.names)) node.names = 1:dim(net.clustering$Mu)[1]
	
	p.corrs = sapply(1:length(group.names), 
									function(i) invcov2parcor(net.clustering$SigInv[,,i]), 
									                          simplify='array')
	
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
