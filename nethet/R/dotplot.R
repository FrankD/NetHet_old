# Code for creating a dot plot of the strongest edges from mixglasso output

require(ggplot2)


# Build up dataframe for plotting dot plot with ggplot2
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

# Main function for creating dotplot of edges with the highest partial correlation
# (in any cluster).
dotPlot <- function(net.clustering, pcor.cutoff=0.25, hard.limit=50,
										display=TRUE, node.names=rownames(net.clustering$Mu),
										cluster.names=sort(unique(net.clustering$comp)),
										dot.size.range=c(3,12)) {
	
  cluster.assignments = net.clustering$comp
	
  results.frame = buildDotPlotDataFrame(net.clustering, cluster.names, 
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

	