% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/auxfunc.R
\name{export_network}
\alias{export_network}
\title{Export networks as a CSV table.}
\usage{
export_network(net.clustering, file = "network_table.csv",
  node.names = rownames(net.clustering$Mu),
  group.names = sort(unique(net.clustering$comp)), p.corrs.thresh = 0.2,
  ...)
}
\arguments{
\item{net.clustering}{A network clustering object as returned by
\code{\link{screen_cv.glasso}} or \code{\link{mixglasso}}.}

\item{file}{Filename to save the network table under.}

\item{node.names}{Names for the nodes in the network. If NULL, names from
net.clustering will be used.}

\item{group.names}{Names for the clusters or groups. If NULL, names from
net.clustering will be used (by default these are integets 1:numClusters).}

\item{p.corrs.thresh}{Threshold applied to the absolute partial correlations.
Edges that are below the threshold in all of the groups are not exported.
Using a negative value will export all possible edges (including those
with zero partial correlation).}

\item{...}{Further parameters passed to \link{write.csv}.}
}
\value{
Function does not return anything.
}
\description{
This function takes the output of \code{\link{het_cv_glasso}} or
\code{\link{mixglasso}} and exports it as a text table in CSV format, where each
entry in the table records an edge in one group and its partial correlation.
}
\examples{
n = 500
p = 10
s = 0.9
n.comp = 3

# Create different mean vectors
Mu = matrix(0,p,n.comp)

# Define non-zero means in each group (non-overlapping)
nonzero.mean = split(sample(1:p),rep(1:n.comp,length=p))

# Set non-zero means to fixed value
for(k in 1:n.comp){
	Mu[nonzero.mean[[k]],k] = -2/sqrt(ceiling(p/n.comp))
}

# Generate data
sim.result = sim_mix_networks(n, p, n.comp, s, Mu=Mu)
mixglasso.result = mixglasso(sim.result$data, n.comp=3)
mixglasso.clustering = mixglasso.result$models[[mixglasso.result$bic.opt]]

\dontrun{
# Save network in CSV format suitable for Cytoscape import
export_network(mixglasso.clustering, file='nethet_network.csv',
							 p.corrs.thresh=0.25, quote=FALSE)
}
}
\author{
Frank Dondelinger
}

