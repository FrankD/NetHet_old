##Packages
library(mvtnorm)
library(glasso)
library(huge)


#' Cross-validated glasso on heterogeneous dataset with grouping
#' 
#'
#' Run glasso on a heterogeneous dataset to obtain networks (inverse covariance 
#' matrices) of the variables in the dataset for each pre-specified group of 
#' samples.
#'
#' This function runs the graphical lasso with 
#' cross-validation to determine the best parameter lambda for each group of 
#' samples. Note that this function defaults to using package huge (rather than
#' package glasso) unless otherwise specified, as it tends to be more 
#' numerically stable.
#'
#'  
#' @param data The heterogenous network data. Needs to be 
#' a num.samples by dim.samples matrix or dataframe.
#' @param grouping The grouping of samples; a vector of length num.samples,
#' with num.groups unique elements.
#' @param mc.flag Whether to use parallel processing via package mclapply to
#' distribute the glasso estimation over different groups.
#' @param use.package 'glasso' for glasso package, or 'huge' for huge package 
#' (default)
#' @param normalise If TRUE, normalise the columns of the data matrix before 
#' running glasso.
#' @param verbose If TRUE, output progress.
#' @param ... Further parameters to be passed to \code{screen_cv.lasso}.
#' @export
#' @return Returns a list with named elements 'Sig', 'SigInv', 'Mu', 'Sigma.diag', 
#' 'group.names' and 'var.names. 
#' The variables Sig and SigInv are arrays of size dim.samples by dim.samples 
#' by num.groups, where the first two dimensions contain the (inverse)
#' covariance matrix for the network obtained by running glasso on group k. Variables 
#' Mu and Sigma.diag contain the mean and variance of the input data,
#' and group.names and var.names contains the names for the groups and
#' variables in the data (if specified as colnames of the input data matrix).
#' 
het_cv_glasso <- function(data, grouping=rep(1, dim(data)[1]), mc.flag=FALSE,
                          use.package='huge', normalise=FALSE, verbose=FALSE, ...) {
  
  group.names = sort(unique(grouping))
  
  mu = matrix(0, dim(data)[2], length(group.names))
  Sigma.diag = matrix(0, dim(data)[2], length(group.names))
  
  data.list = list()
  
  # Scale data if necessary and split into list
  for(group.i in 1:length(group.names)) {
    group.name = group.names[group.i]
    group.data = data[grouping==group.name,]
    
    scaled.group.data = scale(group.data)
    
    data.list[[group.i]] =  if(normalise) scaled.group.data
                            else group.data 
      
    mu[,group.i] = attributes(scaled.group.data)[['scaled:center']]
    Sigma.diag[,group.i] = attributes(scaled.group.data)[['scaled:scale']]
    
    rownames(mu) = colnames(group.data)
    colnames(mu) = group.names
    rownames(Sigma.diag) = colnames(group.data)
    colnames(Sigma.diag) = group.names
  }
  
  # Run glasso on each group
  if(mc.flag) {
    results = mclapply(data.list, screen_cv.glasso, use.package=use.package, 
                       verbose=verbose, include.mean=TRUE, trunc.method='none', ...)
  } else {
    results = lapply(data.list, screen_cv.glasso, use.package=use.package, 
                     verbose=verbose, include.mean=TRUE, trunc.method='none', ...)
  }
  
  
  # Collect results
  results.all = list()
  
  wi = array(0, dim=c(dim(data)[2], dim(data)[2], length(group.names)))
  w = array(0, dim=c(dim(data)[2], dim(data)[2], length(group.names)))
      
  for(group.i in 1:length(group.names)) {
    wi[,,group.i] = results[[group.i]]$wi
    w[,,group.i] = solve(wi[,,group.i])
  }
      
  results.all$SigInv = wi
  results.all$Sig = w
  results.all$Mu = mu
  results.all$Sigma.diag = Sigma.diag
  results.all$group.names = group.names
  results.all$proteins = colnames(data)
  results.all$comp = grouping # For consistency with mixglasso
  
  class(results.all) = 'nethetclustering'
  
  return(results.all)    
  
}
