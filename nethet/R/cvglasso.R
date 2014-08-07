##Packages
library(mvtnorm)
library(glasso)
library(huge)

#############################
##-------Screening---------##
#############################


#' Convert inverse covariance to partial correlation
#' 
#'  
#' @param invcov Inverse covariance matrix
#' @export
#' @return The partial correlation matrix.
#' 
invcov2parcor <- function(invcov){
  return(-invcov*tcrossprod(sqrt(1/diag(invcov))))
}

##Lambda-grid
lambdagrid_mult <- function(lambda.min,lambda.max,nr.gridpoints){
    mult.const <- (lambda.max/lambda.min)^(1/(nr.gridpoints-1))
    return(lambda.min*mult.const^((nr.gridpoints-1):0))
}
lambdagrid_lin <- function(lambda.min,lambda.max,nr.gridpoints){
    mult.const <- (lambda.max/lambda.min)^(1/(nr.gridpoints-1))
    return(seq(lambda.max,lambda.min,length=nr.gridpoints))
}

##Make-grid
make_grid <- function(lambda.min,lambda.max,nr.gridpoints,method='lambdagrid_mult'){
  eval(as.name(method))(lambda.min,lambda.max,nr.gridpoints)
}

##Lambda-max
lambda.max <- function(x){
  n <- nrow(x)
  s.var <- var(x)
  diag(s.var) <- 0
  return(n*max(abs(s.var))/2)
}

##Crossvalidation plot
error.bars <- function (x, upper, lower, width = 0.02, ...)
{
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}
plotCV <- function(lambda,cv,cv.error,se=TRUE,type='b',...){
  if (se==TRUE){ylim <- range(cv,cv+cv.error,cv-cv.error)}
  else {ylim <- range(cv)}
  plot(lambda,cv,type=type,ylim=ylim,...)
  if (se)
    error.bars(lambda,cv+cv.error,cv-cv.error,width=1/length(lambda))
  invisible()
}
cv.fold <- function(n,folds=10){
  split(sample(1:n),rep(1:folds,length=n))
}

##Hugepath
hugepath <- function(s,rholist,penalize.diagonal=NULL,trace=NULL){
  #fit.huge <- huge(s,method = "glasso",cov.output =TRUE,verbose = FALSE)
  fit.huge <- huge(s,lambda=sort(rholist,decreasing=TRUE),method = "glasso",cov.output =TRUE,verbose = FALSE)
  wi <- sapply(fit.huge$icov,as.matrix,simplify='array')
  w <- sapply(fit.huge$cov,as.matrix,simplify='array')
  #return(list(wi=wi[,,length(fit.huge$lambda):1],w=w[,,length(fit.huge$lambda):1]))
  return(list(rholist=rholist,wi=wi[,,length(rholist):1],w=w[,,length(rholist):1]))
}

##Truncating concentration matrix
mytrunc.method <- function(n,wi,method='linear.growth',trunc.k=5){

  p <- ncol(wi)
  if (method=='none'){
    return(list(wi=wi))
  }
  if(method=='linear.growth'){
    wi.trunc <- wi
    diag(wi.trunc) <- 0
    nonzero <- min(2*ceiling(p*ceiling(n/trunc.k)/2),sum(wi.trunc!=0))
    wi.trunc[-order(abs(wi.trunc),decreasing=TRUE)[1:nonzero]] <- 0
    diag(wi.trunc) <- diag(wi)
    return(list(wi=wi.trunc))
  }
  if(method=='sqrt.growth'){
    wi.trunc <- wi
    diag(wi.trunc) <- 0
    nonzero <- min(2*ceiling(p*sqrt(n)/2),sum(wi.trunc!=0))
    wi.trunc[-order(abs(wi.trunc),decreasing=TRUE)[1:nonzero]] <- 0
    diag(wi.trunc) <- diag(wi)
    return(list(wi=wi.trunc))
  }
}

#' Cross-validated glasso on single dataset
#' 
#'
#' Run glasso on a single dataset, using cross-validation to estimate the
#' penalty parameter lambda.
#'
#'  
#' @param x The input data. Needs to be a num.samples by dim.samples matrix.
#' @param include.mean 
#' @param fold Number of folds in the cross-validation (default=10).
#' @param length.lambda Length of lambda path to consider (default=20).
#' @param lambdamin.ratio 
#' @param penalize.diagonal If TRUE apply penalization to diagonal of inverse
#' covariance as well. (default=FALSE)
#' @param trunc.method
#' @param trunc.k
#' @param plot.it
#' @param se
#' @param use.package 'glasso' (default) or 'huge'
#' @param verbose If TRUE, output progress.
#' @export
#' @return Returns a list with named elements 'rho.opt', 'wi', 'wi.orig', 'mu', 
#' Variable rho.opt is the scaled penalization parameter (?). 
#' The variables wi and wi.orig are matrices of size dim.samples by dim.samples 
#' containing the truncated and untruncated inverse covariance matrix. Variable 
#' Mu mean of the input data.
#' 
screen_cv.glasso <- function(x,include.mean=TRUE,folds=10,length.lambda=20,
                             lambdamin.ratio=ifelse(ncol(x)>nrow(x),0.01,0.001),penalize.diagonal=FALSE,
                             trunc.method='none',trunc.k=5,plot.it=FALSE,se=FALSE,use.package='glasso',verbose=TRUE)
{
  
  gridmax <- lambda.max(x)
  gridmin <- lambdamin.ratio*gridmax
  lambda <- make_grid(gridmin,gridmax,length.lambda)[length.lambda:1]
  
  colnames(x)<-paste('x',1:ncol(x),sep='')  
  all.folds <- cv.fold(nrow(x),folds)
  residmat <- matrix(NA,folds,length(lambda))
  for (cvfold in 1:folds){
    omit <- all.folds[[cvfold]]
    s <- var(x[-omit,])
    if (include.mean==TRUE){
      mu <- colMeans(x[-omit,,drop=FALSE])
    }else{
      mu <- rep(0,ncol(x))
    }
    fit.path <- eval(as.name(paste(use.package,'path',sep='')))(s,rholist=2*lambda/nrow(x),penalize.diagonal=penalize.diagonal,trace=0)
    if(length(omit)==1){
      residmat[cvfold,] <- -2*apply(fit.path$w,3,dmvnorm,log=TRUE,mean=mu,x=x[omit,,drop=FALSE])
    }else{
      residmat[cvfold,] <- apply(-2*apply(fit.path$w,3,dmvnorm,log=TRUE,mean=mu,x=x[omit,,drop=FALSE]),2,sum)
    }
  }
  cv <- apply(residmat,2,mean)
  cv.error <- sqrt(apply(residmat,2,var)/folds)
  gl.opt<-glasso(var(x),rho=2*lambda[which.min(cv)]/nrow(x),penalize.diagonal=penalize.diagonal)
  if(verbose){
    cat('la.min:',gridmin,'\n')
    cat('la.max:',gridmax,'\n')
    cat('la.opt:',lambda[which.min(cv)],'\n')
  }
  w<-gl.opt$w
  wi<-gl.opt$wi
  wi[abs(wi)<10^{-3}]<-0
  wi <- (wi+t(wi))/2
  colnames(w)<-rownames(w)<-colnames(wi)<-rownames(wi)<-colnames(x)
  wi.trunc <- mytrunc.method(n=nrow(x),wi=wi,method=trunc.method,trunc.k=trunc.k)$wi
  
  if (plot.it){
    plotCV(lambda,cv,cv.error,se=se)
  }
  list(rho.opt=2*lambda[which.min(cv)]/nrow(x),wi=wi.trunc,wi.orig=wi,mu=colMeans(x))
}

#' Cross-validated glasso on heterogeneous dataset with grouping
#' 
#'
#' Run glasso on a heterogeneous dataset to obtain networks (inverse covariance 
#' matrices) of the variables in the dataset for each pre-specified group of 
#' samples.
#'
#' This function uses code{\link{screen_cv.lasso}} to run glasso with 
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
#' @param ... Further parameters to be passed to code{\link{screen_cv.lasso}}.
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
het.cv.glasso <- function(data, grouping=rep(1, dim(data)[1]), mc.flag=FALSE,
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
                       verbose=verbose, ...)
  } else {
    results = lapply(data.list, screen_cv.glasso, use.package=use.package, 
                     verbose=verbose, ...)
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
  
  class(results.all) = 'nethetclustering'
  
  return(results.all)    
  
}