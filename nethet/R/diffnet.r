###############################################################################
# High-Dimensional Two-Sample Testing for Gaussian Graphical Model
#------------------------------------------------------------------
#
#
###############################################################################

#####################
# Required Packages 
#####################
library(mvtnorm)
library(glasso)
library(parcor)
library(GeneNet)
library(huge)
library(CompQuadForm)
library(ggm)


#############################
##-------Screening---------##
#############################

##' Lambdamax
##'
##' 
##' @title Lambdamax
##' @param x no descr
##' @return no descr
##' @author n.stadler
lambda.max <- function(x){
  n <- nrow(x)
  s.var <- var(x)
  diag(s.var) <- 0
  return(n*max(abs(s.var))/2)
}

##' Lambda-grid (log scale)
##'
##' 
##' @title Lambda-grid
##' @param lambda.min no descr
##' @param lambda.max no descr
##' @param nr.gridpoints no descr
##' @return no descr
##' @author n.stadler
lambdagrid_mult <- function(lambda.min,lambda.max,nr.gridpoints){
    mult.const <- (lambda.max/lambda.min)^(1/(nr.gridpoints-1))
    return(lambda.min*mult.const^((nr.gridpoints-1):0))
}

##' Lambda-grid (linear scale)
##'
##' 
##' @title Lambda-grid
##' @param lambda.min no descr
##' @param lambda.max no descr
##' @param nr.gridpoints no descr
##' @return no descr
##' @author n.stadler
lambdagrid_lin <- function(lambda.min,lambda.max,nr.gridpoints){
    return(seq(lambda.max,lambda.min,length=nr.gridpoints))
}

##' Make grid
##'
##' 
##' @title Make grid
##' @param lambda.min no descr
##' @param lambda.max no descr
##' @param nr.gridpoints no descr
##' @param method no descr
##' @return no descr
##' @author n.stadler
make_grid <- function(lambda.min,lambda.max,nr.gridpoints,method='lambdagrid_mult'){
  eval(as.name(method))(lambda.min,lambda.max,nr.gridpoints)
}

##' Error bars for plotCV
##'
##' 
##' @title Error bars for plotCV
##' @param x no descr
##' @param upper no descr
##' @param lower no descr
##' @param width no descr
##' @param ... no descr
##' @return no descr
##' @author n.stadler
error.bars <- function (x, upper, lower, width = 0.02, ...)
{
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}

##' plotCV
##'
##' 
##' @title plotCV
##' @param lambda no descr
##' @param cv no descr
##' @param cv.error no descr
##' @param se no descr
##' @param type no descr
##' @param ... no descr
##' @return no descr
##' @author n.stadler
plotCV <- function(lambda,cv,cv.error,se=TRUE,type='b',...){
  if (se==TRUE){ylim <- range(cv,cv+cv.error,cv-cv.error)}
  else {ylim <- range(cv)}
  plot(lambda,cv,type=type,ylim=ylim,...)
  if (se)
    error.bars(lambda,cv+cv.error,cv-cv.error,width=1/length(lambda))
  invisible()
}

##' Make folds
##'
##' 
##' @title Make folds
##' @param n no descr
##' @param folds no descr
##' @return no descr
##' @author n.stadler
cv.fold <- function(n,folds=10){
  split(sample(1:n),rep(1:folds,length=n))
}

##' Graphical Lasso path with huge package
##'
##' 
##' @title Graphical Lasso path with huge package
##' @param s no descr
##' @param rholist no descr
##' @param penalize.diagonal no descr
##' @param trace no descr
##' @return no descr
##' @author n.stadler
hugepath <- function(s,rholist,penalize.diagonal=NULL,trace=NULL){
  #fit.huge <- huge(s,method = "glasso",cov.output =TRUE,verbose = FALSE)
  fit.huge <- huge(s,lambda=sort(rholist,decreasing=TRUE),method = "glasso",cov.output =TRUE,verbose = FALSE)
  wi <- sapply(fit.huge$icov,as.matrix,simplify='array')
  w <- sapply(fit.huge$cov,as.matrix,simplify='array')
  #return(list(wi=wi[,,length(fit.huge$lambda):1],w=w[,,length(fit.huge$lambda):1]))
  return(list(rholist=rholist,wi=wi[,,length(rholist):1],w=w[,,length(rholist):1]))
}

##' Additional thresholding
##'
##' 
##' @title Additional thresholding
##' @param n no descr
##' @param wi no descr
##' @param method no descr
##' @param trunc.k no descr
##' @return no descr
##' @author n.stadler
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

##' Crossvalidation for GLasso 
##'
##' 8! lambda-grid has to be increasing (see glassopath)
##' @title Crossvalidation for GLasso
##' @param x no descr
##' @param folds no descr
##' @param lambda lambda-grid (increasing!)
##' @param penalize.diagonal no descr
##' @param plot.it no descr
##' @param se no descr
##' @param include.mean no descr
##' @return no descr
##' @author n.stadler
cv.glasso <- function(x,folds=10,lambda,penalize.diagonal=FALSE,plot.it=FALSE,se=TRUE,include.mean=FALSE)
{
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
    fit.path <- glassopath(s,rholist=2*lambda/nrow(x),penalize.diagonal=penalize.diagonal,trace=0)
    if(length(omit)==1){
      residmat[cvfold,] <- -2*apply(fit.path$w,3,dmvnorm,log=TRUE,mean=mu,x=x[omit,,drop=FALSE])
    }else{
      residmat[cvfold,] <- apply(-2*apply(fit.path$w,3,dmvnorm,log=TRUE,mean=mu,x=x[omit,,drop=FALSE]),2,sum)
    }
  }
  cv <- apply(residmat,2,mean)
  cv.error <- sqrt(apply(residmat,2,var)/folds)
  gl.opt<-glasso(var(x),rho=2*lambda[which.min(cv)]/nrow(x),penalize.diagonal=penalize.diagonal)
  cat('la.opt:',lambda[which.min(cv)],'\n')
  w<-gl.opt$w
  wi<-gl.opt$wi
  wi[abs(wi)<10^{-3}]<-0
  colnames(w)<-rownames(w)<-colnames(wi)<-rownames(wi)<-colnames(x)  

  object <- list(rho.opt=2*lambda[which.min(cv)]/nrow(x),lambda=lambda,residmat=residmat,cv=cv,cv.error=cv.error,w=w,wi=wi,mu=colMeans(x))
  if (plot.it){
    plotCV(lambda,cv,cv.error,se=se)
  }
  invisible(object)
}

##' Cross-validated glasso with additional thresholding
##'
##'
##' Run glasso on a single dataset, using cross-validation to estimate the
##' penalty parameter lambda. Performs additional thresholding (optionally).
##' 
##' @title Cross-validated glasso with additional thresholding
##' @param x The input data. Needs to be a num.samples by dim.samples matrix.
##' @param include.mean Include mean in likelihood. TRUE / FALSE (default).
##' @param folds Number of folds in the cross-validation (default=10).
##' @param length.lambda Length of lambda path to consider (default=20).
##' @param lambdamin.ratio Ratio lambda.min/lambda.max.
##' @param penalize.diagonal If TRUE apply penalization to diagonal of inverse
##'        covariance as well. (default=FALSE)
##' @param trunc.method None / linear.growth (default) / sqrt.growth
##' @param trunc.k truncation constant, number of samples per predictor (default=5)
##' @param plot.it TRUE / FALSE (default)
##' @param se default=FALSE.
##' @param use.package 'glasso' or 'huge' (default).
##' @param verbose If TRUE, output la.min, la.max and la.opt (default=FALSE).
##' @return  Returns a list with named elements 'rho.opt', 'wi', 'wi.orig', 'mu', 
##'          Variable rho.opt is the optimal (scaled) penalization parameter (rho.opt=2*la.opt/n). 
##'          The variables wi and wi.orig are matrices of size dim.samples by dim.samples 
##'          containing the truncated and untruncated inverse covariance matrix. Variable 
##'          Mu mean of the input data.
##' @author n.stadler
##' @export
screen_cv.glasso <- function(x,include.mean=FALSE,
                             folds=10,length.lambda=20,lambdamin.ratio=ifelse(ncol(x)>nrow(x),0.01,0.001),penalize.diagonal=FALSE,
                             trunc.method='linear.growth',trunc.k=5,plot.it=FALSE,se=FALSE,use.package='huge',verbose=FALSE)
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

##' BIC.glasso
##'
##' 
##' @title BIC.glasso
##' @param x no descr
##' @param lambda no descr
##' @param penalize.diagonal no descr
##' @param plot.it no descr
##' @param use.package no descr
##' @param include.mean no descr
##' @return no descr
##' @author n.stadler
bic.glasso <- function(x,lambda,penalize.diagonal=FALSE,plot.it=TRUE,use.package='huge',include.mean=FALSE)
{
  ##glasso; lambda.opt with bic
  if(include.mean==TRUE){
    Mu <- colMeans(x)
  }
  if(include.mean==FALSE){
    Mu <-  rep(0,ncol(x))
  }
  samplecov <- var(x)

  if(is.null(lambda)){
    la <- lambda
  }else{
    la <- 2*lambda/nrow(x)
  }

  fit.path <- eval(as.name(paste(use.package,'path',sep='')))(samplecov,rholist=la,penalize.diagonal=penalize.diagonal,trace=0)

  loglik <- lapply(seq(length(fit.path$rholist)),
                   function(i){
                     return(sum(dmvnorm(x,mean=Mu,sigma=fit.path$w[,,i],log=TRUE)))
                   }
                   )
  degfree <- lapply(seq(length(fit.path$rholist)),
                    function(i){
                      wi <- fit.path$wi[,,i]
                      wi[abs(wi)<10^{-3}]<-0
                      p <- ncol(wi)
                      n.zero<-sum(wi==0)
                      if(include.mean==TRUE){
                        dfree <- p+(p*(p+1)/2)-n.zero/2
                      }
                      if(include.mean==FALSE){
                        dfree <- (p*(p+1)/2)-n.zero/2
                      }
                      return(dfree)
                    }
                    )
  loglik <- simplify2array(loglik,higher=TRUE)
  degfree <- simplify2array(degfree,higher=TRUE)
  myscore <- -loglik+log(nrow(x))*degfree/2

  if(is.null(lambda)){
    lambda <- 0.5*nrow(x)*fit.path$rholist
  }
  index.opt <- which.min(myscore)
  wi <- fit.path$wi[,,index.opt]
  wi[abs(wi)<10^{-3}]<-0
  w <- fit.path$w[,,index.opt]
  
  if (plot.it){
    plot(lambda,myscore,type='b',xlab='lambda')
  }
  
  list(rho.opt=2*lambda[which.min(myscore)]/nrow(x),lambda=lambda,la.opt=lambda[index.opt],bic.score=myscore,Mu=Mu,wi=wi,w=w)
}

##' AIC.glasso
##'
##' 
##' @title AIC.glasso
##' @param x no descr
##' @param lambda no descr
##' @param penalize.diagonal no descr
##' @param plot.it no descr
##' @param use.package no descr
##' @param include.mean no descr
##' @return no descr
##' @author n.stadler
aic.glasso <- function(x,lambda,penalize.diagonal=FALSE,plot.it=TRUE,use.package='huge',include.mean=FALSE)
{
  ##glasso; lambda.opt with aicc
  if(include.mean==TRUE){
    Mu <- colMeans(x)
  }
  if(include.mean==FALSE){
    Mu <-  rep(0,ncol(x))
  }
  
  samplecov <- var(x)

  if(is.null(lambda)){
    la <- lambda
  }else{
    la <- 2*lambda/nrow(x)
  }

  fit.path <- eval(as.name(paste(use.package,'path',sep='')))(samplecov,rholist=la,penalize.diagonal=penalize.diagonal,trace=0)

  loglik <- lapply(seq(length(fit.path$rholist)),
                   function(i){
                     return(sum(dmvnorm(x,mean=Mu,sigma=fit.path$w[,,i],log=TRUE)))
                   }
                   )
  degfree <- lapply(seq(length(fit.path$rholist)),
                    function(i){
                      wi <- fit.path$wi[,,i]
                      wi[abs(wi)<10^{-3}]<-0
                      p <- ncol(wi)
                      n.zero<-sum(wi==0)
                      if(include.mean==TRUE){
                        dfree <- p+(p*(p+1)/2)-n.zero/2
                      }
                      if(include.mean==FALSE){
                        dfree <- (p*(p+1)/2)-n.zero/2
                      }
                      return(dfree)
                    }
                    )
  loglik <- simplify2array(loglik,higher=TRUE)
  degfree <- simplify2array(degfree,higher=TRUE)
  myscore <- -loglik+2*degfree/2
  if(is.null(lambda)){
    lambda <- 0.5*nrow(x)*fit.path$rholist
  }
  index.opt <- which.min(myscore)
  wi <- fit.path$wi[,,index.opt]
  wi[abs(wi)<10^{-3}]<-0
  w <- fit.path$w[,,index.opt]

  if (plot.it){
    plot(lambda,myscore,type='b',xlab='lambda')
  }
  
  list(rho.opt=2*lambda[which.min(myscore)]/nrow(x),lambda=lambda,la.opt=lambda[index.opt],bic.score=myscore,Mu=Mu,wi=wi,w=w)
}

##' BIC-tuned glasso with additional thresholding
##'
##' 
##' @title BIC-tuned glasso with additional thresholding
##' @param x  The input data. Needs to be a num.samples by dim.samples matrix.
##' @param include.mean Include mean in likelihood. TRUE / FALSE (default).
##' @param length.lambda Length of lambda path to consider (default=20).
##' @param lambdamin.ratio Ratio lambda.min/lambda.max.
##' @param penalize.diagonal  If TRUE apply penalization to diagonal of inverse
##'        covariance as well. (default=FALSE)
##' @param plot.it TRUE / FALSE (default)
##' @param trunc.method None / linear.growth (default) / sqrt.growth
##' @param trunc.k truncation constant, number of samples per predictor (default=5)
##' @param use.package 'glasso' or 'huge' (default).
##' @param verbose If TRUE, output la.min, la.max and la.opt (default=FALSE).
##' @return  Returns a list with named elements 'rho.opt', 'wi', 'wi.orig', 
##'          Variable rho.opt is the optimal (scaled) penalization parameter (rho.opt=2*la.opt/n). 
##'          The variables wi and wi.orig are matrices of size dim.samples by dim.samples 
##'          containing the truncated and untruncated inverse covariance matrix. 
##' @author n.stadler
##' @export
screen_bic.glasso <- function(x,include.mean=TRUE,
                              length.lambda=20,lambdamin.ratio=ifelse(ncol(x)>nrow(x),0.01,0.001),penalize.diagonal=FALSE,
                              plot.it=FALSE,trunc.method='linear.growth',trunc.k=5,use.package='huge',verbose=FALSE){

  gridmax <- lambda.max(x)
  gridmin <- gridmax*lambdamin.ratio
  my.grid <- make_grid(gridmin,gridmax,length.lambda)[length.lambda:1]
  
  fit.bicgl <- bic.glasso(x,lambda=my.grid,penalize.diagonal=penalize.diagonal,plot.it=plot.it,
                          use.package=use.package,include.mean=include.mean)
  if(verbose){
    cat('la.min:',gridmin,'\n')
    cat('la.max:',gridmax,'\n')
    cat('la.opt:',fit.bicgl$la.opt,'\n')
  }
  wi <- fit.bicgl$wi
  wi.trunc <- mytrunc.method(n=nrow(x),wi=wi,method=trunc.method,trunc.k=trunc.k)$wi
  list(rho.opt=fit.bicgl$rho.opt,wi=wi.trunc,wi.orig=wi) 
}

##' AIC-tuned glasso with additional thresholding
##'
##' 
##' @title AIC-tuned glasso with additional thresholding
##' @param x  The input data. Needs to be a num.samples by dim.samples matrix.
##' @param include.mean Include mean in likelihood. TRUE / FALSE (default).
##' @param length.lambda Length of lambda path to consider (default=20).
##' @param lambdamin.ratio Ratio lambda.min/lambda.max.
##' @param penalize.diagonal  If TRUE apply penalization to diagonal of inverse
##'        covariance as well. (default=FALSE)
##' @param plot.it TRUE / FALSE (default)
##' @param trunc.method None / linear.growth (default) / sqrt.growth
##' @param trunc.k truncation constant, number of samples per predictor (default=5)
##' @param use.package 'glasso' or 'huge' (default).
##' @param verbose If TRUE, output la.min, la.max and la.opt (default=FALSE).
##' @return  Returns a list with named elements 'rho.opt', 'wi', 'wi.orig'.
##'          Variable rho.opt is the optimal (scaled) penalization parameter (rho.opt=2*la.opt/n). 
##'          The variables wi and wi.orig are matrices of size dim.samples by dim.samples 
##'          containing the truncated and untruncated inverse covariance matrix. 
##' @author n.stadler
##' @export
screen_aic.glasso <- function(x,include.mean=TRUE,length.lambda=20,lambdamin.ratio=ifelse(ncol(x)>nrow(x),0.01,0.001),
                              penalize.diagonal=FALSE,plot.it=FALSE,
                              trunc.method='linear.growth',trunc.k=5,use.package='huge',verbose=FALSE){

  gridmax <- lambda.max(x)
  gridmin <- gridmax*lambdamin.ratio
  my.grid <- make_grid(gridmin,gridmax,length.lambda)[length.lambda:1]
  
  fit.aicgl <- aic.glasso(x,lambda=my.grid,penalize.diagonal=penalize.diagonal,plot.it=plot.it,
                          use.package=use.package,include.mean=include.mean)
  if(verbose){
    cat('la.min:',gridmin,'\n')
    cat('la.max:',gridmax,'\n')
    cat('la.opt:',fit.aicgl$la.opt,'\n')
  }
  wi <- fit.aicgl$wi
  wi.trunc <- mytrunc.method(n=nrow(x),wi=wi,method=trunc.method,trunc.k=trunc.k)$wi
  list(rho.opt=fit.aicgl$rho.opt,wi=wi.trunc,wi.orig=wi) 
}

##' Screen_lasso
##'
##' 
##' @title Screen_lasso
##' @param x no descr
##' @param include.mean no descr
##' @param trunc.method no descr
##' @param trunc.k no descr
##' @return no descr
##' @author n.stadler
screen_lasso <- function(x,include.mean=NULL,
                         trunc.method='linear.growth',trunc.k=5){
  
  wi <- adalasso.net(x, k = 10,use.Gram=FALSE,both=FALSE,verbose=FALSE)$pcor.lasso
  wi.trunc <- mytrunc.method(n=nrow(x),wi=wi,method=trunc.method,trunc.k=trunc.k)$wi
  list(rho.opt=NULL,wi=wi.trunc,wi.orig=wi)
}

##' Shrinkage approach for estimating Gaussian graphical model
##'
##' 
##' @title Shrinkage approach for estimating Gaussian graphical model
##' @param x  The input data. Needs to be a num.samples by dim.samples matrix.
##' @param include.mean Include mean in likelihood. TRUE / FALSE (default).
##' @param trunc.method None / linear.growth (default) / sqrt.growth
##' @param trunc.k truncation constant, number of samples per predictor (default=5)
##' @return Returns a list with named elements 'rho.opt', 'wi', 'wi.orig'.
##'          Variable rho.opt=NULL (no tuning parameter involved). 
##'          The variables wi and wi.orig are matrices of size dim.samples by dim.samples 
##'          containing the truncated and untruncated inverse covariance matrix.
##' @author n.stadler
##' @export
screen_shrink <- function(x,include.mean=NULL,
                          trunc.method='linear.growth',trunc.k=5){
  wi <- ggm.estimate.pcor(x)
  adj <- performance.pcor(wi, fdr=TRUE,verbose=FALSE,plot.it=FALSE)$adj
  wi[adj==0] <- 0
  wi.trunc <- mytrunc.method(n=nrow(x),wi=wi,method=trunc.method,trunc.k=trunc.k)$wi
  list(rho.opt=NULL,wi=wi.trunc,wi.orig=wi)
}

##' Node-wise Lasso-regressions for GGM estimation
##' 
##' (Meinshausen-Buehlmann approach)
##' 
##' @title Node-wise Lasso-regressions for GGM estimation
##' @param x The input data. Needs to be a num.samples by dim.samples matrix.
##' @param include.mean Include mean in likelihood. TRUE / FALSE (default).
##' @param folds Number of folds in the cross-validation (default=10).
##' @param length.lambda  Length of lambda path to consider (default=20).
##' @param lambdamin.ratio  Ratio lambda.min/lambda.max.
##' @param penalize.diagonal If TRUE apply penalization to diagonal of inverse
##'        covariance as well. (default=FALSE)
##' @param trunc.method None / linear.growth (default) / sqrt.growth
##' @param trunc.k truncation constant, number of samples per predictor (default=5)
##' @param plot.it TRUE / FALSE (default)
##' @param se default=FALSE.
##' @param verbose If TRUE, output la.min, la.max and la.opt (default=FALSE).
##' @return  Returns a list with named elements 'rho.opt', 'wi'.
##'          Variable rho.opt is the optimal (scaled) penalization parameter (rho.opt=2*la.opt/n). 
##'          The variables wi is a matrix of size dim.samples by dim.samples 
##'          containing the truncated inverse covariance matrix. Variable 
##'          Mu mean of the input data.
##' @author n.stadler
##' @export
screen_mb <- function(x,include.mean=NULL,
                      folds=10,length.lambda=20,lambdamin.ratio=ifelse(ncol(x)>nrow(x),0.01,0.001),
                      penalize.diagonal=FALSE,trunc.method='linear.growth',trunc.k=5,
                      plot.it=FALSE,se=FALSE,verbose=FALSE)
{
  p <- ncol(x)
  gridmax <- lambda.max(x)
  gridmin <- gridmax*lambdamin.ratio
  lambda <- make_grid(gridmin,gridmax,length.lambda)[length.lambda:1]
  
  colnames(x)<-paste('x',1:ncol(x),sep='')  
  all.folds <- cv.fold(nrow(x),folds)
  residmat <- matrix(NA,folds,length.lambda)
  
  for (cvfold in 1:folds){
    omit <- all.folds[[cvfold]]
    s <- var(x[-omit,])
    fit.path <- glassopath(s,rholist=2*lambda/nrow(x),penalize.diagonal=penalize.diagonal,trace=0,approx=TRUE)
    myres <- sapply(1:p,function(j){colMeans((x[omit,j]-x[omit,-j,drop=FALSE]%*%fit.path$wi[-j,j,])^2)})
    residmat[cvfold,] <- rowSums(myres)
  }
  cv <- apply(residmat,2,mean)
  cv.error <- sqrt(apply(residmat,2,var)/folds)
  gl.opt<-glasso(var(x),rho=2*lambda[which.min(cv)]/nrow(x),penalize.diagonal=penalize.diagonal,approx=TRUE)
  if(verbose){
    cat('la.min:',gridmin,'\n')
    cat('la.max:',gridmax,'\n')
    cat('la.opt:',lambda[which.min(cv)],'\n')
  }
  
  wi<-Beta2parcor(gl.opt$wi)
  #wi[abs(wi)<10^{-3}]<-0
  colnames(wi)<-rownames(wi)<-colnames(x)
  wi <- mytrunc.method(n=nrow(x),wi=wi,method=trunc.method,trunc.k=trunc.k)$wi
  
  if (plot.it){
    plotCV(lambda,cv,cv.error,se=se)
  }
  
  list(rho.opt=2*lambda[which.min(cv)]/nrow(x),wi=wi)
}

##' Screen_mb2
##'
##' 
##' @title Screen_mb2
##' @param x no descr
##' @param include.mean no descr
##' @param length.lambda no descr
##' @param trunc.method no descr
##' @param trunc.k no descr
##' @param plot.it no descr
##' @param verbose no descr
##' @return no descr
##' @author n.stadler
screen_mb2 <- function(x,include.mean=NULL,length.lambda=20,
                       trunc.method='linear.growth',trunc.k=5,plot.it=FALSE,verbose=FALSE)
{
  p <- ncol(x)
  beta <- rep(0,p)
  Beta <- sapply(1:p,function(j){
    fit.cv <- cv.glmnet(x[,-j],x[,j],standardize=FALSE,nlambda=length.lambda)
    if(verbose){
      cat('vertex no ',j,' laopt ',0.5*nrow(x)*fit.cv$lambda.min,'\n')
    }
    if(plot.it){
      if(j==1){
        plot(fit.cv$lambda,fit.cv$cvm/max(fit.cv$cvm),type='b',cex=0.5,ylim=c(0,1))
      }else{
        lines(fit.cv$lambda,fit.cv$cvm/max(fit.cv$cvm),type='b',cex=0.5)
      }
    }
    beta[-j] <- as.numeric(coef(fit.cv,s='lambda.min')[-1])
    return(beta)})
                 
  wi<-Beta2parcor(Beta)
  #wi[abs(wi)<10^{-3}]<-0
  colnames(wi)<-rownames(wi)<-colnames(x)
  wi <- mytrunc.method(n=nrow(x),wi=wi,method=trunc.method,trunc.k=trunc.k)$wi
  
  list(rho.opt=NULL,wi=wi)
}

##' Screen_full
##'
##' 
##' @title Screen_full
##' @param x no descr
##' @param include.mean no descr
##' @param length.lambda no descr
##' @param trunc.method no descr
##' @param trunc.k no descr
##' @return no descr
##' @author n.stadler
screen_full <- function(x,include.mean=NULL,length.lambda=NULL,trunc.method=NULL,trunc.k=NULL){
 wi <- diag(1,ncol(x))
 list(rho.opt=NULL,wi=wi)
}


###################################################################
##-------covMethod------
###################################################################

##' Compute covariance matrix
##'
##' 
##' @title Compute covariance matrix
##' @param x no descr
##' @param include.mean no descr
##' @param covMethod no descr
##' @return no descr
##' @author n.stadler
mcov <- function(x,include.mean,covMethod='ML'){
  if(covMethod=='ML'){
    return(cov.wt(x,method='ML',center=include.mean)$cov)
  }
  if(covMethod=='var'){
    return(var(x))
  }
}
  
#################################################################
##--------P-VALUES----------
#################################################################

##' MLE in GGM
##'
##' 
##' @title MLE in GGM
##' @param x no descr
##' @param wi no descr
##' @param algorithm no descr
##' @param rho no descr
##' @param include.mean no descr
##' @return no descr
##' @author n.stadler
mle.ggm <- function(x,wi,algorithm='glasso_rho0',rho=NULL,include.mean){
  if(is.null(rho)){
    algorithm <- 'glasso_rho0'
    cat('screening method with rho.opt=NULL; set algorithm="glasso_rho0" in mle.ggm')
  }
  if (algorithm=='glasso'){
    fit.mle <- glasso(mcov(x=x,include.mean=include.mean),rho=rho,zero=which(wi==0,arr.ind=TRUE))
    return(list(w=fit.mle$w,wi=fit.mle$wi))
  }
  if (algorithm=='glasso_rho0'){
    fit.mle <- glasso(mcov(x=x,include.mean=include.mean),rho=10^{-10},zero=which(wi==0,arr.ind=TRUE))
    return(list(w=fit.mle$w,wi=fit.mle$wi))
  }
  if (algorithm=='fitcongraph'){
    s <- mcov(x=x,include.mean=include.mean)
    adj <- wi!=0
    colnames(s) <- rownames(s) <- colnames(adj) <- rownames(adj) <- 1:ncol(x)
    fit.mle <- fitConGraph(adj,s,nrow(x))
    w <- fit.mle$Shat
    wi <- solve(w)
    wi[!adj] <- 0
    return(list(w=fit.mle$Shat,wi=wi))
  }
}

##' Log-likelihood-ratio statistics used in Differential Network
##'
##' 
##' @title Log-likelihood-ratio statistics used in DiffNet
##' @param x1 data-matrix sample 1
##' @param x2 data-matrix sample 2
##' @param x pooled data-matrix
##' @param sig1 covariance sample 1
##' @param sig2 covariance sample 2
##' @param sig pooled covariance
##' @param mu1 mean sample 1
##' @param mu2 mean sample 2
##' @param mu pooled mean
##' @return Returns a list with named elements 'twiceLR', 'sig1', 'sig2', 'sig'.
##'         'twiceLR' is twice the log-likelihood-ratio statistic.
##' @author n.stadler
##' @export
logratio <- function(x1,x2,x,sig1,sig2,sig,mu1,mu2,mu){
  twiceLR <- 2*(sum(dmvnorm(x1,mean=mu1,sigma=sig1,log=TRUE))+sum(dmvnorm(x2,mean=mu2,sigma=sig2,log=TRUE))-sum(dmvnorm(x,mean=mu,sigma=sig,log=TRUE)))
  list(twiceLR=twiceLR,sig1=sig1,sig2=sig2,sig=sig)
}

##' Compute Information Matrix of Gaussian Graphical Model
##'
##' computes E_0[s(Y;Omega)s(Y;Omega)'] where s(Y;Omega)=(d/dOmega) LogLik
##' @title Information Matrix of Gaussian Graphical Model
##' @param Sig Sig=solve(SigInv) true covariance matrix under H0
##' @param include.mean no descr
##' @return no descr
##' @author n.stadler
inf.mat<-function(Sig,include.mean=FALSE){
  k <- ncol(Sig)
  uptri.rownr <- row(Sig)[upper.tri(Sig,diag=TRUE)]
  uptri.colnr <- col(Sig)[upper.tri(Sig,diag=TRUE)]

  if (include.mean==TRUE){
    p <- k+k*(k+1)/2
    infmat <- matrix(0,p,p)
    for (i in seq(p-k)){
      colnr.i <- uptri.colnr[i]
      rownr.i <- uptri.rownr[i]
      for (j in seq(p-k)){
        colnr.j <- uptri.colnr[j]
        rownr.j <- uptri.rownr[j]
        infmat[i,j] <- Sig[colnr.i,rownr.j]*Sig[rownr.i,colnr.j]+Sig[colnr.i,colnr.j]*Sig[rownr.i,rownr.j]
      }
    }
    infmat[(p-k+1):p,(p-k+1):p] <- solve(Sig)
  }else{
    p <- k*(k+1)/2
    infmat <- matrix(NA,p,p)
    for (i in seq(p)){
      colnr.i <- uptri.colnr[i]
      rownr.i <- uptri.rownr[i]
      for (j in seq(p)){
        colnr.j <- uptri.colnr[j]
        rownr.j <- uptri.rownr[j]
        infmat[i,j] <- Sig[colnr.i,rownr.j]*Sig[rownr.i,colnr.j]+Sig[colnr.i,colnr.j]*Sig[rownr.i,rownr.j]
      }
    }
  }
  return(infmat)
}

##' Calculates weight-matrix and eigenvalues
##'
##' calculation based on true information matrix
##' @title Weight-matrix and eigenvalues
##' @param imat no descr
##' @param act I_uv
##' @param act1 I_u
##' @param act2 I_v
##' @return no descr
##' @author n.stadler
ww.mat <- function(imat,act,act1,act2){
 
  bfg <- rbind(imat[act1,act,drop=FALSE],imat[act2,act,drop=FALSE])
  bf <- matrix(0,length(act1)+length(act2),length(act1)+length(act2))
  bf[1:length(act1),1:length(act1)] <- imat[act1,act1]
  bf[length(act1)+(1:(length(act2))),length(act1)+(1:(length(act2)))] <- imat[act2,act2]
  bg <- 2*imat[act,act,drop=FALSE]
  mat <- rbind(cbind(diag(1,length(act1)+length(act2)),bfg%*%solve(bg)),cbind(-t(bfg)%*%solve(bf),diag(-1,length(act))))
  eval <- as.double(eigen(mat)$values)
  eval[abs(eval)<10^{-6}] <- 0
  
  return(list(ww.mat=mat,eval=eval))
}

##' Calculates eigenvalues of weight-matrix (using 1st order simplification)
##'
##' calculation based on true information matrix
##' @title Calculates eigenvalues of weight-matrix (using 1st order simplification)
##' @param imat no descr
##' @param act I_uv
##' @param act1 I_u
##' @param act2 I_v
##' @return no descr
##' @author n.stadler
ww.mat2 <- function(imat,act,act1,act2){

  dimf <- length(act1)+length(act2)
  dimg <- length(act)

  bfg <- rbind(imat[act1,act,drop=FALSE],imat[act2,act,drop=FALSE])
  bgf <- t(bfg)
  bf <- matrix(0,length(act1)+length(act2),length(act1)+length(act2))
  bf[1:length(act1),1:length(act1)] <- imat[act1,act1]
  bf[length(act1)+(1:(length(act2))),length(act1)+(1:(length(act2)))] <- imat[act2,act2]
  bg <- 2*imat[act,act]

  if (dimf>=dimg){
    mat <- bgf%*%solve(bf)%*%bfg%*%solve(bg)
    eval<-rep(1,dimf-dimg)
  }
  if (dimf<dimg){
    mat <- bfg%*%solve(bg)%*%bgf%*%solve(bf)
    eval<-rep(-1,dimg-dimf)
  }
  eval.mu <- as.double(eigen(mat)$values)
  eval2 <- 1-eval.mu
  eval2[abs(eval2)<10^{-6}] <- 0
  eval <- c(eval,sqrt(eval2),-sqrt(eval2))
  return(list(ww.mat=mat,eval=eval))
}

##' Compute beta-matrix 
##'
##' beta-matrix=E[s_ind1(Y;sig1) s_ind2(Y;sig2)'|sig]
##' @title Compute beta-matrix 
##' @param ind1 no descr
##' @param ind2 no descr
##' @param sig1 no descr
##' @param sig2 no descr
##' @param sig no descr
##' @return no descr
##' @author n.stadler
beta.mat<-function(ind1,ind2,sig1,sig2,sig){
  uptri.rownr <- row(sig)[upper.tri(sig,diag=TRUE)]
  uptri.colnr <- col(sig)[upper.tri(sig,diag=TRUE)]
  betamat <- .C('betamat_diffnet',betamat=matrix(0,length(ind1),length(ind2)),
                ind1=as.integer(ind1-1),ind2=as.integer(ind2-1),##8!indexing
                uptrirownr=uptri.rownr,uptricolnr=uptri.colnr,
                lind1=length(ind1),lind2=length(ind2),
                sig1=sig1,sig2=sig2,sig=sig,k=ncol(sig))$betamat
  
  return(betamat)
}

## beta.mat<-function(ind1,ind2,sig1,sig2,sig){
##   k <- ncol(sig)
##   p <- k*(k+1)/2
##   uptri.rownr <- row(sig)[upper.tri(sig,diag=TRUE)]
##   uptri.colnr <- col(sig)[upper.tri(sig,diag=TRUE)]

##   bmat <- matrix(NA,p,p)
##   for (i in ind1){
##     colnr.i <- uptri.colnr[i]
##     rownr.i <- uptri.rownr[i]
##     for (j in ind2){
##       colnr.j <- uptri.colnr[j]
##       rownr.j <- uptri.rownr[j]
##       bmat[i,j] <- sig[colnr.i,rownr.i]*sig[colnr.j,rownr.j]-sig1[colnr.i,rownr.i]*sig[colnr.j,rownr.j]-sig[colnr.i,rownr.i]*sig2[colnr.j,rownr.j]+sig1[colnr.i,rownr.i]*sig2[colnr.j,rownr.j]+sig[colnr.i,rownr.j]*sig[rownr.i,colnr.j]+sig[colnr.i,colnr.j]*sig[rownr.i,rownr.j]
##     }
##   }
##   return(bmat[ind1,ind2,drop=FALSE])
## }

##' Compute Q-matrix 
##'
##' 
##' @title Compute Q-matrix 
##' @param sig no descr
##' @param sig.a no descr
##' @param sig.b no descr
##' @param act.a no descr
##' @param act.b no descr
##' @param ss no descr
##' @return no descr
##' @author n.stadler
q.matrix3 <- function(sig,sig.a,sig.b,act.a,act.b,ss){
  b.ab<-beta.mat(act.a,act.b,sig.a,sig.b,sig)
  aa<-seq(1,length(act.a))[!(act.a%in%ss)]
  bb<-seq(1,length(act.b))[!(act.b%in%ss)]
  s.a<-seq(1,length(act.a))[(act.a%in%ss)]
  s.b<-seq(1,length(act.b))[(act.b%in%ss)]
  
  if(length(ss)==0){
    return(b.ab[aa,bb,drop=FALSE])
  }else{
    return(b.ab[aa,bb,drop=FALSE]-(b.ab[aa,s.b,drop=FALSE]%*%solve(b.ab[s.a,s.b,drop=FALSE])%*%b.ab[s.a,bb,drop=FALSE]))
  }
}

##' q.matrix4
##'
##' 
##' @title q.matrix4
##' @param b.mat no descr
##' @param act.a no descr
##' @param act.b no descr
##' @param ss no descr
##' @return no descr
##' @author n.stadler
q.matrix4 <- function(b.mat,act.a,act.b,ss){
 ##Estimate Q
    aa<-seq(1,length(act.a))[!(act.a%in%ss)]
    bb<-seq(1,length(act.b))[!(act.b%in%ss)]
    s.a<-seq(1,length(act.a))[(act.a%in%ss)]
    s.b<-seq(1,length(act.b))[(act.b%in%ss)]

    if(length(ss)==0){
        return(b.mat[aa,bb,drop=FALSE])
    }else{
        return(b.mat[aa,bb,drop=FALSE]-(b.mat[aa,s.b,drop=FALSE]%*%solve(b.mat[s.a,s.b,drop=FALSE])%*%b.mat[s.a,bb,drop=FALSE]))
    }
}

##' Compute weights of sum-w-chi2 (1st order simplification)
##'
##' 
##' @title Weights of sum-w-chi2
##' @param sig1 no descr
##' @param sig2 no descr
##' @param sig no descr
##' @param act1 no descr
##' @param act2 no descr
##' @param act no descr
##' @param include.mean no descr
##' @return no descr
##' @author n.stadler
est2.ww.mat2 <- function(sig1,sig2,sig,act1,act2,act,include.mean=FALSE){

  k <- nrow(sig1)
  
  #######################
  ##dimension of models##
  #######################
  dimf1 <- length(act1)
  dimf2 <- length(act2)
  dimf <- dimf1+dimf2
  dimg <- length(act)
  if(include.mean==TRUE){
    dimf <- dimf+2*k
    dimg <- dimg+k
  }

  ###############
  ##Beta-matrix##
  ###############
  b1.act<- beta.mat(act,act,sig,sig,sig1)
  b2.act<- beta.mat(act,act,sig,sig,sig2)
  b1.act1<- beta.mat(act1,act1,sig1,sig1,sig1)
  b2.act2<- beta.mat(act2,act2,sig2,sig2,sig2)
  b1.act1.act<- beta.mat(act1,act,sig1,sig,sig1)
  b2.act2.act<- beta.mat(act2,act,sig2,sig,sig2)
  
  #if (any(c(dimf1,dimf2,dimg)>min(nrow(xx1),nrow(xx2)))){
  #  warning('|active-set| > n')
  #  return(list(eval=NA))
  #}else{
  bfg <- rbind(b1.act1.act,b2.act2.act)
  bgf <- t(bfg)
  bf <- matrix(0,dimf1+dimf2,dimf1+dimf2)
  bf[1:dimf1,1:dimf1] <- b1.act1
  bf[dimf1+(1:dimf2),dimf1+(1:dimf2)] <- b2.act2
  bg <- b1.act+b2.act
  if (dimf>=dimg){
    mat <- bgf%*%solve(bf)%*%bfg%*%solve(bg)
    eval<-rep(1,dimf-dimg)
  }
  if (dimf<dimg){
    mat <- bfg%*%solve(bg)%*%bgf%*%solve(bf)
    eval<-rep(-1,dimg-dimf)
  }
  eval.mu.complex<-eigen(mat)$values
  onemineval.mu <- 1-as.double(eval.mu.complex)
  onemineval.mu[abs(onemineval.mu)<10^{-6}]<-0
  eval <- c(eval,sqrt(onemineval.mu),-sqrt(onemineval.mu))
  if(include.mean==TRUE){
    eval <- c(rep(0,2*k),eval)
  }
  return(list(ww.mat=mat,eval=eval,eval.mu.complex=eval.mu.complex))
}

##' Compute weights of sum-w-chi2 (2nd order simplification)
##' 
##' *expansion of W in two directions ("dimf>dimg direction" & "dimf>dimg direction") 
##' *simplified computation of weights is obtained by assuming H0 and that X_u~X_v holds
##' @title Weights of sum-w-chi2
##' @param sig1 no descr
##' @param sig2 no descr
##' @param sig no descr
##' @param act1 no descr
##' @param act2 no descr
##' @param act no descr
##' @param include.mean no descr
##' @return no descr
##' @author n.stadler
est2.my.ev2 <- function(sig1,sig2,sig,act1,act2,act,include.mean=FALSE){

  k <- nrow(sig1)

  #######################
  ##dimension of models##
  #######################
  dimf1 <- length(act1)
  dimf2 <- length(act2)
  dimf <- dimf1+dimf2
  dimg <- length(act)
  if(include.mean==TRUE){
    dimf <- dimf+2*k
    dimg <- dimg+k
  }
  ##########################
  ##intersection of models##
  ##########################
  ss <- intersect(act,intersect(act1,act2))
  
  if(length(ss)==0){warning('no intersection between models')}
  
  aa <- setdiff(act1,ss)
  bb <- setdiff(act2,ss)
  cc <- setdiff(act,ss)
 
  ev.aux <- ev.aux.complex <- numeric(0)
  if (dimf>=dimg){
    if (length(cc)!=0){
      qcc <- q.matrix3(sig1,sig,sig,act,act,ss)+q.matrix3(sig2,sig,sig,act,act,ss)
      aux.mat <- matrix(0,length(cc),length(cc))
      if(length(aa)!=0){
        qac <- q.matrix3(sig1,sig1,sig,act1,act,ss)
        qaa <- q.matrix3(sig1,sig1,sig1,act1,act1,ss)
        aux.mat <- aux.mat+(t(qac)%*%solve(qaa)%*%qac)%*%solve(qcc)
      }
      if(length(bb)!=0){
        qbc <- q.matrix3(sig2,sig2,sig,act2,act,ss)
        qbb <- q.matrix3(sig2,sig2,sig2,act2,act2,ss)
        aux.mat <- aux.mat+(t(qbc)%*%solve(qbb)%*%qbc)%*%solve(qcc)
      }
      ev.aux.complex <- eigen(aux.mat)$values
      ev.aux <- as.double(ev.aux.complex)
      ev.aux <-  sqrt(1-ev.aux)
    }
    eval<-rep(1,dimf-dimg)
    eval <- c(eval,rep(0,2*length(ss)),ev.aux,-ev.aux)
  }## end if (dimf>=dimg){
  if (dimf<dimg){
    if (length(cc)!=0){
      qcc <- q.matrix3(sig1,sig,sig,act,act,ss)+q.matrix3(sig2,sig,sig,act,act,ss)
      if(length(aa)!=0){
        qac <- q.matrix3(sig1,sig1,sig,act1,act,ss)
        qaa <- q.matrix3(sig1,sig1,sig1,act1,act1,ss)
        oaa <- qac%*%solve(qcc)%*%t(qac)
        aux.mat.aa<- oaa%*%solve(qaa)
        aux.mat <- aux.mat.aa
      }
      if(length(bb)!=0){
        qbc <- q.matrix3(sig2,sig2,sig,act2,act,ss)
        qbb <- q.matrix3(sig2,sig2,sig2,act2,act2,ss)
        obb <- qbc%*%solve(qcc)%*%t(qbc)
        aux.mat.bb<- obb%*%solve(qbb)
        aux.mat <- aux.mat.bb
      }
      if((length(aa)!=0)&(length(bb)!=0)){
        oab <- qac%*%solve(qcc)%*%t(qbc)
        aux.mat<- rbind(cbind(aux.mat.aa,oab%*%solve(qbb)),cbind(t(oab)%*%solve(qaa),aux.mat.bb))
      }
    }
    if ((length(aa)!=0)|(length(bb)!=0)){ ##if (length(aa)==0)&(length(bb)==0) 'MI included in MJ; therefore -X^2(dim(MJ)-dim(MI)) distributed'
      ev.aux.complex <- eigen(aux.mat)$values
      ev.aux <- as.double(ev.aux.complex)
      ev.aux <-  sqrt(1-ev.aux)
    }
    eval<-rep(-1,dimg-dimf)
    eval <- c(eval,rep(-1,length(ss)),rep(1,length(ss)),rep(0,2*length(ss)),ev.aux,-ev.aux)
  }## end if (dimf<dimg){
  
  if(include.mean==TRUE){
    eval <- c(rep(0,2*k),eval)
  }
  return(list(eval=eval,ev.aux.complex=ev.aux.complex))
}

##' Compute weights of sum-of-weighted-chi2s
##'
##' *'2nd order simplification':
##'   1) Factor out (1-vi)^(d1+d2) "expansion in dimf>dimg direction (old terminology)"
##'   2) Factor out (1-mu)^d0 
##' *simplified computation of weights is obtained without further invoking H0, or assuming X_u~X_v
##' 
##' @title Compute weights of sum-of-weighted-chi2s
##' @param sig1 MLE (covariance matrix) sample 1
##' @param sig2 MLE (covariance matrix) sample 2
##' @param sig Pooled MLE (covariance matrix)
##' @param act1 Active-set sample 1
##' @param act2 Active-set sample 2
##' @param act Pooled active-set
##' @param include.mean Should the mean be in cluded in the likelihood?
##' @return Eigenvalues of M, respectively the weights.
##' @author n.stadler
est2.my.ev3 <- function(sig1,sig2,sig,act1,act2,act,include.mean=FALSE){

  show.warn <- FALSE
  k <- nrow(sig1)

  #######################
  ##dimension of models##
  #######################
  dimf1 <- length(act1)
  dimf2 <- length(act2)
  dimf <- dimf1+dimf2
  dimg <- length(act)
  if(include.mean==TRUE){
    dimf <- dimf+2*k
    dimg <- dimg+k
  }
  ##########################
  ##intersection of models##
  ##########################
  ss <- intersect(act,intersect(act1,act2))
  
  if((length(ss)==0)&show.warn){warning('no intersection between models')}
  
  aa <- setdiff(act1,ss)
  bb <- setdiff(act2,ss)
  cc <- setdiff(act,ss)
 
  ev.aux <- ev.aux.complex <- numeric(0)
  no.zero.ev.aux <- 0
  #if (dimf>=dimg){
    if (length(cc)!=0){
      qcc1 <- q.matrix3(sig1,sig,sig,act,act,ss)
      qcc2 <- q.matrix3(sig2,sig,sig,act,act,ss)
      bmat <- beta.mat(act,act,sig,sig,sig1)+beta.mat(act,act,sig,sig,sig2)
      qcc12 <- q.matrix4(bmat,act,act,ss)
      aux.mat <- diag(1,length(cc))-qcc1%*%solve(qcc12)-qcc2%*%solve(qcc12)
      if(length(aa)!=0){
        qac <- q.matrix3(sig1,sig1,sig,act1,act,ss)
        qaa <- q.matrix3(sig1,sig1,sig1,act1,act1,ss)
        aux.mat <- aux.mat+(t(qac)%*%solve(qaa)%*%qac)%*%solve(qcc12)
      }
      if(length(bb)!=0){
        qbc <- q.matrix3(sig2,sig2,sig,act2,act,ss)
        qbb <- q.matrix3(sig2,sig2,sig2,act2,act2,ss)
        aux.mat <- aux.mat+(t(qbc)%*%solve(qbb)%*%qbc)%*%solve(qcc12)
      }
      ev.aux.complex <- eigen(aux.mat)$values
      ev.aux <- as.double(ev.aux.complex)
      zero.ev.aux <- abs(ev.aux)<10^{-10}
      no.zero.ev.aux <- sum(zero.ev.aux)
      ev.aux <- ev.aux[!zero.ev.aux]
      ev.aux <-  sqrt(1-ev.aux)
    }
  eval<-c(rep(1,dimf-dimg+no.zero.ev.aux),rep(-1,no.zero.ev.aux))
  eval <- c(eval,rep(0,2*length(ss)),ev.aux,-ev.aux)
  #}## end if (dimf>=dimg){
  
  if(include.mean==TRUE){
    eval <- c(rep(0,2*k),eval)
  }
  return(list(eval=eval,ev.aux.complex=ev.aux.complex))
}

##' P-value aggregation
##'
##' 
##' @title P-value aggregation (Meinshausen et al 2009)
##' @param gamma see Meinshausen et al 2009
##' @param pval vector of p-values
##' @return inf-quantile aggregated p-value
##' @author n.stadler
agg.pval <- function(gamma,pval){
    min(quantile(pval/gamma,probs=gamma),1)
}

##' P-value calculation 
##'
##' 
##' @title P-value calculation
##' @param x1 no descr
##' @param x2 no descr
##' @param x no descr
##' @param sig1 no descr
##' @param sig2 no descr
##' @param sig no descr
##' @param mu1 no descr
##' @param mu2 no descr
##' @param mu no descr
##' @param act1 no descr
##' @param act2 no descr
##' @param act no descr
##' @param compute.evals no descr
##' @param include.mean no descr
##' @param method.compquadform no descr
##' @param acc no descr
##' @param epsabs no descr
##' @param epsrel no descr
##' @param show.warn no descr
##' @return no descr
##' @author n.stadler
diffnet_pval <- function(x1,x2,x,sig1,sig2,sig,mu1,mu2,mu,act1,act2,act,compute.evals,include.mean,method.compquadform,acc,epsabs,epsrel,show.warn){
  ##########################
  ##compute test-statistic##
  ##########################
  teststat <- logratio(x1,x2,x,sig1,sig2,sig,mu1,mu2,mu)$twiceLR
  #################################
  ##compute weights of sum-w-chi2##
  #################################
  weights.nulldistr <- eval(as.name(compute.evals))(sig1,sig2,sig,act1,act2,act,include.mean)$eval
  weights.nulldistr <- weights.nulldistr[weights.nulldistr!=0]
  if (any(is.na(weights.nulldistr))){
    cat('warning: weight with value NA; pval=NA','\n')
    pval.onesided <- pval.twosided <- NA
  }else{
    if(method.compquadform=='davies'){
      pval.onesided <- davies(teststat,lambda=weights.nulldistr,acc=acc)$Qq;if(show.warn){cat('ifault(davies):',davies(teststat,lambda=weights.nulldistr,acc=acc)$ifault,'\n')}
    }
    if(method.compquadform=='imhof'){
      pval.onesided <-imhof(teststat, lambda=weights.nulldistr,epsabs = epsabs, epsrel = epsrel, limit = 10000)$Qq
      if(pval.onesided<0){pval.onesided <- 0}
      if(pval.onesided>1){pval.onesided <- 1}
    }
    pval.twosided <- 2*min(pval.onesided,1-pval.onesided)
  }
  return(list(pval.onesided=pval.onesided,pval.twosided=pval.twosided,weights.nulldistr=weights.nulldistr,teststat=teststat))
}

##' Differential Network for user specified data splits
##'
##' Remark:
##' 
##' * If include.mean=FALSE, then x1 and x2 have mean zero and DiffNet tests
##'   the hypothesis H0: \Omega_1=\Omega_2. You might need to center x1 and x2.
##' * If include.mean=TRUE, then DiffNet tests the hypothesis
##'   H0: \mu_1=\mu_2 & \Omega_1=\Omega_2
##' * However, we recommend to set include.mean=FALSE and to test equality of the means
##'   separately. 
##' * You might also want to scale x1 and x2, if you are only interested in 
##'   differences due to (partial) correlations.
##'   
##' 
##' @title Differential Network for user specified data splits
##' @param x1 Data-matrix sample 1.
##'           You might need to center and scale your data-matrix.
##' @param x2 Data-matrix sample 2.
##'           You might need to center and scale your data-matrix.
##' @param split1 Samples (condition 1) used in screening step. 
##' @param split2 Samples (condition 2) used in screening step. 
##' @param screen.meth Screening procedure. Options: 'screen_bic.glasso' (default),
##'                    'screen_cv.glasso', 'screen_shrink' (not recommended), 'screen_mb'.
##' @param compute.evals Method to estimate the weights in the weighted-sum-of-chi2s distribution.
##'                      The default and (currently) the only available option 
##'                      is the method 'est2.my.ev3'.
##' @param algorithm.mleggm Algorithm to compute MLE of GGM. The algorithm 'glasso_rho' is the
##'                         default and (currently) the only available option.
##' @param include.mean Should sample specific means be included in hypothesis?
##'                     Use include.mean=FALSE (default and recommended) which assumes mu1=mu2=0
##'                     and tests the hypothesis H0: \Omega_1=\Omega_2.
##' @param method.compquadform Method to compute distribution function of weighted-sum-of-chi2s
##'                            (default='imhof').
##' @param acc See ?davies (default 1e-04).
##' @param epsabs See ?imhof (default 1e-10).
##' @param epsrel See ?imhof (default 1e-10).
##' @param show.warn Should warnings be showed (default=FALSE)?
##' @param save.mle Should MLEs be in the output list (default=FALSE)?
##' @param ... Additional arguments for screen.meth.
##' @return list consisting of
##' \item{pval.onesided}{p-value}
##' \item{pval.twosided}{ignore this output}
##' \item{teststat}{log-likelihood-ratio test statistic}
##' \item{weights.nulldistr}{estimated weights}
##' \item{active}{active-sets obtained in screening-step}
##' \item{sig}{constrained mle (covariance) obtained in cleaning-step}
##' \item{wi}{constrained mle (inverse covariance) obtained in cleaning-step}
##' \item{mu}{mle (mean) obtained in cleaning-step}
##' @author n.stadler
##' @export
diffnet_singlesplit<- function(x1,x2,split1,split2,screen.meth='screen_bic.glasso',
                               compute.evals='est2.my.ev3',algorithm.mleggm='glasso_rho0',include.mean=FALSE,
                               method.compquadform='imhof',acc=1e-04,epsabs=1e-10,epsrel=1e-10,
                               show.warn=FALSE,save.mle=FALSE,...){
  
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  k <- ncol(x1)

  df.param <- k*(k+1)/2
  
  est.mu <- est.sig <- est.wi <- active <- list()#save est.mu;est.sig;est.wi;active
  
  ###############
  ##Joint Model##
  ###############
  xx.train <- rbind(x1[split1,],x2[split2,])
  xx.valid <- rbind(x1[-split1,],x2[-split2,])
  fit.screen <- eval(as.name(screen.meth))(xx.train,include.mean=include.mean,...)
  act <- which(fit.screen$wi[upper.tri(diag(1,k),diag=TRUE)]!=0)
  active[['modJ']] <- act
  if(include.mean==TRUE){
    est.mu[['modJ']] <- colMeans(xx.valid)
  }else{
    est.mu[['modJ']] <- rep(0,k)
  }
  if (length(active[['modJ']])==df.param){
    w <- mcov(x=xx.valid,include.mean=include.mean)
    est.sig[['modJ']] <- w
    est.wi[['modJ']] <- solve(w)
  }
  if (length(active[['modJ']])!=df.param){
    fit.mle <- mle.ggm(xx.valid,fit.screen$wi,algorithm=algorithm.mleggm,rho=fit.screen$rho.opt,include.mean=include.mean)
    w <- fit.mle$w
    est.sig[['modJ']] <- w
    est.wi[['modJ']] <- fit.mle$wi
  }
  #####################  
  ##Individual Models##
  #####################
  for (j in c('1','2')){
    split.train <- eval(as.name(paste('split',j,sep='')))
    xx.train <- eval(as.name(paste('x',j,sep='')))[split.train,]
    xx.valid <- eval(as.name(paste('x',j,sep='')))[-split.train,]
    fit.screen <- eval(as.name(screen.meth))(xx.train,include.mean=include.mean,...)
    act <- which(fit.screen$wi[upper.tri(diag(1,k),diag=TRUE)]!=0)
    active[[paste('modIpop',j,sep='')]] <- act
    if(include.mean==TRUE){
      est.mu[[paste('modIpop',j,sep='')]] <- colMeans(xx.valid)
    }else{
      est.mu[[paste('modIpop',j,sep='')]] <- rep(0,k)
    }
    if (length(active[[paste('modIpop',j,sep='')]])==df.param){
      w <- mcov(xx.valid,include.mean=include.mean)
      est.sig[[paste('modIpop',j,sep='')]] <- w
      est.wi[[paste('modIpop',j,sep='')]] <- solve(w)
    }
    if (length(active[[paste('modIpop',j,sep='')]])!=df.param){
      fit.mle <- mle.ggm(xx.valid,fit.screen$wi,algorithm=algorithm.mleggm,rho=fit.screen$rho.opt,include.mean=include.mean)
      w <- fit.mle$w
      est.sig[[paste('modIpop',j,sep='')]] <- w
      est.wi[[paste('modIpop',j,sep='')]] <- fit.mle$wi
    }
  }

  ###############
  ##Some Checks##
  ###############
  l.act <- lapply(active,length)
  n1.valid <- nrow(x1[-split1,])
  n2.valid <- nrow(x2[-split2,])
  if (any(l.act==0)){ if(show.warn){cat('warning:at least one active-set is empty','\n')}}
  if (all(l.act>= c(n1.valid+n2.valid,n1.valid,n2.valid))){if(show.warn){cat('warning:dim(model) > n-1','\n')}}

  ###########
  ##Pvalue ##
  ###########
  res.pval <- diffnet_pval(x1=x1[-split1,,drop=FALSE],x2=x2[-split2,,drop=FALSE],x=rbind(x1[-split1,,drop=FALSE],x2[-split2,,drop=FALSE]),
                           sig1=est.sig[['modIpop1']],sig2=est.sig[['modIpop2']],sig=est.sig[['modJ']],
                           mu1=est.mu[['modIpop1']],mu2=est.mu[['modIpop2']],mu=est.mu[['modJ']],
                           active[['modIpop1']],active[['modIpop2']],active[['modJ']],
                           compute.evals,include.mean,method.compquadform,acc,epsabs,epsrel,show.warn)
  
  if(save.mle==FALSE){
    est.sig <- est.wi <- est.mu <- NULL
  }
  return(list(pval.onesided=res.pval$pval.onesided,pval.twosided=res.pval$pval.twosided,
              teststat=res.pval$teststat,weights.nulldistr=res.pval$weights.nulldistr,
              active=active,sig=est.sig,wi=est.wi,mu=est.mu))
}

##' Differential Network 
##'
##' Remark:
##' 
##' * If include.mean=FALSE, then x1 and x2 have mean zero and DiffNet tests
##'   the hypothesis H0: \Omega_1=\Omega_2. You might need to center x1 and x2.
##' * If include.mean=TRUE, then DiffNet tests the hypothesis
##'   H0: \mu_1=\mu_2 & \Omega_1=\Omega_2
##' * However, we recommend to set include.mean=FALSE and to test equality of the means
##'   separately. 
##' * You might also want to scale x1 and x2, if you are only interested in 
##'   differences due to (partial) correlations.
##' 
##' 
##' @title Differential Network 
##' @param x1 Data-matrix sample 1.
##'           You might need to center and scale your data-matrix.
##' @param x2 Data-matrix sample 1.
##'           You might need to center and scale your data-matrix.
##' @param b.splits Number of splits (default=50).
##' @param frac.split Fraction train-data (screening) / test-data (cleaning) (default=0.5).
##' @param screen.meth Screening procedure. Options: 'screen_bic.glasso' (default),
##'                    'screen_cv.glasso', 'screen_shrink' (not recommended), 'screen_mb'.
##' @param include.mean Should sample specific means be included in hypothesis?
##'                     Use include.mean=FALSE (default and recommended) which assumes mu1=mu2=0
##'                     and tests the hypothesis H0: \Omega_1=\Omega_2.
##' @param gamma.min Tuning parameter in p-value aggregation of Meinshausen et al (2009). (Default=0.05).
##' @param compute.evals Method to estimate the weights in the weighted-sum-of-chi2s distribution.
##'                      The default and (currently) the only available option 
##'                      is the method 'est2.my.ev3'.
##' @param algorithm.mleggm Algorithm to compute MLE of GGM. The algorithm 'glasso_rho' is the
##'                         default and (currently) the only available option.
##' @param method.compquadform Method to compute distribution function of weighted-sum-of-chi2s
##'                            (default='imhof').
##' @param acc See ?davies (default 1e-04).
##' @param epsabs See ?imhof (default 1e-10).
##' @param epsrel See ?imhof (default 1e-10).
##' @param show.warn Should warnings be showed (default=FALSE)?
##' @param save.mle If TRUE, MLEs (inverse covariance matrices for samples 1 and 2)
##'                 are saved for all b.splits. The median aggregated inverse covariance matrix
##'                 is provided in the output as 'medwi'. The default is save.mle=FALSE.
##' @param ... Additional arguments for screen.meth.
##' @return list consisting of
##' \item{ms.pval}{p-values for all b.splits}
##' \item{ss.pval}{single-split p-value}
##' \item{medagg.pval}{median aggregated p-value}
##' \item{meinshagg.pval}{meinshausen aggregated p-value (meinshausen et al 2009)}
##' \item{teststat}{test statistics for b.splits}
##' \item{weights.nulldistr}{estimated weights}
##' \item{active.last}{active-sets obtained in last screening-step}
##' \item{medwi}{median of inverse covariance matrices over b.splits}
##' \item{sig.last}{constrained mle (covariance matrix) obtained in last cleaning-step}
##' \item{wi.last}{constrained mle (inverse covariance matrix) obtained in last cleaning-step}
##' @author n.stadler
##' @export
##' @example ../../diffnet_ex.R
diffnet_multisplit<- function(x1,x2,b.splits=50,frac.split=1/2,screen.meth='screen_bic.glasso',include.mean=FALSE,
                              gamma.min=0.05,compute.evals='est2.my.ev3',algorithm.mleggm='glasso_rho0',
                              method.compquadform='imhof',acc=1e-04,epsabs=1e-10,epsrel=1e-10,
                              show.warn=FALSE,save.mle=FALSE,...){

  ##????Important Notes: Pval can be NA, because...
  ##????
  ##????                 1) Eval=sqrt(neg. value)=NA
  ##????                 2) W cannot be estimated (|active-set| > n, problems matrix inversion) 
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  
  res.multisplit <- lapply(seq(b.splits),
                           function(i){
                               cat(' split: ',i,'\n\n')
                               split1 <- sample(1:n1,round(n1*frac.split),replace=FALSE)
                               split2 <- sample(1:n2,round(n2*frac.split),replace=FALSE)
                               res.singlesplit <- diffnet_singlesplit(x1,x2,split1,split2,screen.meth,
                                                                      compute.evals,algorithm.mleggm,include.mean,
                                                                      method.compquadform,acc,epsabs,epsrel,show.warn,save.mle,...)
                               
                           })
  pval.onesided <- sapply(res.multisplit,function(x){x[['pval.onesided']]},simplify='array')
  pval.twosided <- sapply(res.multisplit,function(x){x[['pval.twosided']]},simplify='array')
  teststat <- sapply(res.multisplit,function(x){x[['teststat']]},simplify='array')
  weights.nulldistr <- sapply(res.multisplit,function(x){x[['weights.nulldistr']]},simplify='array')
  aggpval.onesided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.onesided[!is.na(pval.onesided)])$objective)
  aggpval.twosided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.twosided[!is.na(pval.twosided)])$objective)
  k <- ncol(x1)
  if(save.mle==TRUE){
    medwi <- list(modJ=matrix(apply(sapply(res.multisplit,function(x){x[['wi']][['modJ']]}),1,median),k,k),
                  modIpop1=matrix(apply(sapply(res.multisplit,function(x){x[['wi']][['modIpop1']]}),1,median),k,k),
                  modIpop2=matrix(apply(sapply(res.multisplit,function(x){x[['wi']][['modIpop2']]}),1,median),k,k))
  }else{
    medwi <- NULL
  }
  
 result <- list(ms.pval=pval.onesided,
                ss.pval=pval.onesided[1],
                medagg.pval=median(pval.onesided,na.rm=TRUE),
                meinshagg.pval=aggpval.onesided,
                teststat=teststat,weights.nulldistr=weights.nulldistr,
                active.last=res.multisplit[[b.splits]]$active,
                medwi=medwi,
                sig.last=res.multisplit[[b.splits]]$sig,wi.last=res.multisplit[[b.splits]]$wi)
  class(result) <- 'diffnet'
  return(result)
}

##' Plotting function for object of class 'diffnet' 
##'
##' 
##' @title Plotting function for object of class 'diffnet' 
##' @param object of class 'diffnet'
##' @return Histogram over multi-split p-values.
##' @author nicolas
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


summary.diffnet <- function(x){
    out <- data.frame(medagg.pval=x$medagg.pval,meinshagg.pval=x$meinshagg.pval)
    rownames(out) <- 'aggregated p-values'
    return(out)
    print(out)
}
