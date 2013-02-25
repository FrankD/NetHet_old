##Packages
library(mvtnorm)
library(glasso)
library(parcor)
library(GeneNet)

#############################
##-------Screening---------##
#############################

##Lambda-grid
lambdagrid_mult <- function(lambda.min,lambda.max,nr.gridpoints){
    mult.const <- (lambda.max/lambda.min)^(1/(nr.gridpoints-1))
    return(lambda.min*mult.const^((nr.gridpoints-1):0))
}

##Crossvalidation Plot
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

##' Crossvalidation for GLasso 
##'
##' 8! lambda-grid has to be increasing (see glassopath)
##' @title Crossvalidation for GLasso
##' @param x 
##' @param folds 
##' @param lambda lambda-grid (increasing!)
##' @param penalize.diagonal 
##' @param plot.it 
##' @param se 
##' @param include.mean 
##' @return 
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

  object <- list(lambda=lambda,residmat=residmat,cv=cv,cv.error=cv.error,w=w,wi=wi,mu=colMeans(x))
  if (plot.it){
    plotCV(lambda,cv,cv.error,se=se)
  }
  invisible(object)
}

lambda.max <- function(x){
  n <- nrow(x)
  s.var <- var(x)
  diag(s.var) <- 0
  return(n*max(abs(s.var))/2)
}

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
    
screen_cv.glasso <- function(x,include.mean=FALSE,
                             folds=10,length.lambda=15,penalize.diagonal=FALSE,
                             trunc.method='linear.growth',trunc.k=5,plot.it=FALSE,se=FALSE)
{
  
  gridmax <- lambda.max(x)
  gridmin <- gridmax/length.lambda
  lambda <- lambdagrid_mult(gridmin,gridmax,length.lambda)[length.lambda:1]
  
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
  wi <- (wi+t(wi))/2
  colnames(w)<-rownames(w)<-colnames(wi)<-rownames(wi)<-colnames(x)
  wi.trunc <- mytrunc.method(n=nrow(x),wi=wi,method=trunc.method,trunc.k=trunc.k)$wi
  
  if (plot.it){
    plotCV(lambda,cv,cv.error,se=se)
  }
  list(wi=wi.trunc,wi.orig=wi)
}

bic.glasso <- function(x,lambda,penalize.diagonal=FALSE,plot.it=TRUE)
{
  ##glasso; lambda.opt with bic
  Mu <- colMeans(x)
  samplecov <- var(x)

  if(is.null(lambda)){
    la <- lambda
  }else{
    la <- 2*lambda/nrow(x)
  }

  fit.path <- glassopath(samplecov,rholist=la,penalize.diagonal=penalize.diagonal,trace=0)

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
                      return(p+(p*(p+1)/2)-n.zero/2)
                    }
                    )
  loglik <- simplify2array(loglik,higher=TRUE)
  degfree <- simplify2array(degfree,higher=TRUE)
  myscore <- -loglik+log(nrow(x))*degfree/2

  if(is.null(lambda)){
    lambda <- 0.5*nrow(x)*fit.path$rholist
  }
  
  if(ncol(x)<nrow(x)){
    loglik.la0 <- sum(dmvnorm(x,mean=Mu,sigma=var(x),log=TRUE))
    degfree.la0 <- ncol(x)+(ncol(x)*(ncol(x)+1)/2)
    myscore.la0 <- -loglik.la0+log(nrow(x))*degfree.la0/2
    myscore <- c(myscore.la0,myscore)
    index.opt <- which.min(myscore)
    if(index.opt==1){
      w <- var(x)
      wi <- solve(var(x))
    }else{
      wi <- fit.path$wi[,,index.opt-1]
      wi[abs(wi)<10^{-3}]<-0
      w <- fit.path$w[,,index.opt-1]
    }
    lambda <- c(0,lambda)
  }else{
    index.opt <- which.min(myscore)
    wi <- fit.path$wi[,,index.opt]
    wi[abs(wi)<10^{-3}]<-0
    w <- fit.path$w[,,index.opt]
  }
  if (plot.it){
    plot(lambda,myscore,type='b',xlab='lambda')
  }
  
  list(lambda=lambda,bic.score=myscore,Mu=Mu,wi=wi,w=w)
}

aic.glasso <- function(x,lambda,penalize.diagonal=FALSE,plot.it=TRUE)
{
  ##glasso; lambda.opt with aicc
  aic.score <-rep(NA,length(lambda))
  Mu <- colMeans(x)
  samplecov <- var(x)

  if(is.null(lambda)){
    la <- lambda
  }else{
    la <- 2*lambda/nrow(x)
  }

  fit.path <- glassopath(samplecov,rholist=la,penalize.diagonal=penalize.diagonal,trace=0)

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
                      return(p+(p*(p+1)/2)-n.zero/2)
                    }
                    )
  loglik <- simplify2array(loglik,higher=TRUE)
  degfree <- simplify2array(degfree,higher=TRUE)
  myscore <- -loglik+2*degfree/2

  if(is.null(lambda)){
    lambda <- 0.5*nrow(x)*fit.path$rholist
  }
  
  if(ncol(x)<nrow(x)){
    loglik.la0 <- sum(dmvnorm(x,mean=Mu,sigma=var(x),log=TRUE))
    degfree.la0 <- ncol(x)+(ncol(x)*(ncol(x)+1)/2)
    myscore.la0 <- -loglik.la0+2*degfree.la0/2
    myscore <- c(myscore.la0,myscore)
    index.opt <- which.min(myscore)
    if(index.opt==1){
      w <- var(x)
      wi <- solve(var(x))
    }else{
      wi <- fit.path$wi[,,index.opt-1]
      wi[abs(wi)<10^{-3}]<-0
      w <- fit.path$w[,,index.opt-1]
    }
    lambda <- c(0,lambda)
  }else{
    index.opt <- which.min(myscore)
    wi <- fit.path$wi[,,index.opt]
    wi[abs(wi)<10^{-3}]<-0
    w <- fit.path$w[,,index.opt]
  }
  if (plot.it){
    plot(lambda,myscore,type='b',xlab='lambda')
  }
  
  list(lambda=lambda,bic.score=myscore,Mu=Mu,wi=wi,w=w)
}


screen_bic.glasso <- function(x,include.mean=TRUE,
                              length.lambda=20,plot.it=FALSE,
                              trunc.method='linear.growth',trunc.k=5){

  gridmax <- lambda.max(x)
  gridmin <- gridmax/length.lambda
  my.grid <- lambdagrid_mult(gridmin,gridmax,length.lambda)[length.lambda:1]
  
  fit.bicgl <- bic.glasso(x,lambda=my.grid,penalize.diagonal=FALSE,plot.it=plot.it)
  wi <- fit.bicgl$wi
  wi.trunc <- mytrunc.method(n=nrow(x),wi=wi,method=trunc.method,trunc.k=trunc.k)$wi
  list(wi=wi.trunc,wi.orig=wi) 
}

screen_aic.glasso <- function(x,include.mean=TRUE,
                              length.lambda=20,plot.it=FALSE,
                              trunc.method='linear.growth',trunc.k=5){

  gridmax <- lambda.max(x)
  gridmin <- gridmax/length.lambda
  my.grid <- lambdagrid_mult(gridmin,gridmax,length.lambda)[length.lambda:1]
  
  fit.aicgl <- aic.glasso(x,lambda=my.grid,penalize.diagonal=FALSE,plot.it=plot.it)
  wi <- fit.aicgl$wi
  wi.trunc <- mytrunc.method(n=nrow(x),wi=wi,method=trunc.method,trunc.k=trunc.k)$wi
  list(wi=wi.trunc,wi.orig=wi) 
}

screen_lasso <- function(x,include.mean=NULL,
                         trunc.method='linear.growth',trunc.k=5){
  
  wi <- adalasso.net(x, k = 10,use.Gram=FALSE,both=FALSE,verbose=FALSE)$pcor.lasso
  wi.trunc <- mytrunc.method(n=nrow(x),wi=wi,method=trunc.method,trunc.k=trunc.k)$wi
  list(wi=wi.trunc,wi.orig=wi)
}

screen_shrink <- function(x,include.mean=NULL,
                          trunc.method='linear.growth',trunc.k=5){
  wi <- ggm.estimate.pcor(x)
  adj <- performance.pcor(wi, fdr=TRUE,verbose=FALSE,plot=FALSE)$adj
  wi[adj==0] <- 0
  wi.trunc <- mytrunc.method(n=nrow(x),wi=wi,method=trunc.method,trunc.k=trunc.k)$wi
  list(wi=wi.trunc,wi.orig=wi)
}

## cv.glasso.approx.trunc <- function(x,include.mean=FALSE,
##                             folds=10,length.lambda=15,penalize.diagonal=FALSE,
##                             trunc.method='linear.growth',trunc.k=5,plot.it=FALSE,se=FALSE)
## {
  
##   gridmax <- lambda.max(x)
##   gridmin <- gridmax/length.lambda
##   lambda <- lambdagrid_mult(gridmin,gridmax,length.lambda)[length.lambda:1]
  
##   colnames(x)<-paste('x',1:ncol(x),sep='')  
##   all.folds <- cv.fold(nrow(x),folds)
##   residmat <- matrix(NA,folds,length(lambda))
##   for (cvfold in 1:folds){
##     omit <- all.folds[[cvfold]]
##     s <- var(x[-omit,])
##     if (include.mean==TRUE){
##       mu <- colMeans(x[-omit,,drop=FALSE])
##     }else{
##       mu <- rep(0,ncol(x))
##     }
##     fit.path <- glassopath(s,rholist=2*lambda/nrow(x),penalize.diagonal=penalize.diagonal,trace=0,approx=TRUE)
##     my.w <- array(apply(fit.path$wi,3,Beta2parcor),dim=dim(fit.path$wi));print(dim(my.w))
##     #my.w <- array(sapply(seq(dim(my.pc)[3]),function(x){solve(my.pc[,,x]*tcrossprod(sqrt(1/diag(s))))}),dim=dim(fit.path$wi));print(dim(my.w))
##     if(length(omit)==1){
##       residmat[cvfold,] <- -2*apply(my.w,3,dmvnorm,log=TRUE,mean=mu,x=x[omit,,drop=FALSE])
##     }else{
##       residmat[cvfold,] <- apply(-2*apply(my.w,3,dmvnorm,log=TRUE,mean=mu,x=x[omit,,drop=FALSE]),2,sum)
##     }
##   }
##   cv <- apply(residmat,2,mean);print(cv)
##   cv.error <- sqrt(apply(residmat,2,var)/folds)
##   gl.opt<-glasso(var(x),rho=2*lambda[which.min(cv)]/nrow(x),penalize.diagonal=penalize.diagonal)
##   cat('la.opt:',lambda[which.min(cv)],'\n')
##   w<-gl.opt$w
##   wi<-gl.opt$wi
##   wi[abs(wi)<10^{-3}]<-0
##   wi <- (wi+t(wi))/2
##   colnames(w)<-rownames(w)<-colnames(wi)<-rownames(wi)<-colnames(x)
##   wi <- mytrunc.method(n=nrow(x),wi=wi,method=trunc.method,trunc.k=trunc.k)$wi
  
##   if (plot.it){
##     plotCV(lambda,cv,cv.error,se=se)
##   }
##   list(wi=wi)
## }


#screen_full <- function(x,include.mean=NULL,trunc.method='linear.growth',trunc.k=5){
#  wi <- diag(1,ncol(x))
#  list(wi=wi)
#}

#aicc.glasso <- function(x,lambda,penalize.diagonal=FALSE,plot.it=TRUE)
#{
#  ##glasso; lambda.opt with aicc
#  aic.score <-rep(NA,length(lambda))
#  Mu <- colMeans(x)
#  samplecov <- var(x)
#
#  if(is.null(lambda)){
#    la <- lambda
#  }else{
#    la <- 2*lambda/nrow(x)
#  }
#
#  fit.path <- glassopath(samplecov,rholist=la,penalize.diagonal=penalize.diagonal,trace=0)
#
#  loglik <- lapply(seq(length(fit.path$rholist)),
#                   function(i){
#                     return(sum(dmvnorm(x,mean=Mu,sigma=fit.path$w[,,i],log=TRUE)))
#                   }
#                   )
#  degfree <- lapply(seq(length(fit.path$rholist)),
#                    function(i){
#                      wi <- fit.path$wi[,,i]
#                      wi[abs(wi)<10^{-3}]<-0
#                      p <- ncol(wi)
#                      n.zero<-sum(wi==0)
#                      return(p+(p*(p+1)/2)-n.zero/2)
#                    }
#                    )
#  loglik <- simplify2array(loglik,higher=TRUE)
#  degfree <- simplify2array(degfree,higher=TRUE)
#  myscore <- -loglik+2*degfree/2+(degfree*(degfree+1)/(nrow(x)-degfree-1))
#
#  if(is.null(lambda)){
#    lambda <- 0.5*nrow(x)*fit.path$rholist
#  }
#  
#  if(ncol(x)<nrow(x)){
#    loglik.la0 <- sum(dmvnorm(x,mean=Mu,sigma=var(x),log=TRUE))
#    degfree.la0 <- ncol(x)+(ncol(x)*(ncol(x)+1)/2)
#    myscore.la0 <- -loglik.la0+2*degfree.la0/2+(degfree.la0*(degfree.la0+1)/(nrow(x)-degfree.la0-1))
#    myscore <- c(myscore.la0,myscore)
#    index.opt <- which.min(myscore)
#    if(index.opt==1){
#      w <- var(x)
#      wi <- solve(var(x))
#    }else{
#      wi <- fit.path$wi[,,index.opt-1]
#      wi[abs(wi)<10^{-3}]<-0
#      w <- fit.path$w[,,index.opt-1]
#    }
#    lambda <- c(0,lambda)
#  }else{
#    index.opt <- which.min(myscore)
#    wi <- fit.path$wi[,,index.opt]
#    wi[abs(wi)<10^{-3}]<-0
#    w <- fit.path$w[,,index.opt]
#  }
#  if (plot.it){
#    plot(lambda,myscore,type='b',xlab='lambda')
#  }
#  
#  list(lambda=lambda,bic.score=myscore,Mu=Mu,wi=wi,w=w)
#}
