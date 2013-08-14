###TwoSampleTest: based on Permutation-Test
###Date: 01/08/2012
###
###Changes: -add additional option: lambda.cv={cv.glasso,cv.glasso.1se}

##Packages
library(mvtnorm)
library(glasso)
library(glmnet)
library(parcor)
library(GeneNet)
library(huge)

#############################
##-------Screening---------##
#############################

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

screen_cv.glasso <- function(x,include.mean=FALSE,covMethod=NULL,
                             folds=10,length.lambda=20,lambdamin.ratio=ifelse(ncol(x)>nrow(x),0.01,0.001),penalize.diagonal=FALSE,
                             trunc.method='none',trunc.k=5,plot.it=FALSE,se=FALSE,use.package='huge',verbose=TRUE)
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
  
  if (include.mean==TRUE){
    mu <- colMeans(x[,,drop=FALSE])
  }else{
    mu <- rep(0,ncol(x))
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
  list(rho.opt=2*lambda[which.min(cv)]/nrow(x),mu=mu,w=w,wi=wi.trunc,wi.orig=wi)
}

####################
##Permutation-Test##
####################

##Regr
kldistregr <- function(beta1,beta2,sd1,sd2,Sig){
  return(-0.5+0.5*(sd1^2)/(sd2^2)+t(beta1-beta2)%*%Sig%*%(beta1-beta2)/(2*sd2^2)+0.5*log((sd2^2)/(sd1^2)))
}

perm.regr.teststatistic <- function(y,x,n1,n2,lambda.cv){

  x1 <- x[1:n1,]
  y1 <- y[1:n1]
  fit.cv.1<- cv.glmnet(x1,y1)
  beta1 <- as.numeric(coef(fit.cv.1,s=lambda.cv))
  mu1<-as.numeric(x1%*%beta1[-1]+beta1[1])
  sig1 <- (sum((y1-mu1)^2)/n1)

  x2 <- x[n1+(1:n2),]
  y2 <- y[n1+(1:n2)]
  fit.cv.2 <- cv.glmnet(x2,y2)
  beta2 <- as.numeric(coef(fit.cv.2,s=lambda.cv))
  mu2<-as.numeric(x2%*%beta2[-1]+beta2[1])
  sig2 <- (sum((y2-mu2)^2)/n2)
  
  kldist <- kldistregr(beta1,beta2,sqrt(sig1),sqrt(sig2),var(cbind(rep(1,n1),x1)))+kldistregr(beta2,beta1,sqrt(sig2),sqrt(sig1),var(cbind(rep(1,n2),x2)))
  return(list(kldist=as.numeric(kldist)))
}

perm.regr.pval <- function(y1,y2,x1,x2,lambda.cv='lambda.min',nr.perm){
  n1 <- nrow(x1);n2 <- nrow(x2)
  x <- rbind(x1,x2)
  y <- c(y1,y2)
  tobs.kldist <- perm.regr.teststatistic(y,x,n1,n2,lambda.cv=lambda.cv)$kldist
  tperm.kldist<- rep(NA,nr.perm)
  for (i in 1:nr.perm){
    perm <- sample(1:(n1+n2))
    x.perm <- x[perm,]
    y.perm <- y[perm]
    tperm.kldist[i] <- perm.regr.teststatistic(y.perm,x.perm,n1,n2,lambda.cv=lambda.cv)$kldist
  }
  return(list(pval.kldist=(1+sum(tperm.kldist>=tobs.kldist))/nr.perm))
}
              
##DiffNet
tr <- function(m){sum(diag(m))}

symmkldistmvn <- function(mu1,mu2,sig1,sig2){
    symmkl <- 0.5*tr((sig1-sig2)%*%(solve(sig2)-solve(sig1)))+0.5*t(mu1-mu2)%*%(solve(sig1)+solve(sig2))%*%(mu1-mu2)
    return(as.numeric(symmkl))
}

perm.diffnet.teststatistic <- function(x,n1,n2,screen.meth='screen_cv.glasso',plot.it=TRUE,se=TRUE,include.mean=FALSE){
  
  fit.cv.1<- eval(as.name(screen.meth))(x[1:n1,],include.mean=include.mean)
  fit.cv.2 <- eval(as.name(screen.meth))(x[n1+(1:n2),],include.mean=include.mean)
  nz1 <- fit.cv.1$wi;nz1[nz1!=0] <- 1
  nz2 <- fit.cv.2$wi;nz2[nz2!=0] <- 1
  kldist <- symmkldistmvn(mu1=fit.cv.1$mu,mu2=fit.cv.2$mu,sig1=fit.cv.1$w,sig2=fit.cv.2$w)
  return(list(edgediff=sum(nz1!=nz2)/2,kldist=kldist))
}

perm.diffnet.pval <- function(x1,x2,nr.perm,screen.meth='screen_cv.glasso',include.mean=FALSE){
  n1 <- nrow(x1);n2 <- nrow(x2)
  x <- rbind(x1,x2)
  fit.tobs <- perm.diffnet.teststatistic(x,n1,n2,screen.meth=screen.meth,include.mean=include.mean)
  tobs.edgediff <- fit.tobs$edgediff
  tobs.kldist <- fit.tobs$kldist
  tperm.edgediff <- tperm.kldist<- rep(NA,nr.perm)
  for (i in 1:nr.perm){
    x.perm <- x[sample(1:(n1+n2)),]
    fit.tperm <-  perm.diffnet.teststatistic(x.perm,n1,n2,screen.meth=screen.meth,include.mean=include.mean)
    tperm.edgediff[i] <-fit.tperm$edgediff
    tperm.kldist[i] <- fit.tperm$kldist
  }
  return(list(pval.edgediff=(1+sum(tperm.edgediff>=tobs.edgediff))/nr.perm,
              pval.kldist=(1+sum(tperm.kldist>=tobs.kldist))/nr.perm))
}
              
