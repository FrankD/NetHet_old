###TwoSampleTest: based on Permutation-Test
###Date: 01/08/2012
###
###Changes: -add additional option: lambda.cv={cv.glasso,cv.glasso.1se}

##Packages
library(mvtnorm)
library(glasso)
library(glmnet)

###########
##GLasso ##
###########

##Lambda-grid
lambdagrid_mult <- function(lambda.min,lambda.max,nr.gridpoints){
    mult.const <- (lambda.max/lambda.min)^(1/(nr.gridpoints-1))
    return(lambda.min*mult.const^((nr.gridpoints-1):0))
}

##Crossvalidation glasso
plotCV <- function(lambda,cv,cv.error,se=TRUE,type='b',...){
  if (se==TRUE){ylim <- range(cv,cv+cv.error,cv-cv.error)}
  else {ylim <- range(cv)}
  plot(lambda,cv,type=type,ylim=ylim,...)
  if (se)
    error.bars(lambda,cv+cv.error,cv-cv.error,width=1/length(lambda))
  invisible()
}

error.bars <- function (x, upper, lower, width = 0.02, ...)
{
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}

cv.fold <- function(n,folds=10){
  split(sample(1:n),rep(1:folds,length=n))
}

cv.glasso <- function(x,folds=10,lambda,penalize.diagonal=FALSE,plot.it=FALSE,se=TRUE)
##8! lambda has to be increasing for glassopath
{
  colnames(x)<-paste('x',1:ncol(x),sep='')  
  all.folds <- cv.fold(nrow(x),folds)
  residmat <- matrix(NA,folds,length(lambda))
  for (cvfold in 1:folds){
    omit <- all.folds[[cvfold]]
    s <- var(x[-omit,])###!!!!!!!!!!!!!
    fit.path <- glassopath(s,rholist=2*lambda/nrow(x),penalize.diagonal=penalize.diagonal,trace=0)
    if(length(omit)==1){
      residmat[cvfold,] <- -2*apply(fit.path$w,3,dmvnorm,log=TRUE,mean=rep(0,ncol(x)),x=x[omit,,drop=FALSE])
    }else{
      residmat[cvfold,] <- apply(-2*apply(fit.path$w,3,dmvnorm,log=TRUE,mean=rep(0,ncol(x)),x=x[omit,,drop=FALSE]),2,sum)
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

cv.glasso.1se <- function(x,folds=10,lambda,penalize.diagonal=FALSE,plot.it=FALSE,se=TRUE)
##8! lambda has to be increasing for glassopath
{
  colnames(x)<-paste('x',1:ncol(x),sep='')  
  all.folds <- cv.fold(nrow(x),folds)
  residmat <- matrix(NA,folds,length(lambda))

  for (cvfold in 1:folds){
    omit <- all.folds[[cvfold]]
    s <- var(x[-omit,])###!!!!!!!!!!!!!
    fit.path <- glassopath(s,rholist=2*lambda/nrow(x),penalize.diagonal=penalize.diagonal,trace=0)
    if(length(omit)==1){
      residmat[cvfold,] <- -2*apply(fit.path$w,3,dmvnorm,log=TRUE,mean=rep(0,ncol(x)),x=x[omit,,drop=FALSE])
    }else{
      residmat[cvfold,] <- apply(-2*apply(fit.path$w,3,dmvnorm,log=TRUE,mean=rep(0,ncol(x)),x=x[omit,,drop=FALSE]),2,sum)
    }
  }
  cv <- apply(residmat,2,mean)
  cv.error <- sqrt(apply(residmat,2,var)/folds)
  la.1se <- max(lambda[which(cv < min(cv)+cv.error[which.min(cv)])])
  gl.opt<-glasso(var(x),rho=2*la.1se/nrow(x),penalize.diagonal=penalize.diagonal)
  w<-gl.opt$w
  wi<-gl.opt$wi
  wi[abs(wi)<10^{-3}]<-0
  colnames(w)<-rownames(w)<-colnames(wi)<-rownames(wi)<-colnames(x)  

  object <- list(lambda=lambda,residmat=residmat,cv=cv,cv.error=cv.error,w=w,wi=wi,mu=colMeans(x),la.1se=la.1se)
  if (plot.it){
    plotCV(lambda,cv,cv.error,se=se)
  }
  invisible(object)
}

## glasso.parcor.launi <- function(x,maxiter=1000,term=10^{-3}){
##   s <- var(x)
##   rho.uni <- sqrt(2*log(ncol(x))/nrow(x))

##   ww <- rep(1,p)#sqrt(diag(s))
##   iter <- 0
##   err <- Inf #convergence of parameters
##   param <- as.vector(diag(ww))
##   while((err>term)&(iter<maxiter)){
##     gl <- glasso(s,rho=rho.uni*ww,penalize.diagonal=FALSE)
##     ww <- 1/(diag(gl$wi))
##     param.old <- param
##     param <- as.vector(gl$w)
##     err <- max(abs(param-param.old)/(1+abs(param)))
##     iter <- iter+1
##   }
##   list(w=gl$w,wi=gl$wi,iter=iter)
## }

## glasso.launi <- function(x){
##   s <- var(x)
##   rho.uni <- sqrt(2*log(ncol(x))/nrow(x))
##   gl <- glasso(s,rho=rho.uni,penalize.diagonal=FALSE)
##   list(w=gl$w,wi=gl$wi)
## }


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

perm.regr.pval <- function(y1,y2,x1,x2,lambda.cv='lambda.min',nr.perm=100){
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
  return(list(pval.kldist=mean(tperm.kldist>tobs.kldist)))
}
              
##DiffNet
tr <- function(m){sum(diag(m))}

symmkldistmvn <- function(mu1,mu2,sig1,sig2){
    symmkl <- 0.5*tr((sig1-sig2)%*%(solve(sig2)-solve(sig1)))+0.5*t(mu1-mu2)%*%(solve(sig1)+solve(sig2))%*%(mu1-mu2)
    return(as.numeric(symmkl))
}

perm.diffnet.teststatistic <- function(x,n1,n2,folds=10,lambda,lambda.cv='cv.glasso',plot.it=TRUE,se=TRUE){
  
  fit.cv.1<- eval(as.name(lambda.cv))(x[1:n1,],folds=folds,lambda=lambda,plot.it=plot.it,se=se)
  fit.cv.2 <- eval(as.name(lambda.cv))(x[n1+(1:n2),],folds=folds,lambda=lambda,plot.it=plot.it,se=se)
  nz1 <- fit.cv.1$wi;nz1[nz1!=0] <- 1
  nz2 <- fit.cv.2$wi;nz2[nz2!=0] <- 1
  kldist <- symmkldistmvn(mu1=fit.cv.1$mu,mu2=fit.cv.2$mu,sig1=fit.cv.1$w,sig2=fit.cv.2$w)
  return(list(edgediff=sum(nz1!=nz2)/2,kldist=kldist))
}

perm.diffnet.pval <- function(x1,x2,nr.perm=100,lambda,lambda.cv='cv.glasso'){
  n1 <- nrow(x1);n2 <- nrow(x2)
  x <- rbind(x1,x2)
  fit.tobs <- perm.diffnet.teststatistic(x,n1,n2,lambda=lambda,lambda.cv=lambda.cv)
  tobs.edgediff <- fit.tobs$edgediff
  tobs.kldist <- fit.tobs$kldist
  tperm.edgediff <- tperm.kldist<- rep(NA,nr.perm)
  for (i in 1:nr.perm){
    x.perm <- x[sample(1:(n1+n2)),]
    fit.tperm <-  perm.diffnet.teststatistic(x.perm,n1,n2,lambda=lambda,lambda.cv=lambda.cv)
    tperm.edgediff[i] <-fit.tperm$edgediff
    tperm.kldist[i] <- fit.tperm$kldist
  }
  return(list(pval.edgediff=mean(tperm.edgediff>tobs.edgediff),
              pval.kldist=mean(tperm.kldist>tobs.kldist)))
}
              

