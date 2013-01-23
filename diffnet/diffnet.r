###High-Dimensional Two-Sample Testing for Gaussian Graphical Model
###Date: 27/11/2012
###
###Changes: - align with twosample_highdimregr-07062012.r
###         - add function for multi-split pvals
###         - 07062012: additional option diag.invcov=TRUE/FALSE
###         - 12062012: oracle-versions of multisplit-pvals
###         - 15062012: set Pval=NA if "p>n" (changes only in est.my.ev2)
###         - 21062012: new function: beta.mat, q.matrix3, est2.ww.mat2, est2.my.ev2
###         - 25062012: twosample_diffnet2 with additional glasso.launi option
###         - 12072012: small change in cv.glasso function
###         - 20072012: cleaning up and get rid of old functions; single-split method added on 13082012
###         - 27112012: diffnet.r originates from twosample_diffnet-20072012.R

##Packages
library(mvtnorm)
library(glasso)
library(CompQuadForm)

##load C-code
#dyn.load("../code/betamat_diffnet.so")

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
cv.glasso <- function(x,folds=10,lambda,penalize.diagonal=FALSE,plot.it=FALSE,se=TRUE,include.mean=TRUE)
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

cvtrunc.glasso <- function(x,folds=10,lambda,trunc.k=5,penalize.diagonal=FALSE,include.mean=TRUE)
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
  wi <- (wi+t(wi))/2
  colnames(w)<-rownames(w)<-colnames(wi)<-rownames(wi)<-colnames(x)

  wi.trunc <- wi
  diag(wi.trunc) <- 0
  nonzero <- min(2*ceiling(ncol(x)*ceiling(nrow(x)/trunc.k)/2),sum(wi.trunc!=0))
  wi.trunc[-order(abs(wi.trunc),decreasing=TRUE)[1:nonzero]] <- 0
  diag(wi.trunc) <- diag(wi)

  list(wi=wi.trunc)
  
}

##' Crossvalidation for GLasso (1-se rule)
##'
##' 8! lambda-grid has to be increasing (see glassopath)
##' @title Crossvalidation for GLasso (1-se rule)
##' @param x 
##' @param folds 
##' @param lambda lambda-grid (increasing!)
##' @param penalize.diagonal 
##' @param plot.it 
##' @param se 
##' @param include.mean 
##' @return 
##' @author n.stadler
cv.glasso.1se <- function(x,folds=10,lambda,penalize.diagonal=FALSE,plot.it=FALSE,se=TRUE,include.mean=TRUE)
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

glasso.parcor.launi <- function(x,maxiter=1000,term=10^{-3},include.mean=NULL,lambda=NULL){
  p <- ncol(x)
  s <- var(x)
  rho.uni <- sqrt(2*log(ncol(x))/nrow(x))

  ww <- rep(1,p)#sqrt(diag(s))
  iter <- 0
  err <- Inf #convergence of parameters
  param <- as.vector(diag(ww))
  while((err>term)&(iter<maxiter)){
    gl <- glasso(s,rho=rho.uni*ww,penalize.diagonal=FALSE)
    ww <- 1/(diag(gl$wi))
    param.old <- param
    param <- as.vector(gl$w)
    err <- max(abs(param-param.old)/(1+abs(param)))
    iter <- iter+1
  }
  list(w=gl$w,wi=gl$wi,mu=colMeans(x),iter=iter)
}

glasso.parcor.launi.trunc <- function(x,trunc.k=5,maxiter=1000,term=10^{-3},include.mean=NULL,lambda=NULL){
  p <- ncol(x)
  s <- var(x)
  rho.uni <- sqrt(2*log(ncol(x))/nrow(x))

  ww <- rep(1,p)#sqrt(diag(s))
  iter <- 0
  err <- Inf #convergence of parameters
  param <- as.vector(diag(ww))
  while((err>term)&(iter<maxiter)){
    gl <- glasso(s,rho=rho.uni*ww,penalize.diagonal=FALSE)
    ww <- 1/(diag(gl$wi))
    param.old <- param
    param <- as.vector(gl$w)
    err <- max(abs(param-param.old)/(1+abs(param)))
    iter <- iter+1
  }
  wi <- gl$wi
  wi <- (wi+t(wi))/2

  wi.trunc <- wi
  diag(wi.trunc) <- 0
  nonzero <- min(2*ceiling(ncol(x)*ceiling(nrow(x)/trunc.k)/2),sum(wi.trunc!=0))
  wi.trunc[-order(abs(wi.trunc),decreasing=TRUE)[1:nonzero]] <- 0
  diag(wi.trunc) <- diag(wi)

  list(wi=wi.trunc)
 
}



##############################
##--------P-VALUES----------##
##############################
##' LogLikelihood-Ratio 
##'
##' .. content for \details{} ..
##' @title LogLikelihood-Ratio
##' @param x1 
##' @param x2 
##' @param x 
##' @param sig1 
##' @param sig2 
##' @param sig 
##' @param mu1 
##' @param mu2 
##' @param mu 
##' @return 
##' @author n.stadler
logratio <- function(x1,x2,x,sig1,sig2,sig,mu1,mu2,mu){
  twiceLR <- 2*(sum(dmvnorm(x1,mean=mu1,sigma=sig1,log=TRUE))+sum(dmvnorm(x2,mean=mu2,sigma=sig2,log=TRUE))-sum(dmvnorm(x,mean=mu,sigma=sig,log=TRUE)))
  list(twiceLR=twiceLR,sig1=sig1,sig2=sig2,sig=sig)
}

##' Compute Information Matrix of Gaussian Graphical Model
##'
##' computes E_0[s(Y;Omega)s(Y;Omega)'] where s(Y;Omega)=(d/dOmega) LogLik
##' @title Information Matrix of Gaussian Graphical Model
##' @param Sig Sig=solve(SigInv) true covariance matrix under H0
##' @param include.mean 
##' @return 
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
##' @param imat 
##' @param act I_uv
##' @param act1 I_u
##' @param act2 I_v
##' @return 
##' @author n.stadler
ww.mat <- function(imat,act,act1,act2){
 
  bfg <- rbind(imat[act1,act,drop=FALSE],imat[act2,act,drop=FALSE])
  bf <- matrix(0,length(act1)+length(act2),length(act1)+length(act2))
  bf[1:length(act1),1:length(act1)] <- imat[act1,act1]
  bf[length(act1)+(1:(length(act2))),length(act1)+(1:(length(act2)))] <- imat[act2,act2]
  bg <- 2*imat[act,act,drop=FALSE]
  mat <- rbind(cbind(diag(1,length(act1)+length(act2)),bfg%*%solve(bg)),cbind(-t(bfg)%*%solve(bf),diag(-1,length(act))))
  eval <- as.real(eigen(mat)$values)
  eval[abs(eval)<10^{-6}] <- 0
  
  return(list(ww.mat=mat,eval=eval))
}

##' Calculates eigenvalues of weight-matrix (using 1st order simplification)
##'
##' calculation based on true information matrix
##' @title 
##' @param imat 
##' @param act I_uv
##' @param act1 I_u
##' @param act2 I_v
##' @return 
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
  eval.mu <- as.real(eigen(mat)$values)
  eval2 <- 1-eval.mu
  eval2[abs(eval2)<10^{-6}] <- 0
  eval <- c(eval,sqrt(eval2),-sqrt(eval2))
  return(list(ww.mat=mat,eval=eval))
}

##' Compute beta-matrix 
##'
##' beta-matrix=E[s_ind1(Y;sig1) s_ind2(Y;sig2)'|sig]
##' @title Compute beta-matrix 
##' @param ind1 
##' @param ind2 
##' @param sig1 
##' @param sig2 
##' @param sig 
##' @return 
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
##' .. content for \details{} ..
##' @title Compute Q-matrix 
##' @param sig 
##' @param sig.a 
##' @param sig.b 
##' @param act.a 
##' @param act.b 
##' @param ss 
##' @return 
##' @author n.stadler
q.matrix3 <- function(sig,sig.a,sig.b,act.a,act.b,ss){
  b.ab<-beta.mat(act.a,act.b,sig.a,sig.b,sig)
  aa<-seq(1,length(act.a))[!(act.a%in%ss)]
  bb<-seq(1,length(act.b))[!(act.b%in%ss)]
  s.a<-seq(1,length(act.a))[(act.a%in%ss)]
  s.b<-seq(1,length(act.b))[(act.b%in%ss)]
  
  if(length(ss)==0){
    return(b.ab[aa,bb])
  }else{
    return(b.ab[aa,bb,drop=FALSE]-(b.ab[aa,s.b,drop=FALSE]%*%solve(b.ab[s.a,s.b,drop=FALSE])%*%b.ab[s.a,bb,drop=FALSE]))
  }
}

##' Compute weights of sum-w-chi2 (1st order simplification)
##'
##' .. content for \details{} ..
##' @title Weights of sum-w-chi2
##' @param sig1 
##' @param sig2 
##' @param sig 
##' @param act1 
##' @param act2 
##' @param act 
##' @param include.mean 
##' @return 
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
  onemineval.mu <- 1-as.real(eval.mu.complex)
  onemineval.mu[abs(onemineval.mu)<10^{-6}]<-0
  eval <- c(eval,sqrt(onemineval.mu),-sqrt(onemineval.mu))
  if(include.mean==TRUE){
    eval <- c(rep(0,2*k),eval)
  }
  return(list(ww.mat=mat,eval=eval,eval.mu.complex=eval.mu.complex))
}


##' Compute weights of sum-w-chi2 (2nd order simplification)
##'
##' .. content for \details{} ..
##' @title Weights of sum-w-chi2
##' @param sig1 
##' @param sig2 
##' @param sig 
##' @param act1 
##' @param act2 
##' @param act 
##' @param include.mean 
##' @return 
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
      ev.aux <- as.real(ev.aux.complex)
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
      ev.aux <- as.real(ev.aux.complex)
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
   
##' P-value aggregation
##'
##' .. content for \details{} ..
##' @title P-value aggregation
##' @param gamma 
##' @param pval 
##' @return 
##' @author n.stadler
agg.pval <- function(gamma,pval){
    min(quantile(pval/gamma,probs=gamma),1)
}

##' P-value calculation
##'
##' .. content for \details{} ..
##' @title P-value calculation
##' @param x1 
##' @param x2 
##' @param x 
##' @param sig1 
##' @param sig2 
##' @param sig 
##' @param mu1 
##' @param mu2 
##' @param mu 
##' @param act1 
##' @param act2 
##' @param act 
##' @param compute.evals 
##' @param include.mean 
##' @param acc 
##' @return 
##' @author n.stadler
diffnet_pval <- function(x1,x2,x,sig1,sig2,sig,mu1,mu2,mu,act1,act2,act,compute.evals,include.mean,acc,show.trace){
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
  }else{
    pval.onesided <- davies(teststat,lambda=weights.nulldistr,acc=acc)$Qq;if(show.trace){cat('ifault(davies):',davies(teststat,lambda=weights.nulldistr,acc=acc)$ifault,'\n')}
    pval.twosided <- 2*min(pval.onesided,1-pval.onesided)
  }
  return(list(pval.onesided=pval.onesided,pval.twosided=pval.twosided,weights.nulldistr=weights.nulldistr,teststat=teststat))
}

##' Differential network for given split
##'
##' if include.mean=FALSE then x1 and x2 have to be scaled to var=1 & mean=0
##' @title DiffNet for given split
##' @param x1 
##' @param x2 
##' @param split1 
##' @param split2 
##' @param screen.meth 
##' @param compute.evals 
##' @param include.mean 
##' @param diag.invcov 
##' @param acc 
##' @param ... 
##' @return 
##' @author n.stadler
diffnet_singlesplit<- function(x1,x2,split1,split2,screen.meth='cv.glasso',
                               compute.evals='est2.my.ev2',include.mean=TRUE,
                               diag.invcov=TRUE,acc=1e-04,show.trace=FALSE,...){
  
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  k <- ncol(x1)

  if(diag.invcov){df.param <- k*(k+1)/2}
  if(!diag.invcov){df.param <- k*(k-1)/2}

  est.mu <- est.sig <- est.wi <- active <- list()#save est.mu;est.sig;est.wi;active
  
  ###############
  ##Joint Model##
  ###############
  xx.train <- rbind(x1[split1,],x2[split2,])
  xx.valid <- rbind(x1[-split1,],x2[-split2,])
  fit.screen <- eval(as.name(screen.meth))(xx.train,include.mean=include.mean,...)
  act <- which(fit.screen$wi[upper.tri(diag(1,k),diag=TRUE)]!=0)
  if(!diag.invcov){
    ind.diag <- which(diag(1,k)[upper.tri(diag(1,k),diag=TRUE)]!=0)
    act <- setdiff(act,ind.diag)
  }
  active[['modJ']] <- act
  if(include.mean==TRUE){
    est.mu[['modJ']] <- colMeans(xx.valid)
  }else{
    est.mu[['modJ']] <- rep(0,k)
  }
  if (length(active[['modJ']])==df.param){
    w <- var(xx.valid)
    if(!diag.invcov){
      w <-w*sqrt(tcrossprod(diag(solve(w))))
    }
    est.sig[['modJ']] <- w
    est.wi[['modJ']] <- solve(w)
  }
  if (length(active[['modJ']])!=df.param){
    fit.mle <- glasso(var(xx.valid),rho=10^{-10},zero=which(fit.screen$wi==0,arr.in=TRUE))
    w <- fit.mle$w
    if(!diag.invcov){
      w <-w*sqrt(tcrossprod(diag(solve(fit.mle$w))))
    }
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
    if(!diag.invcov){
      ind.diag <- which(diag(1,k)[upper.tri(diag(1,k),diag=TRUE)]!=0)
      act <- setdiff(act,ind.diag)
    }
    active[[paste('modIpop',j,sep='')]] <- act
    if(include.mean==TRUE){
      est.mu[[paste('modIpop',j,sep='')]] <- colMeans(xx.valid)
    }else{
      est.mu[[paste('modIpop',j,sep='')]] <- rep(0,k)
    }
    if (length(active[[paste('modIpop',j,sep='')]])==df.param){
      w <- var(xx.valid)
      if(!diag.invcov){
        w <-w*sqrt(tcrossprod(diag(solve(w))))
      }
      est.sig[[paste('modIpop',j,sep='')]] <- w
      est.wi[[paste('modIpop',j,sep='')]] <- solve(w)
    }
    if (length(active[[paste('modIpop',j,sep='')]])!=df.param){
      fit.mle <- glasso(var(xx.valid),rho=10^{-10},zero=which(fit.screen$wi==0,arr.in=TRUE))
      w <- fit.mle$w
      if(!diag.invcov){
        w <-w*sqrt(tcrossprod(diag(solve(fit.mle$w))))
      }
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
  if (any(l.act==0)){ if(show.trace){cat('warning:at least one active-set is empty','\n')}}
  if (all(l.act>= c(n1.valid+n2.valid,n1.valid,n2.valid))){if(show.trace){cat('warning:dim(model) > n-1','\n')}}

  ###########
  ##Pvalue ##
  ###########
  res.pval <- diffnet_pval(x1=x1[-split1,,drop=FALSE],x2=x2[-split2,,drop=FALSE],x=rbind(x1[-split1,,drop=FALSE],x2[-split2,,drop=FALSE]),
                           sig1=est.sig[['modIpop1']],sig2=est.sig[['modIpop2']],sig=est.sig[['modJ']],
                           mu1=est.mu[['modIpop1']],mu2=est.mu[['modIpop2']],mu=est.mu[['modJ']],
                           active[['modIpop1']],active[['modIpop2']],active[['modJ']],
                           compute.evals,include.mean,acc,show.trace)
  

  return(list(pval.onesided=res.pval$pval.onesided,pval.twosided=res.pval$pval.twosided,
              teststat=res.pval$teststat,weights.nulldistr=res.pval$weights.nulldistr,
              active=active,sig=est.sig,wi=est.wi,mu=est.mu))
}

##' Differential Network (mulitsplit)
##'
##' if include.mean=FALSE then x1 and x2 have to be scaled to var=1 & mean=0
##' @title DiffNet
##' @param x1 if include.mean=FALSE: scale data (var=1,mean=0) !
##' @param x2 if include.mean=FALSE: scale data (var=1,mean=0) !
##' @param b.splits number of splits
##' @param frac.split  fraction train-data (screening) / test-data (cleaning)
##' @param screen.meth 
##' @param include.mean 
##' @param gamma.min see p-value aggregation (Meinshausen&Meier&Buehlmann)
##' @param compute.evals 'est2.my.ev2'/'est2.ww.mat2 '
##' @param diag.invcov  TRUE (parameter of interest is invcov; diag(invcov) also considered) / FALSE (param of interest is par.cor; no diagonal)
##' @param acc 
##' @param ... additional arguments for screen.meth
##' @return 
##' @author n.stadler
diffnet_multisplit<- function(x1,x2,b.splits=50,frac.split=1/2,screen.meth='cv.glasso',include.mean=FALSE,
                              gamma.min=0.05,compute.evals='est2.my.ev2',diag.invcov=TRUE,acc=1e-04,show.trace=FALSE,...){

  ##????Important Notes: Pval can be NA, because...
  ##????
  ##????                 1) Eval=sqrt(neg. value)=NA
  ##????                 2) W cannot be estimated (|active-set| > n, problems matrix inversion) 
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  
  res.multisplit <- lapply(seq(b.splits),
                          function(i){
                            split1 <- sample(1:n1,round(n1*frac.split),replace=FALSE)
                            split2 <- sample(1:n2,round(n2*frac.split),replace=FALSE)
                            res.singlesplit <- diffnet_singlesplit(x1,x2,split1,split2,screen.meth,
                                                                   compute.evals,include.mean,diag.invcov,acc,show.trace,...)
                            
                          })
  pval.onesided <- sapply(res.multisplit,function(x){x[['pval.onesided']]},simplify='array')
  pval.twosided <- sapply(res.multisplit,function(x){x[['pval.twosided']]},simplify='array')
  teststat <- sapply(res.multisplit,function(x){x[['teststat']]},simplify='array')
  weights.nulldistr <- sapply(res.multisplit,function(x){x[['weights.nulldistr']]},simplify='array')
  aggpval.onesided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.onesided[!is.na(pval.onesided)])$objective)
  aggpval.twosided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.twosided[!is.na(pval.twosided)])$objective)
    
  return(list(pval.onesided=pval.onesided,pval.twosided=pval.twosided,
              sspval.onesided=pval.onesided[1],sspval.twosided=pval.twosided[1],
              aggpval.onesided=aggpval.onesided,aggpval.twosided=aggpval.twosided,
              teststat=teststat,weights.nulldistr=weights.nulldistr,
              active.last=res.multisplit[[b.splits]]$active,sig.last=res.multisplit[[b.splits]]$sig,wi.last=res.multisplit[[b.splits]]$wi))
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

