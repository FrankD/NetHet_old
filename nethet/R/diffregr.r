##########################################################################################
# Differential Regression: Two-sample testing for high-dimensional regression
#-------------------------------------------------------------------------------
# * Intercepts are assumed to be zero in log-likelihood (mu1=mu2).
#   Therefore, center input data (y1,y2,x1,x2).
#
#
#
##########################################################################################


#####################
##Required Packages##
#####################
library(mvtnorm)
library(glmnet)
library(CompQuadForm)


#######################
##-----Screening-----##
#######################

##' Cross-validated Lasso screening (lambda.min-rule)
##'
##' 
##' @title Cross-validation lasso screening (lambda.min-rule)
##' @param x Predictor matrix
##' @param y Response vector
##' @return Active-set
##' @author n.stadler
##' @export
screen_cvmin.lasso <- function(x,y){
  fit.cv <- cv.glmnet(x,y)
  beta <- as.numeric(coef(fit.cv,s='lambda.min')[-1])
  p <- length(beta)
  n <- nrow(x)
  beta[-(order(abs(beta),decreasing=TRUE)[1:min(p,n)])] <- 0
  return(which(beta!=0))
}

##' Cross-validated Lasso screening (lambda.1se-rule)
##'
##' 
##' @title Cross-validated Lasso screening (lambda.1se-rule)
##' @param x Predictor matrix
##' @param y Response vector 
##' @return Active-set
##' @author n.stadler
##' @export
screen_cv1se.lasso <- function(x,y){
  fit.cv <- cv.glmnet(x,y)
  beta <- as.numeric(coef(fit.cv,s='lambda.1se')[-1])
  p <- length(beta)
  n <- nrow(x)
  d <- min(p,n)
  beta[-order(abs(beta),decreasing=TRUE)[1:d]] <- 0
  return(which(beta!=0))
}

##' Cross-validated Lasso screening and additional truncation.
##'
##' Computes Lasso coefficients (cross-validation optimal lambda). Truncates
##' smallest coefficients to zero, such that there are no more than n/k.trunc
##' non-zero coefficients.
##' 
##' @title Cross-validated Lasso screening and additional truncation.
##' @param x Predictor matrix.
##' @param y Response vector.
##' @param k.trunc Truncation constant="number of samples per predictor" (default=5).
##' @return Active-set.
##' @author n.stadler
##' @export
screen_cvtrunc.lasso <- function(x,y,k.trunc=5){
  n <- nrow(x)
  fit.cv <- cv.glmnet(x,y)
  beta <- as.numeric(coef(fit.cv,s='lambda.min')[-1])
  p <- length(beta)
  d <- min(floor(n/k.trunc),p,n)
  beta[-order(abs(beta),decreasing=TRUE)[1:d]] <- 0
  return(which(beta!=0))
}

##' Cross-validated Lasso screening and sqrt-truncation. 
##'
##' Computes Lasso coefficients (cross-validation optimal lambda). Truncates
##' smallest coefficients to zero, such that there are no more than sqrt(n)
##' non-zero coefficients.
##' 
##' @title Cross-validated Lasso screening and sqrt-truncation.
##' @param x Predictor matrix.
##' @param y Response vector. 
##' @return Active-set.
##' @author n.stadler
##' @export
screen_cvsqrt.lasso <- function(x,y){
  n <- nrow(x)
  fit.cv <- cv.glmnet(x,y)
  beta <- as.numeric(coef(fit.cv,s='lambda.min')[-1])
  p <- length(beta)
  d <- min(floor(sqrt(n)),p,n)
  beta[-order(abs(beta),decreasing=TRUE)[1:d]] <- 0
  return(which(beta!=0))
}

##' Cross-validated Lasso screening and upper bound on number of predictors
##'
##' Computes Lasso coefficients (cross-validation optimal lambda). Truncates
##' smalles coefficients to zero such that there are no more than no.predictors
##' non-zero coefficients
##' 
##' @title Cross-validated Lasso screening and upper bound on number of predictors.
##' @param x Predictor matrix.
##' @param y Response vector. 
##' @param no.predictors Upper bound on number of active predictors,
##' @return Active-set.
##' @author n.stadler
##' @export
screen_cvfix.lasso <- function(x,y,no.predictors=10){
  n <- nrow(x)
  fit.cv <- cv.glmnet(x,y)
  beta <- as.numeric(coef(fit.cv,s='lambda.min')[-1])
  p <- length(beta)
  d <- min(no.predictors,p,n)
  beta[-order(abs(beta),decreasing=TRUE)[1:d]] <- 0
  return(which(beta!=0))
}

##############################
##--------P-VALUES----------##
##############################

##' Log-likelihood ratio statistics for Differential Regression.
##'
##' 
##' @title Log-likelihood ratio statistics for Differential Regression.
##' @param y1 Response vector condition 1.
##' @param y2 Response vector condition 2.
##' @param y Pooled response vector.
##' @param xx1 Predictor matrix condition 1.
##' @param xx2 Predictor matrix condition 2.
##' @param xx Pooled predictor matrix
##' @param beta1 Regression coefficients condition 1.
##' @param beta2 Regression coefficients condition 2.
##' @param beta Pooled regression coefficients.
##' @return 2 times log-likelihood ratio statistics.
##' @author n.stadler
logratio.diffregr <- function(y1,y2,y,xx1,xx2,xx,beta1,beta2,beta){
  ##Compute 2*log-likelihood ratio
  ##Input:
  ##-y1,y2,y
  ##-xx1,xx2,xx (include only selected features)
  ##-MLE estimates beta1,beta2,beta

  n1 <- length(y1)
  n2 <- length(y2)
  n <- length(y)
  mu1<-mu2<-mu<-0
  sig1 <- (sum(y1^2)/n1);sig2 <- (sum(y2^2)/n2);sig <- (sum(y^2)/n)
  if(length(beta1)!=0){
    mu1<-xx1%*%beta1
    sig1 <- sum((y1-mu1)^2)/n1
  }
  if(length(beta2)!=0){
    mu2<-xx2%*%beta2
    sig2 <- sum((y2-mu2)^2)/n2
  }
  if(length(beta)!=0){
    mu<-xx%*%beta
    sig <- sum((y-mu)^2)/n
  }

  2*(sum(dnorm(y1,mean=mu1,sd=sqrt(sig1),log=TRUE))+sum(dnorm(y2,mean=mu2,sd=sqrt(sig2),log=TRUE))-sum(dnorm(y,mean=mu,sd=sqrt(sig),log=TRUE)))
}

##' Computation M matrix and eigenvalues
##'
##' 
##' @title Computation M matrix and eigenvalues
##' @param Sig no descr
##' @param act no descr
##' @param act1 no descr
##' @param act2 no descr
##' @return no descr
##' @author n.stadler
ww.mat.diffregr <- function(Sig,act,act1,act2){
  ##Compute W and Eval(W) (without simplification of W)
  ##
  ##Input:
  ##  Sig: E_0[s(X)s(X)´] (information matrix)
  ##  act: active variables for joint model
  ##  act1,act2: active variables for individual models
  bfg <- rbind(Sig[act1,act,drop=FALSE],Sig[act2,act,drop=FALSE])
  bf <- matrix(0,length(act1)+length(act2),length(act1)+length(act2))
  bf[1:length(act1),1:length(act1)] <- Sig[act1,act1]
  bf[length(act1)+(1:(length(act2))),length(act1)+(1:(length(act2)))] <- Sig[act2,act2]
  bg <- 2*Sig[act,act,drop=FALSE]
  mat <- rbind(cbind(diag(1,length(act1)+length(act2)),bfg%*%solve(bg)),cbind(-t(bfg)%*%solve(bf),diag(-1,length(act))))
  eval <- as.double(eigen(mat)$values)
  eval[abs(eval)<10^{-6}] <- 0
  return(list(ww.mat=mat,eval=eval))
}

##' Computation M matrix and eigenvalues
##'
##' 
##' @title Computation M matrix and eigenvalues
##' @param Sig no descr
##' @param act no descr
##' @param act1 no descr
##' @param act2 no descr
##' @return no descr
##' @author n.stadler
ww.mat2.diffregr <- function(Sig,act,act1,act2){
  ##Compute W and Eval(W) ('1st order' simplification of W)
  ##
  ##Input:
  ##  Sig: E_0[s(X)s(X)´]
  ##  act: active variables for joint model (including sigma^2)
  ##  act1,act2: active variables for individual models (including sigma_1^2, sigma_2^2)
  dimf <- length(act1)+length(act2)
  dimg <- length(act)
  bfg <- rbind(Sig[act1,act,drop=FALSE],Sig[act2,act,drop=FALSE])
  bgf <- t(bfg)
  bf <- matrix(0,length(act1)+length(act2),length(act1)+length(act2))
  bf[1:length(act1),1:length(act1)] <- Sig[act1,act1]
  bf[length(act1)+(1:(length(act2))),length(act1)+(1:(length(act2)))] <- Sig[act2,act2]
  bg <- 2*Sig[act,act]
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

##' Computation Q matrix
##'
##' 
##' @title Computation Q matrix
##' @param Sig no descr
##' @param a no descr
##' @param b no descr
##' @param s no descr
##' @return no descr
##' @author n.stadler
q.matrix.diffregr <- function(Sig,a,b,s){
  ##Compute Q: based on true information matrix (Sig)
  if(length(s)==0){
    return(Sig[a,b,drop=FALSE])
  }else{
    return(Sig[a,b,drop=FALSE]-Sig[a,s,drop=FALSE]%*%solve(Sig[s,s,drop=FALSE])%*%Sig[s,b,drop=FALSE])
  }
}

##' Computation eigenvalues
##'
##' 
##' @title Computation eigenvalues
##' @param Sig no descr
##' @param act no descr
##' @param act1 no descr
##' @param act2 no descr
##' @return no descr
##' @author n.stadler
my.ev2.diffregr <- function(Sig,act,act1,act2){
  ##Compute Eval (with 2nd order simplification of W)
  ##
  ##Input:
  ##  Sig: E_0[s(X)s(X)´]
  ##  act: active beta's for joint model (without sigma^2)
  ##  act1,act2: active beta's for individual models (without sigma_1^2, sigma_2^2)
  
  #dimension of models
  dimf1 <- length(act1)+1
  dimf2 <- length(act2)+1
  dimf <- dimf1+dimf2
  dimg <- length(act)+1
  #intersection of models
  ss <- intersect(act,intersect(act1,act2))
  if(length(ss)==0){warning('no intersection between models')}
  aa <- setdiff(act1,ss)
  bb <- setdiff(act2,ss)
  cc <- setdiff(act,ss)
  ev.aux <- ev.aux.complex<-numeric(0)
  if (dimf>=dimg){
    if (length(cc)!=0){
      qcc <- q.matrix.diffregr(Sig,cc,cc,ss)
      aux.mat <- matrix(0,length(cc),length(cc))
      if(length(aa)!=0){
        qac <- q.matrix.diffregr(Sig,aa,cc,ss)
        qaa <- q.matrix.diffregr(Sig,aa,aa,ss)
        aux.mat <- aux.mat+(t(qac)%*%solve(qaa)%*%qac)%*%solve(qcc)
      }
      if(length(bb)!=0){
        qbc <- q.matrix.diffregr(Sig,bb,cc,ss)
        qbb <- q.matrix.diffregr(Sig,bb,bb,ss)
        aux.mat <- aux.mat+(t(qbc)%*%solve(qbb)%*%qbc)%*%solve(qcc)
      }
      ev.aux.complex <- eigen(aux.mat)$values
      ev.aux <- as.double(ev.aux.complex)
      ev.aux <-  sqrt(1-ev.aux/2)
    }
    eval<-rep(1,dimf-dimg)
    eval <- c(eval,rep(0,2*(length(ss)+1)),ev.aux,-ev.aux)
  }# end if (dimf>=dimg){
  if (dimf<dimg){
    if (length(cc)!=0){
      qcc <- q.matrix.diffregr(Sig,cc,cc,ss)
      if(length(aa)!=0){
        qac <- q.matrix.diffregr(Sig,aa,cc,ss)
        qaa <- q.matrix.diffregr(Sig,aa,aa,ss)
        oaa <- qac%*%solve(qcc)%*%t(qac)
        aux.mat.aa<- oaa%*%solve(qaa)
        aux.mat <- aux.mat.aa
      }
      if(length(bb)!=0){
        qbc <- q.matrix.diffregr(Sig,bb,cc,ss)
        qbb <- q.matrix.diffregr(Sig,bb,bb,ss)
        obb <- qbc%*%solve(qcc)%*%t(qbc)
        aux.mat.bb<- obb%*%solve(qbb)
        aux.mat <- aux.mat.bb
      }
      if((length(aa)!=0)&(length(bb)!=0)){
        oab <- qac%*%solve(qcc)%*%t(qbc)
        aux.mat<- rbind(cbind(aux.mat.aa,oab%*%solve(qbb)),cbind(t(oab)%*%solve(qaa),aux.mat.bb))
      }
    }
    if ((length(aa)!=0)|(length(bb)!=0)){##if (length(aa)==0)&(length(bb)==0) 'MI included in MJ; therefore -X^2(dim(MJ)-dim(MI)) distributed'
      ev.aux.complex <- eigen(aux.mat)$values
      ev.aux <- as.double(ev.aux.complex)
      ev.aux <-  sqrt(1-ev.aux/2)
    }
    eval<-rep(-1,dimg-dimf)
    eval <- c(eval,rep(-1,(length(ss)+1)),rep(1,(length(ss)+1)),rep(0,2*(length(ss)+1)),ev.aux,-ev.aux)
  }
  return(list(eval=eval,ev.aux.complex=ev.aux.complex))
}

##' Computation beta matrix
##'
##' 
##' @title Computation beta matrix
##' @param ind1 no descr
##' @param ind2 no descr
##' @param beta1 no descr
##' @param beta2 no descr
##' @param beta no descr
##' @param sig1 no descr
##' @param sig2 no descr
##' @param sig no descr
##' @param Sig no descr
##' @return no descr
##' @author n.stadler
beta.mat.diffregr<-function(ind1,ind2,beta1,beta2,beta,sig1,sig2,sig,Sig){
  ##compute Beta-Matrix
  ##Input:
  ##-beta1,beta2
  ## beta (for expectation)
  ##-sig1,sig2
  ## sig (for expectation)
  ##-Sig (variance of predictors=var(x))

  beta.mat <- Sig*(sig+as.numeric(t(beta-beta1)%*%Sig%*%(beta-beta2)))/(sig1*sig2)
  return(beta.mat[ind1,ind2,drop=FALSE])
}

##' Computation Q matrix
##'
##' 
##' @title Computation Q matrix
##' @param beta.a no descr
##' @param beta.b no descr
##' @param beta no descr
##' @param sig.a no descr
##' @param sig.b no descr
##' @param sig no descr
##' @param Sig no descr
##' @param act.a no descr
##' @param act.b no descr
##' @param ss no descr
##' @return no descr
##' @author n.stadler
q.matrix.diffregr3 <- function(beta.a,beta.b,beta,sig.a,sig.b,sig,Sig,act.a,act.b,ss){

 ##Estimate Q
    b.ab<-beta.mat.diffregr(act.a,act.b,beta.a,beta.b,beta,sig.a,sig.b,sig,Sig)
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

##' Computation Q matrix
##'
##' 
##' @title Computation Q matrix
##' @param b.mat no descr
##' @param act.a no descr
##' @param act.b no descr
##' @param ss no descr
##' @return no descr
##' @author n.stadler
q.matrix.diffregr4 <- function(b.mat,act.a,act.b,ss){

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

##' Estimate weights 
##'
##' estimate W-matrix (using plug-in estimates of Beta-matrix); calculate eigenvalues(W-matrix)
##' 
##' @title Estimate weights
##' @param y1 no descr
##' @param y2 no descr
##' @param x1 no descr
##' @param x2 no descr
##' @param beta1 no descr
##' @param beta2 no descr
##' @param beta no descr
##' @param act1 no descr
##' @param act2 no descr
##' @param act no descr
##' @return no descr
##' @author n.stadler
est2.ww.mat.diffregr <- function(y1,y2,x1,x2,beta1,beta2,beta,act1,act2,act){

  n1 <- length(y1)
  n2 <- length(y2)
  y <- c(y1,y2)
  x <- rbind(x1,x2)
  n <- length(y)

  ##get mle's for sigma
  mu1<-mu2<-mu<-0
  sig1 <- (sum(y1^2)/n1);sig2 <- (sum(y2^2)/n2);sig <- (sum(y^2)/n)
  if(length(beta1)!=0){
    mu1<-x1[,act1,drop=FALSE]%*%beta1
    sig1 <- sum((y1-mu1)^2)/n1
  }
  if(length(beta2)!=0){
    mu2<-x2[,act2,drop=FALSE]%*%beta2
    sig2 <- sum((y2-mu2)^2)/n2
  }
  if(length(beta)!=0){
    mu<-x[,act,drop=FALSE]%*%beta
    sig <- sum((y-mu)^2)/n
  }

  ##expand beta's with zeros
  exp.beta <- exp.beta1 <- exp.beta2 <- rep(0,ncol(x1))
  exp.beta[act] <- beta
  exp.beta1[act1] <- beta1
  exp.beta2[act2] <- beta2

  ##compute beta.mat
  bact.pop1<- beta.mat.diffregr(act,act,exp.beta,exp.beta,exp.beta1,sig,sig,sig1,var(x1))
  bact.pop2<- beta.mat.diffregr(act,act,exp.beta,exp.beta,exp.beta2,sig,sig,sig2,var(x2))
  bact1<- beta.mat.diffregr(act1,act1,exp.beta1,exp.beta1,exp.beta1,sig1,sig1,sig1,var(x1))
  bact2<- beta.mat.diffregr(act2,act2,exp.beta2,exp.beta2,exp.beta2,sig2,sig2,sig2,var(x2))
  bact1.act<- beta.mat.diffregr(act1,act,exp.beta1,exp.beta,exp.beta1,sig1,sig,sig1,var(x1))
  bact2.act<- beta.mat.diffregr(act2,act,exp.beta2,exp.beta,exp.beta2,sig2,sig,sig2,var(x2))
  
  bfg <- rbind(bact1.act,bact2.act)
  bgf <- t(bfg)
  bf <- matrix(0,length(act1)+length(act2),length(act1)+length(act2))
  bf[1:length(act1),1:length(act1)] <- bact1
  bf[length(act1)+(1:length(act2)),length(act1)+(1:length(act2))] <- bact2
  bg <- bact.pop1+bact.pop2
  mat <- rbind(cbind(diag(1,length(act1)+length(act2)),bfg%*%solve(bg)),cbind(-t(bfg)%*%solve(bf),diag(-1,length(act))))
  
  eval.complex<-eigen(mat)$values
  eval <- as.double(eval.complex)
  
  return(list(eval=c(1,0,0,eval),eval.complex=eval.complex))## eigenvalues 1,0,0 correspond to sig1, sig2 and sig12
}

##' Estimate weights
##'
##' 
##' @title Estimate weights
##' @param y1 no descr
##' @param y2 no descr
##' @param x1 no descr
##' @param x2 no descr
##' @param beta1 no descr
##' @param beta2 no descr
##' @param beta no descr
##' @param act1 no descr
##' @param act2 no descr
##' @param act no descr
##' @return no descr
##' @author n.stadler
est2.ww.mat2.diffregr <- function(y1,y2,x1,x2,beta1,beta2,beta,act1,act2,act){
  ##Estimate W and Eval(W) ('1st order' simplification of W)
  ##
  ##Input:
  ##  data: y1,y2,x1,x2
  ##  beta: estimate for joint model
  ##  beta1,beta2: estimates for individual model
  ##  act: active beta's for joint model
  ##  act1,act2: active beta's for individual models
  ##
  ##8!8! does not include weights from estimating error variances... (TO DO)

  n1 <- length(y1)
  n2 <- length(y2)
  y <- c(y1,y2)
  x <- rbind(x1,x2)
  n <- length(y)

  ##get mle's for sigma
  mu1<-mu2<-mu<-0
  sig1 <- (sum(y1^2)/n1);sig2 <- (sum(y2^2)/n2);sig <- (sum(y^2)/n)
  if(length(beta1)!=0){
    mu1<-x1[,act1,drop=FALSE]%*%beta1
    sig1 <- sum((y1-mu1)^2)/n1
  }
  if(length(beta2)!=0){
    mu2<-x2[,act2,drop=FALSE]%*%beta2
    sig2 <- sum((y2-mu2)^2)/n2
  }
  if(length(beta)!=0){
    mu<-x[,act,drop=FALSE]%*%beta
    sig <- sum((y-mu)^2)/n
  }

  ##expand beta's with zeros
  exp.beta <- exp.beta1 <- exp.beta2 <- rep(0,ncol(x1))
  exp.beta[act] <- beta
  exp.beta1[act1] <- beta1
  exp.beta2[act2] <- beta2

  ##dimension of models
  dimf1 <- length(act1)
  dimf2 <- length(act2)
  dimf <- dimf1+dimf2
  dimg <- length(act)

  ##compute beta.mat
  bact.pop1<- beta.mat.diffregr(act,act,exp.beta,exp.beta,exp.beta1,sig,sig,sig1,var(x1))
  bact.pop2<- beta.mat.diffregr(act,act,exp.beta,exp.beta,exp.beta2,sig,sig,sig2,var(x2))
  bact1<- beta.mat.diffregr(act1,act1,exp.beta1,exp.beta1,exp.beta1,sig1,sig1,sig1,var(x1))
  bact2<- beta.mat.diffregr(act2,act2,exp.beta2,exp.beta2,exp.beta2,sig2,sig2,sig2,var(x2))
  bact1.act<- beta.mat.diffregr(act1,act,exp.beta1,exp.beta,exp.beta1,sig1,sig,sig1,var(x1))
  bact2.act<- beta.mat.diffregr(act2,act,exp.beta2,exp.beta,exp.beta2,sig2,sig,sig2,var(x2))
  
  bfg <- rbind(bact1.act,bact2.act)
  bgf <- t(bfg)
  bf <- matrix(0,dimf,dimf)
  bf[1:dimf1,1:dimf1] <- bact1
  bf[dimf1+(1:dimf2),dimf1+(1:dimf2)] <- bact2
  bg <- bact.pop1+bact.pop2
  
  if (dimf>=dimg){
    mat <- bgf%*%solve(bf)%*%bfg%*%solve(bg)
    eval<-rep(1,dimf-dimg)
  }
  if (dimf<dimg){
    mat <- bfg%*%solve(bg)%*%bgf%*%solve(bf)
    eval<-rep(-1,dimg-dimf)
  }
  eval.mu.complex<-eigen(mat)$values
  eval.mu <- as.double(eval.mu.complex)
  eval <- c(eval,sqrt(1-eval.mu),-sqrt(1-eval.mu))
  
  return(list(eval=eval,eval.mu.complex=eval.mu.complex))
}

##' Compute weights of sum-w-chi2 (2nd order simplification)
##'
##' 
##' *expansion of W in two directions ("dimf>dimg direction" & "dimf>dimg direction") 
##' *simplified computation of weights is obtained by assuming H0 and that X_u~X_v holds
##' @param y1 no descr
##' @param y2 no descr
##' @param x1 no descr
##' @param x2 no descr
##' @param beta1 no descr
##' @param beta2 no descr
##' @param beta no descr
##' @param act1 no descr
##' @param act2 no descr
##' @param act no descr
##' @return no descr
##' @author n.stadler
est2.my.ev2.diffregr <- function(y1,y2,x1,x2,beta1,beta2,beta,act1,act2,act){

  ##Estimate Evals ('2nd order' simplification of W)
  ##
  ##Input:
  ##  data: y1,y2,x1,x2
  ##  beta: estimate for joint model
  ##  beta1,beta2: estimates for individual model
  ##  act: active beta's for joint model
  ##  act1,act2: active beta's for individual models

  n1 <- length(y1)
  n2 <- length(y2)
  y <- c(y1,y2)
  x <- rbind(x1,x2)
  n <- length(y)

  ##get mle's for sigma
  mu1<-mu2<-mu<-0
  sig1 <- (sum(y1^2)/n1);sig2 <- (sum(y2^2)/n2);sig <- (sum(y^2)/n)
  if(length(beta1)!=0){
    mu1<-x1[,act1,drop=FALSE]%*%beta1
    sig1 <- sum((y1-mu1)^2)/n1
  }
  if(length(beta2)!=0){
    mu2<-x2[,act2,drop=FALSE]%*%beta2
    sig2 <- sum((y2-mu2)^2)/n2
  }
  if(length(beta)!=0){
    mu<-x[,act,drop=FALSE]%*%beta
    sig <- sum((y-mu)^2)/n
  }

  ##expand beta's with zeros
  exp.beta <- exp.beta1 <- exp.beta2 <- rep(0,ncol(x1))
  exp.beta[act] <- beta
  exp.beta1[act1] <- beta1
  exp.beta2[act2] <- beta2

  ##dimension of models
  dimf1 <- length(act1)+1
  dimf2 <- length(act2)+1
  dimf <- dimf1+dimf2
  dimg <- length(act)+1
  
  ##intersection of models
  ss <- intersect(act,intersect(act1,act2))
  if(length(ss)==0){cat('warning! no intersection between models','\n')}
  aa <- setdiff(act1,ss)
  bb <- setdiff(act2,ss)
  cc <- setdiff(act,ss)

  ev.aux <- ev.aux.complex<-numeric(0)
  if (dimf>=dimg){
    if (length(cc)!=0){
      qcc <- q.matrix.diffregr3(exp.beta,exp.beta,exp.beta1,sig,sig,sig1,var(x1),act,act,ss)+q.matrix.diffregr3(exp.beta,exp.beta,exp.beta2,sig,sig,sig2,var(x2),act,act,ss)
      aux.mat <- matrix(0,length(cc),length(cc))
      if(length(aa)!=0){
        qac <- q.matrix.diffregr3(exp.beta1,exp.beta,exp.beta1,sig1,sig,sig1,var(x1),act1,act,ss)
        qaa <- q.matrix.diffregr3(exp.beta1,exp.beta1,exp.beta1,sig1,sig1,sig1,var(x1),act1,act1,ss)
        aux.mat <- aux.mat+(t(qac)%*%solve(qaa)%*%qac)%*%solve(qcc)
      }
      if(length(bb)!=0){
        qbc <- q.matrix.diffregr3(exp.beta2,exp.beta,exp.beta2,sig2,sig,sig2,var(x2),act2,act,ss)
        qbb <- q.matrix.diffregr3(exp.beta2,exp.beta2,exp.beta2,sig2,sig2,sig2,var(x2),act2,act2,ss)
        aux.mat <- aux.mat+(t(qbc)%*%solve(qbb)%*%qbc)%*%solve(qcc)
      }
      ev.aux.complex <- eigen(aux.mat)$values
      ev.aux <- as.double(ev.aux.complex)
      ev.aux <-  sqrt(1-ev.aux)
    }
    eval<-rep(1,dimf-dimg)
    eval <- c(eval,rep(0,2*(length(ss)+1)),ev.aux,-ev.aux)
  }# end if (dimf>=dimg){
  if (dimf<dimg){
    if (length(cc)!=0){
      qcc <- q.matrix.diffregr3(exp.beta,exp.beta,exp.beta1,sig,sig,sig1,var(x1),act,act,ss)+q.matrix.diffregr3(exp.beta,exp.beta,exp.beta2,sig,sig,sig2,var(x2),act,act,ss)
      aux.mat <- matrix(0,length(cc),length(cc))
      if(length(aa)!=0){
        qac <- q.matrix.diffregr3(exp.beta1,exp.beta,exp.beta1,sig1,sig,sig1,var(x1),act1,act,ss)
        qaa <- q.matrix.diffregr3(exp.beta1,exp.beta1,exp.beta1,sig1,sig1,sig1,var(x1),act1,act1,ss)
        oaa <- qac%*%solve(qcc)%*%t(qac)
        aux.mat.aa<- oaa%*%solve(qaa)
        aux.mat <- aux.mat.aa
      }
      if(length(bb)!=0){
        qbc <- q.matrix.diffregr3(exp.beta2,exp.beta,exp.beta2,sig2,sig,sig2,var(x2),act2,act,ss)
        qbb <- q.matrix.diffregr3(exp.beta2,exp.beta2,exp.beta2,sig2,sig2,sig2,var(x2),act2,act2,ss)
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
    eval <- c(eval,rep(-1,(length(ss)+1)),rep(1,(length(ss)+1)),rep(0,2*(length(ss)+1)),ev.aux,-ev.aux)
  }# end if (dimf<dimg){
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
##' @param y1 Response vector sample 1.
##' @param y2 Response vector sample 2.
##' @param x1 Predictor matrix sample 1.
##' @param x2 Predictor matrix sample 2.
##' @param beta1 MLE (regression coefficients) sample 1.
##' @param beta2 MLE (regression coefficients) sample 2.
##' @param beta Pooled MLE (regression coefficients).
##' @param act1 Active-set sample 1
##' @param act2 Active-set sample 2
##' @param act Pooled active-set
##' @return Eigenvalues of M, respectively the weights.
##' @author n.stadler
est2.my.ev3.diffregr <- function(y1,y2,x1,x2,beta1,beta2,beta,act1,act2,act){

  show.warn <- FALSE

  n1 <- length(y1)
  n2 <- length(y2)
  y <- c(y1,y2)
  x <- rbind(x1,x2)
  n <- length(y)

  ##get mle's for sigma
  mu1<-mu2<-mu<-0
  sig1 <- (sum(y1^2)/n1);sig2 <- (sum(y2^2)/n2);sig <- (sum(y^2)/n)
  if(length(beta1)!=0){
    mu1<-x1[,act1,drop=FALSE]%*%beta1
    sig1 <- sum((y1-mu1)^2)/n1
  }
  if(length(beta2)!=0){
    mu2<-x2[,act2,drop=FALSE]%*%beta2
    sig2 <- sum((y2-mu2)^2)/n2
  }
  if(length(beta)!=0){
    mu<-x[,act,drop=FALSE]%*%beta
    sig <- sum((y-mu)^2)/n
  }

  ##expand beta's with zeros
  exp.beta <- exp.beta1 <- exp.beta2 <- rep(0,ncol(x1))
  exp.beta[act] <- beta
  exp.beta1[act1] <- beta1
  exp.beta2[act2] <- beta2

  ##dimension of models
  dimf1 <- length(act1)+1
  dimf2 <- length(act2)+1
  dimf <- dimf1+dimf2
  dimg <- length(act)+1
  
  ##intersection of models
  ss <- intersect(act,intersect(act1,act2))
  if((length(ss)==0)&show.warn){cat('warning! no intersection between models','\n')}
  aa <- setdiff(act1,ss)
  bb <- setdiff(act2,ss)
  cc <- setdiff(act,ss)

  ev.aux <- ev.aux.complex<-numeric(0)
  no.zero.ev.aux <- 0
  #if (dimf>=dimg){
  if (length(cc)!=0){
    qcc1 <- q.matrix.diffregr3(exp.beta,exp.beta,exp.beta1,sig,sig,sig1,var(x1),act,act,ss)
    qcc2 <- q.matrix.diffregr3(exp.beta,exp.beta,exp.beta2,sig,sig,sig2,var(x2),act,act,ss)
    bmat <- beta.mat.diffregr(act,act,exp.beta,exp.beta,exp.beta1,sig,sig,sig1,var(x1))+beta.mat.diffregr(act,act,exp.beta,exp.beta,exp.beta2,sig,sig,sig2,var(x2))
    qcc12 <- q.matrix.diffregr4(bmat,act,act,ss)
    aux.mat <- diag(1,length(cc))-qcc1%*%solve(qcc12)-qcc2%*%solve(qcc12)
    if(length(aa)!=0){
      qac <- q.matrix.diffregr3(exp.beta1,exp.beta,exp.beta1,sig1,sig,sig1,var(x1),act1,act,ss)
      qaa <- q.matrix.diffregr3(exp.beta1,exp.beta1,exp.beta1,sig1,sig1,sig1,var(x1),act1,act1,ss)
      aux.mat <- aux.mat+(t(qac)%*%solve(qaa)%*%qac)%*%solve(qcc12)
    }
    if(length(bb)!=0){
      qbc <- q.matrix.diffregr3(exp.beta2,exp.beta,exp.beta2,sig2,sig,sig2,var(x2),act2,act,ss)
      qbb <- q.matrix.diffregr3(exp.beta2,exp.beta2,exp.beta2,sig2,sig2,sig2,var(x2),act2,act2,ss)
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
  eval <- c(eval,rep(0,2*(length(ss)+1)),ev.aux,-ev.aux)
  #}# end if (dimf>=dimg){
  
  return(list(eval=eval,ev.aux.complex=ev.aux.complex))
}

##' Computation "split-asym"/"split-perm" p-values.
##'
##' 
##' @title Computation "split-asym" p-values.
##' @param y1 Response vector condition 1.
##' @param y2 Response vector condition 2.
##' @param x1 Predictor matrix condition 1.
##' @param x2 Predictor matrix condition 2.
##' @param beta1 Regression coefficients condition 1.
##' @param beta2 Regression coefficients condition 2.
##' @param beta Pooled regression coefficients.
##' @param act1 Active-set condition 1.
##' @param act2 Active-set condition 2.
##' @param act Pooled active-set.
##' @param compute.evals Method for computation of weights.
##' @param method.compquadform Method to compute distribution function of w-sum-of-chi2.
##' @param acc See ?davies.
##' @param epsabs See ?imhof.
##' @param epsrel See ?imhof.
##' @param show.warn Show warnings?
##' @param n.perm Number of permutations.
##' @return P-value, test statistic, estimated weights.
##' @author n.stadler
diffregr_pval <- function(y1,y2,x1,x2,beta1,beta2,beta,act1,act2,act,compute.evals,method.compquadform,acc,epsabs,epsrel,show.warn,n.perm){

  if(is.null(n.perm)){
    ##########################
    ##compute test-statistic##
    ##########################
    teststat <- logratio.diffregr(y1,y2,c(y1,y2),x1[,act1,drop=FALSE],x2[,act2,drop=FALSE],rbind(x1,x2)[,act,drop=FALSE],beta1,beta2,beta)
    #################################
    ##compute weights of sum-w-chi2##
    #################################
    weights.nulldistr <- eval(as.name(compute.evals))(y1,y2,x1,x2,beta1,beta2,beta,act1,act2,act)$eval
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
  }else{
    
    return(perm.diffregr_pval(y1,y2,x1,x2,act1,act2,act,n.perm))
  }
}

##' Auxiliary function for computation of "split-perm" p-value.
##'
##' 
##' @title Auxiliary function for computation of "split-perm" p-value.
##' @param y1 Response vector condition 1.
##' @param y2 Response vector condition 2.
##' @param y12 Pooled response vector.
##' @param x1 Predictor matrix condition 1.
##' @param x2 Predictor matrix condition 2.
##' @param x12 Pooled predictor matrix
##' @return Test statistic (log-likelihood-ratio statistic).
##' @author n.stadler
perm.diffregr_teststat <- function(y1,y2,y12,x1,x2,x12){
  p1 <- ncol(x1)
  p2 <- ncol(x2)
  p12 <- ncol(x12)
  
  if(p1!=0){
    beta1 <- as.numeric(coef(lm(y1~x1-1)))
  }else{
   beta1<-numeric(0)
  }
  if(p2!=0){
    beta2 <- as.numeric(coef(lm(y2~x2-1)))
  }else{
   beta2<-numeric(0)
  }
  if(p12!=0){
    beta <- as.numeric(coef(lm(y12~x12-1)))
  }else{
   beta<-numeric(0)
  }
  return(logratio.diffregr(y1,y2,y12,x1,x2,x12,beta1,beta2,beta))
}

##' Computation "split-perm" p-value.
##'
##' 
##' @title Computation "split-perm" p-value.
##' @param y1 Response vector condition 1.
##' @param y2 Response vector condition 2.
##' @param x1 Predictor matrix condition 1.
##' @param x2 Predictor matrix condition 2.
##' @param act1 Active-set condition 1.
##' @param act2 Active-set condition 2.
##' @param act Pooled active-set.
##' @param n.perm Number of permutations.
##' @return Permutation based p-value.
##' @author n.stadler
perm.diffregr_pval <- function(y1,y2,x1,x2,act1,act2,act,n.perm){
  n1 <- nrow(x1);n2 <- nrow(x2)
  x12 <- rbind(x1,x2)
  x12.act1 <- x12[,act1,drop=FALSE]
  x12.act2 <- x12[,act2,drop=FALSE]
  x12.act <- x12[,act,drop=FALSE]
  y12 <- c(y1,y2)

  tobs <- perm.diffregr_teststat(y1,y2,y12,x12.act1[1:n1,,drop=FALSE],x12.act2[(n1+1):(n1+n2),,drop=FALSE],x12.act)

  tperm <- sapply(1:n.perm,
                  function(i){
                    my.perm <- sample(c(rep(1,n1),rep(2,n2)))
                    x1.p <- x12.act1[my.perm==1,,drop=FALSE]
                    x2.p <- x12.act2[my.perm==2,,drop=FALSE]
                    y1.p <- y12[my.perm==1]
                    y2.p <- y12[my.perm==2]
                    perm.diffregr_teststat(y1.p,y2.p,y12,x1.p,x2.p,x12.act)
                  }
                  )
  pval <- (1+sum(tperm>=tobs))/n.perm
  return(list(pval.onesided=pval,pval.twosided=pval,weights.nulldistr=NULL,teststat=tobs))
}
              
##' Differential Regression (single-split version).
##'
##' Intercepts in regression models are assumed to be zero (mu1=mu2=0).
##' You might need to center the input data prior to running
##' Differential Regression.
##' 
##' @title Differential Regression (single-split version).
##' @param y1 Response vector condition 1.
##' @param y2 Response vector condition 2.
##' @param x1 Predictor matrix condition 1.
##' @param x2 Predictor matrix condition 2.
##' @param split1 Samples condition 1 used in screening-step.
##' @param split2 Samples condition 2 used in screening-step.
##' @param screen.meth Screening method (default='screen_cvtrunc.lasso').
##' @param compute.evals Method to estimate the weights in the weighted-sum-of-chi2s distribution.
##'                      The default and (currently) the only available option 
##'                      is the method 'est2.my.ev3.diffregr'.
##' @param method.compquadform Algorithm for computing distribution function
##'                            of weighted-sum-of-chi2 (default='imhof').
##' @param acc See ?davies (default=1e-4).
##' @param epsabs See ?imhof (default=1e-10).
##' @param epsrel See ?imhof (default=1e-10).
##' @param show.warn Show warnings (default=FALSE)? 
##' @param n.perm Number of permutation for "split-perm" p-value (default=NULL).
##' @param ... Other arguments specific to screen.meth.
##' @return List consisting of
##' \item{pval.onesided}{"One-sided" p-value.}
##' \item{pval.twosided}{"Two-sided" p-value. Ignore all "*.twosided results.}
##' \item{teststat}{2 times Log-likelihood-ratio statistics}
##' \item{weights.nulldistr}{Estimated weights of weighted-sum-of-chi2s.}
##' \item{active}{List of active-sets obtained in screening step.}
##' \item{beta}{Regression coefficients (MLE) obtaind in cleaning-step.}
##' @author n.stadler
##' @export
diffregr_singlesplit<- function(y1,y2,x1,x2,split1,split2,screen.meth='screen_cvtrunc.lasso',
                                compute.evals='est2.my.ev3.diffregr',method.compquadform='imhof',acc=1e-04,
                                epsabs=1e-10,epsrel=1e-10,
                                show.warn=FALSE,n.perm=NULL,...){
  
  n1 <- nrow(x1)
  n2 <- nrow(x2)
 
  est.beta <- active <- list()#save est.beta & active variables

  ###############
  ##Joint Model##
  ###############
  xx.train <- rbind(x1[split1,],x2[split2,])
  yy.train <- c(y1[split1],y2[split2])
  xx.valid <- rbind(x1[-split1,],x2[-split2,])
  yy.valid <- c(y1[-split1],y2[-split2])
  active[['modJ']] <- eval(as.name(screen.meth))(xx.train,yy.train,...)
  if(length(active[['modJ']])!=0){
    est.beta[['modJ']] <- as.numeric(coef(lm(yy.valid~xx.valid[,active[['modJ']]]-1)))
  }else{
    est.beta[['modJ']]<-numeric(0)
  }
  ####################  
  ##Individual Model##
  ####################
  for (j in c('1','2')){
    split.train <- eval(as.name(paste('split',j,sep='')))
    xx.train <- eval(as.name(paste('x',j,sep='')))[split.train,]
    yy.train <- eval(as.name(paste('y',j,sep='')))[split.train]
    xx.valid <- eval(as.name(paste('x',j,sep='')))[-split.train,]
    yy.valid <- eval(as.name(paste('y',j,sep='')))[-split.train]
    active[[paste('modIpop',j,sep='')]] <- eval(as.name(screen.meth))(xx.train,yy.train,...)

    if(length(active[[paste('modIpop',j,sep='')]])!=0){
      est.beta[[paste('modIpop',j,sep='')]] <- as.numeric(coef(lm(yy.valid~xx.valid[,active[[paste('modIpop',j,sep='')]]]-1)))
    }else{
      est.beta[[paste('modIpop',j,sep='')]]<-numeric(0)
    }
  }
  ###############
  ##Some Checks##
  ###############
  l.act <- lapply(active,length)
  n1.valid <- nrow(x1[-split1,])
  n2.valid <- nrow(x2[-split2,])
  if (any(l.act==0)){ if(show.warn){cat('warning: at least one active-set is empty','\n')}}
  if (all(l.act<c(n1.valid+n2.valid,n1.valid,n2.valid))){
    res.pval <- diffregr_pval(y1=y1[-split1],y2=y2[-split2],
                              x1=x1[-split1,,drop=FALSE],x2=x2[-split2,,drop=FALSE],
                              beta1=est.beta[['modIpop1']],beta2=est.beta[['modIpop2']],beta=est.beta[['modJ']],
                              active[['modIpop1']],active[['modIpop2']],active[['modJ']],
                              compute.evals,method.compquadform,acc,epsabs,epsrel,show.warn,n.perm)
  }else{
    if(show.warn){cat('warning: dim(model) > n: pval=NA','\n')}
    res.pval <- list(pval.onesided=NA,pval.twosided=NA,weights.nulldistr=NA,teststat=NA)
  }
  
  return(list(pval.onesided=res.pval$pval.onesided,pval.twosided=res.pval$pval.twosided,
              teststat=res.pval$teststat,weights.nulldistr=res.pval$weights.nulldistr,
              active=active,beta=est.beta))
}

##' Differential Regression (multi-split version).
##'
##' Intercepts in regression models are assumed to be zero (mu1=mu2=0).
##' You might need to center the input data prior to running
##' Differential Regression.
##'
##' 
##' @title Differential Regression (multi-split version).
##' @param y1 Response vector condition 1.
##' @param y2 Response vector condition 2.
##' @param x1 Predictor matrix condition 1.
##' @param x2 Predictor matrix condition 2.
##' @param b.splits Number of splits (default=50).
##' @param frac.split Fraction train-data (screening) / test-data (cleaning) (default=0.5).
##' @param screen.meth Screening method (default='screen_cvtrunc.lasso').
##' @param gamma.min Tuning parameter in p-value aggregation of Meinshausen et al (2009) (default=0.05).
##' @param compute.evals Method to estimate the weights in the weighted-sum-of-chi2s distribution.
##'                      The default and (currently) the only available option 
##'                      is the method 'est2.my.ev3.diffregr'.
##' @param method.compquadform Algorithm for computing distribution function
##'                            of weighted-sum-of-chi2 (default='imhof').
##' @param acc See ?davies (default=1e-4).
##' @param epsabs See ?imhof (default=1e-10).
##' @param epsrel See ?imhof (default=1e-10).
##' @param show.warn Show warnings (default=FALSE)?
##' @param n.perm Number of permutation for "split-perm" p-value. Default=NULL, which means
##'               that the asymptotic approximation is used.
##' @param ... Other arguments specific to screen.meth.
##' @return List consisting of
##' \item{ms.pval}{p-values for all b.splits}
##' \item{ss.pval}{single-split p-value}
##' \item{medagg.pval}{median aggregated p-value}
##' \item{meinshagg.pval}{meinshausen aggregated p-value (meinshausen et al 2009)}
##' \item{teststat}{test statistics for b.splits}
##' \item{weights.nulldistr}{estimated weights}
##' \item{active.last}{active-sets obtained in last screening-step}
##' \item{beta.last}{constrained mle (regression coefficients) obtained in last cleaning-step}
##' @author n.stadler
##' @export
##' @example ../diffregr_ex.R
diffregr_multisplit<- function(y1,y2,x1,x2,b.splits=50,frac.split=1/2,screen.meth='screen_cvtrunc.lasso',
                               gamma.min=0.05,compute.evals='est2.my.ev3.diffregr',
                               method.compquadform='imhof',acc=1e-04,epsabs=1e-10,epsrel=1e-10,
                               show.warn=FALSE,n.perm=NULL,...){

  n1 <- nrow(x1)
  n2 <- nrow(x2)
  
  res.multisplit <- lapply(seq(b.splits),
                          function(i){
                            split1 <- sample(1:n1,floor((n1-1)*frac.split),replace=FALSE)
                            split2 <- sample(1:n2,floor((n2-1)*frac.split),replace=FALSE)
                            res.singlesplit <- diffregr_singlesplit(y1,y2,x1,x2,split1,split2,screen.meth,
                                                                    compute.evals,method.compquadform,
                                                                    acc,epsabs,epsrel,show.warn,n.perm,...)                      
                          })
  pval.onesided <- sapply(res.multisplit,function(x){x[['pval.onesided']]},simplify='array')
  pval.twosided <- sapply(res.multisplit,function(x){x[['pval.twosided']]},simplify='array')
  teststat <- sapply(res.multisplit,function(x){x[['teststat']]},simplify='array')
  weights.nulldistr <- sapply(res.multisplit,function(x){x[['weights.nulldistr']]},simplify='array')
  aggpval.onesided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.onesided[!is.na(pval.onesided)])$objective)
  aggpval.twosided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.twosided[!is.na(pval.twosided)])$objective)

  result <- list(ms.pval=pval.onesided,
                 ss.pval=pval.onesided[1],
                 medagg.pval=median(pval.onesided,na.rm=TRUE),
                 meinshagg.pval=aggpval.onesided,
                 teststat=teststat,weights.nulldistr=weights.nulldistr,
                 active.last=res.multisplit[[b.splits]]$active,
                 beta.last=res.multisplit[[b.splits]]$beta)
  class(result) <- 'diffregr'
  return(result)
}
    

##' Old single-split function for diffregr
##'
##' 
##' @title old single-split function for diffregr
##' @param y1 no descr
##' @param y2 no descr
##' @param x1 no descr
##' @param x2 no descr
##' @param n.screen.pop1 no descr
##' @param n.screen.pop2 no descr
##' @param screen.meth no descr
##' @param compute.evals no descr
##' @return no descr
##' @author n.stadler
twosample_single_regr <- function(y1,y2,x1,x2,n.screen.pop1=100,n.screen.pop2=100,screen.meth='screen_cvmin.lasso',compute.evals='est2.my.ev3.diffregr'){

  ##Single-Split Pvalues
  ##
  ##Input:
  ##
  ## -Data:y1,y2,x1,x2
  ## -n.screen.pop1, n.screen.pop2: no samples for screening
  ## -screen.meth={lasso.cvmin,lasso.cv1se}
  ## -compute.evals: 'est2.my.ev2.diffregr'

  n1 <- nrow(x1)
  n2 <- nrow(x2)

  ##split data
  split1 <- sample(1:n1,n.screen.pop1,replace=FALSE)
  split2 <- sample(1:n2,n.screen.pop2,replace=FALSE)

  est.beta <- active <- list()#save est.beta & active variables

  ##model joint
  xx.train <- rbind(x1[split1,],x2[split2,])
  yy.train <- c(y1[split1],y2[split2])
  xx.valid <- rbind(x1[-split1,],x2[-split2,])
  yy.valid <- c(y1[-split1],y2[-split2])
  active[['modJ']] <- eval(as.name(screen.meth))(xx.train,yy.train)
  if(length(active[['modJ']])!=0){
    est.beta[['modJ']] <- as.numeric(coef(lm(yy.valid~xx.valid[,active[['modJ']]]-1)))
  }else{
    est.beta[['modJ']]<-numeric(0)
  }
    
  ##model individual
  for (j in c('1','2')){
    split.train <- eval(as.name(paste('split',j,sep='')))
    xx.train <- eval(as.name(paste('x',j,sep='')))[split.train,]
    yy.train <- eval(as.name(paste('y',j,sep='')))[split.train]
    xx.valid <- eval(as.name(paste('x',j,sep='')))[-split.train,]
    yy.valid <- eval(as.name(paste('y',j,sep='')))[-split.train]
    active[[paste('modIpop',j,sep='')]] <- eval(as.name(screen.meth))(xx.train,yy.train)

    if(length(active[[paste('modIpop',j,sep='')]])!=0){
      est.beta[[paste('modIpop',j,sep='')]] <- as.numeric(coef(lm(yy.valid~xx.valid[,active[[paste('modIpop',j,sep='')]]]-1)))
    }else{
      est.beta[[paste('modIpop',j,sep='')]]<-numeric(0)
    }
  }
  
  l.act <- lapply(active,length)
  n1.valid <- nrow(x1[-split1,])
  n2.valid <- nrow(x2[-split2,])
  if (any(l.act==0)){ warning('at least one active-set is empty')}
    
  if (all(l.act<c(n1.valid+n2.valid,n1.valid,n2.valid))){
      
    teststat <- logratio.diffregr(y1[-split1],y2[-split2],c(y1[-split1],y2[-split2]),
                         x1[-split1,active[['modIpop1']],drop=FALSE],x2[-split2,active[['modIpop2']],drop=FALSE],
                         rbind(x1[-split1,active[['modJ']],drop=FALSE],x2[-split2,active[['modJ']],drop=FALSE]),
                         est.beta[['modIpop1']],est.beta[['modIpop2']],est.beta[['modJ']])
      
    ev.nulldistr <- eval(as.name(compute.evals))(y1[-split1],y2[-split2],x1[-split1,],x2[-split2,],
                                                 est.beta[['modIpop1']],est.beta[['modIpop2']],est.beta[['modJ']],
                                                 active[['modIpop1']],active[['modIpop2']],active[['modJ']])$eval
    if (any(is.na(ev.nulldistr))){
      warning('Eigenval is NA: pval=NA')
    }else{
      sspval.onesided <- davies(teststat,lambda=ev.nulldistr)$Qq
      sspval.twosided <- 2*min(sspval.onesided,1-sspval.onesided)
    }
  }else{warning('dim(model) > n-1: pval=NA')}

    
  return(list(sspval.onesided=sspval.onesided,
              sspval.twosided=sspval.twosided,
              LR=teststat,active=active,beta=est.beta))
}


##' Plotting function for object of class 'diffregr' 
##'
##' 
##' @title Plotting function for object of class 'diffregr' 
##' @param x object of class 'diffregr'
##' @return Histogram over multi-split p-values.
##' @author nicolas
##' @export
plot.diffregr <- function(x,...){
        hh <- hist(x$ms.pval,
                   main='histogram single-split p-values',xlab='p-values',ylab='frequency',...)
        abline(v=x$medagg.pval,lty=2,col='red')
        abline(v=x$meinshagg.pval,lty=2,col='green')
        legend(x=min(hh$mids),y=max(hh$counts),lty=c(2,2),col=c('red','green'),legend=c('median aggregated','meinshausen aggregated'))
   
}

##' Summary function for object of class 'diffregr'
##'
##' 
##' @title Summary function for object of class 'diffregr'
##' @param x object of class 'diffregr'
##' @return aggregated p-values
##' @author nicolas
summary.diffregr <- function(x){
    out <- data.frame(medagg.pval=x$medagg.pval,meinshagg.pval=x$meinshagg.pval)
    rownames(out) <- 'aggregated p-values'
    return(out)
    print(out)
}
