###TwoSampleTest for HighDim Regression: assumes that mu1=mu2=mu0=0, therefore center input data (y1,y2,x1,x2)
###Date: 11/12/2012
###Changes: - estimate sig1,sig2,sig (see logratio)
###         - new computation of Evals
###         - some minor changes
###         - include sig1, sig2 and sig in Eval Computation
###         - -30052012: change my.ev2 (adapt it as est.my.ev2)
###         - -01062012: output one-sided and two-sided pvals
###         - -07062012: Pval=NA if p>n-1
###         - -03072012: new function: beta.mat, q.matrix3, est2.ww.mat2, est2.my.ev2
###         - 20072012: cleaning up and get rid of old functions
###         - 24072012: mistake in est2.my.ev2: if (dimf<dimg) ....
###         - 01082012: use coef(fit.cv) to extract coefficients from glmnet
###         - 09082012: add single-split method
###         - 11122012: diffregr.r originates from twosample_highdimregr-09082012.r

##Packages
library(mvtnorm)
library(glmnet)
library(CompQuadForm)

#######################
##-----Screening-----##
#######################
lasso.cvmin <- function(x,y){
  fit.cv <- cv.glmnet(x,y)
  beta <- as.numeric(coef(fit.cv,s='lambda.min')[-1])
  return(which(beta!=0))
}

##' Screening using crossvalidation & truncation
##'
##' Computes crossvalidated coefficients; if more nonzeros than n/k then truncate smallest coefficients
##' @title Screening using crossvalidation & truncation
##' @param x 
##' @param y 
##' @param k 
##' @return 
##' @author n.stadler
lasso.cvtrunc <- function(x,y,k.trunc=5){
  n <- nrow(x)
  fit.cv <- cv.glmnet(x,y)
  beta <- as.numeric(coef(fit.cv,s='lambda.min')[-1])
  beta[-order(abs(beta),decreasing=TRUE)[1:ceiling(n/k.trunc)]] <- 0
  return(which(beta!=0))
}

lasso.cv1se <- function(x,y){
  fit.cv <- cv.glmnet(x,y)
  beta <- as.numeric(coef(fit.cv,s='lambda.1se')[-1])
  return(which(beta!=0))
}

##############################
##--------P-VALUES----------##
##############################
logratio <- function(y1,y2,y,xx1,xx2,xx,beta1,beta2,beta){
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

ww.mat <- function(Sig,act,act1,act2){
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
  eval <- as.real(eigen(mat)$values)
  eval[abs(eval)<10^{-6}] <- 0
  return(list(ww.mat=mat,eval=eval))
}

ww.mat2 <- function(Sig,act,act1,act2){
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
  eval.mu <- as.real(eigen(mat)$values)
  eval2 <- 1-eval.mu
  eval2[abs(eval2)<10^{-6}] <- 0
  eval <- c(eval,sqrt(eval2),-sqrt(eval2))
  return(list(ww.mat=mat,eval=eval))
}

q.matrix <- function(Sig,a,b,s){
  ##Compute Q: based on true information matrix (Sig)
  if(length(s)==0){
    return(Sig[a,b,drop=FALSE])
  }else{
    return(Sig[a,b,drop=FALSE]-Sig[a,s,drop=FALSE]%*%solve(Sig[s,s,drop=FALSE])%*%Sig[s,b,drop=FALSE])
  }
}

my.ev2 <- function(Sig,act,act1,act2){
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
      qcc <- q.matrix(Sig,cc,cc,ss)
      aux.mat <- matrix(0,length(cc),length(cc))
      if(length(aa)!=0){
        qac <- q.matrix(Sig,aa,cc,ss)
        qaa <- q.matrix(Sig,aa,aa,ss)
        aux.mat <- aux.mat+(t(qac)%*%solve(qaa)%*%qac)%*%solve(qcc)
      }
      if(length(bb)!=0){
        qbc <- q.matrix(Sig,bb,cc,ss)
        qbb <- q.matrix(Sig,bb,bb,ss)
        aux.mat <- aux.mat+(t(qbc)%*%solve(qbb)%*%qbc)%*%solve(qcc)
      }
      ev.aux.complex <- eigen(aux.mat)$values
      ev.aux <- as.real(ev.aux.complex)
      ev.aux <-  sqrt(1-ev.aux/2)
    }
    eval<-rep(1,dimf-dimg)
    eval <- c(eval,rep(0,2*(length(ss)+1)),ev.aux,-ev.aux)
  }# end if (dimf>=dimg){
  if (dimf<dimg){
    if (length(cc)!=0){
      qcc <- q.matrix(Sig,cc,cc,ss)
      if(length(aa)!=0){
        qac <- q.matrix(Sig,aa,cc,ss)
        qaa <- q.matrix(Sig,aa,aa,ss)
        oaa <- qac%*%solve(qcc)%*%t(qac)
        aux.mat.aa<- oaa%*%solve(qaa)
        aux.mat <- aux.mat.aa
      }
      if(length(bb)!=0){
        qbc <- q.matrix(Sig,bb,cc,ss)
        qbb <- q.matrix(Sig,bb,bb,ss)
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
      ev.aux <- as.real(ev.aux.complex)
      ev.aux <-  sqrt(1-ev.aux/2)
    }
    eval<-rep(-1,dimg-dimf)
    eval <- c(eval,rep(-1,(length(ss)+1)),rep(1,(length(ss)+1)),rep(0,2*(length(ss)+1)),ev.aux,-ev.aux)
  }
  return(list(eval=eval,ev.aux.complex=ev.aux.complex))
}

beta.mat<-function(ind1,ind2,beta1,beta2,beta,sig1,sig2,sig,Sig){
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

q.matrix3 <- function(beta.a,beta.b,beta,sig.a,sig.b,sig,Sig,act.a,act.b,ss){

 ##Estimate Q
    b.ab<-beta.mat(act.a,act.b,beta.a,beta.b,beta,sig.a,sig.b,sig,Sig)
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

est2.my.ev2 <- function(y1,y2,x1,x2,beta1,beta2,beta,act1,act2,act){

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
      qcc <- q.matrix3(exp.beta,exp.beta,exp.beta1,sig,sig,sig1,var(x1),act,act,ss)+q.matrix3(exp.beta,exp.beta,exp.beta2,sig,sig,sig2,var(x2),act,act,ss)
      aux.mat <- matrix(0,length(cc),length(cc))
      if(length(aa)!=0){
        qac <- q.matrix3(exp.beta1,exp.beta,exp.beta1,sig1,sig,sig1,var(x1),act1,act,ss)
        qaa <- q.matrix3(exp.beta1,exp.beta1,exp.beta1,sig1,sig1,sig1,var(x1),act1,act1,ss)
        aux.mat <- aux.mat+(t(qac)%*%solve(qaa)%*%qac)%*%solve(qcc)
      }
      if(length(bb)!=0){
        qbc <- q.matrix3(exp.beta2,exp.beta,exp.beta2,sig2,sig,sig2,var(x2),act2,act,ss)
        qbb <- q.matrix3(exp.beta2,exp.beta2,exp.beta2,sig2,sig2,sig2,var(x2),act2,act2,ss)
        aux.mat <- aux.mat+(t(qbc)%*%solve(qbb)%*%qbc)%*%solve(qcc)
      }
      ev.aux.complex <- eigen(aux.mat)$values
      ev.aux <- as.real(ev.aux.complex)
      ev.aux <-  sqrt(1-ev.aux)
    }
    eval<-rep(1,dimf-dimg)
    eval <- c(eval,rep(0,2*(length(ss)+1)),ev.aux,-ev.aux)
  }# end if (dimf>=dimg){
  if (dimf<dimg){
    if (length(cc)!=0){
      qcc <- q.matrix3(exp.beta,exp.beta,exp.beta1,sig,sig,sig1,var(x1),act,act,ss)+q.matrix3(exp.beta,exp.beta,exp.beta2,sig,sig,sig2,var(x2),act,act,ss)
      aux.mat <- matrix(0,length(cc),length(cc))
      if(length(aa)!=0){
        qac <- q.matrix3(exp.beta1,exp.beta,exp.beta1,sig1,sig,sig1,var(x1),act1,act,ss)
        qaa <- q.matrix3(exp.beta1,exp.beta1,exp.beta1,sig1,sig1,sig1,var(x1),act1,act1,ss)
        oaa <- qac%*%solve(qcc)%*%t(qac)
        aux.mat.aa<- oaa%*%solve(qaa)
        aux.mat <- aux.mat.aa
      }
      if(length(bb)!=0){
        qbc <- q.matrix3(exp.beta2,exp.beta,exp.beta2,sig2,sig,sig2,var(x2),act2,act,ss)
        qbb <- q.matrix3(exp.beta2,exp.beta2,exp.beta2,sig2,sig2,sig2,var(x2),act2,act2,ss)
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
    eval <- c(eval,rep(-1,(length(ss)+1)),rep(1,(length(ss)+1)),rep(0,2*(length(ss)+1)),ev.aux,-ev.aux)
  }# end if (dimf<dimg){
  return(list(eval=eval,ev.aux.complex=ev.aux.complex))
}

agg.pval <- function(gamma,pval){
    min(quantile(pval/gamma,probs=gamma),1)
}

diffregr_pval <- function(y1,y2,x1,x2,beta1,beta2,beta,act1,act2,act,compute.evals,acc,verbose){
  ##########################
  ##compute test-statistic##
  ##########################
  teststat <- logratio(y1,y2,c(y1,y2),x1[,act1,drop=FALSE],x2[,act2,drop=FALSE],rbind(x1,x2)[,act,drop=FALSE],beta1,beta2,beta)
  #################################
  ##compute weights of sum-w-chi2##
  #################################
  weights.nulldistr <- eval(as.name(compute.evals))(y1,y2,x1,x2,beta1,beta2,beta,act1,act2,act)$eval
  weights.nulldistr <- weights.nulldistr[weights.nulldistr!=0]
  if (any(is.na(weights.nulldistr))){
    cat('warning: weight with value NA; pval=NA','\n')
    pval.onesided <- pval.twosided <- NA
  }else{
    pval.onesided <- davies(teststat,lambda=weights.nulldistr,acc=acc)$Qq;if(verbose){cat('ifault(davies):',davies(teststat,lambda=weights.nulldistr,acc=acc)$ifault,'\n')}
    pval.twosided <- 2*min(pval.onesided,1-pval.onesided)
  }
  return(list(pval.onesided=pval.onesided,pval.twosided=pval.twosided,weights.nulldistr=weights.nulldistr,teststat=teststat))
}

diffregr_singlesplit<- function(y1,y2,x1,x2,split1,split2,screen.meth='lasso.cvmin',
                                compute.evals='est2.my.ev2',
                                acc=1e-04,verbose=FALSE){
  
  ##Multisplit Pvalues
  ##
  ##Input:
  ##
  ## -Data:y1,y2,x1,x2 (centered)
  ## -b.splits: number of splits
  ## -frac.split: fraction train-data (for model selection) / test-data (for pval calculations)
  ## -screen.meth={lasso.cvmin,lasso.cv1se}
  ## -gamma.min: see Meinshausen&Meier&Buehlmann
  ## -compute.evals: 'est2.my.ev2'

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
  active[['modJ']] <- eval(as.name(screen.meth))(xx.train,yy.train)
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
    active[[paste('modIpop',j,sep='')]] <- eval(as.name(screen.meth))(xx.train,yy.train)

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
  if (any(l.act==0)){ cat('warning: at least one active-set is empty','\n')}
  if (all(l.act<c(n1.valid+n2.valid,n1.valid,n2.valid))){
    res.pval <- diffregr_pval(y1=y1[-split1],y2=y2[-split2],
                              x1=x1[-split1,,drop=FALSE],x2=x2[-split2,,drop=FALSE],
                              beta1=est.beta[['modIpop1']],beta2=est.beta[['modIpop2']],beta=est.beta[['modJ']],
                              active[['modIpop1']],active[['modIpop2']],active[['modJ']],
                              compute.evals,acc,verbose)
  }else{
    cat('warning: dim(model) > n: pval=NA','\n')
    res.pval <- list(pval.onesided=NA,pval.twosided=NA,weights.nulldistr=NA,teststat=NA)
  }
  
  return(list(pval.onesided=res.pval$pval.onesided,pval.twosided=res.pval$pval.twosided,
              teststat=res.pval$teststat,weights.nulldistr=res.pval$weights.nulldistr,
              active=active,beta=est.beta))
}

diffregr_multisplit<- function(y1,y2,x1,x2,b.splits=50,frac.split=1/2,screen.meth='lasso.cvmin',
                              gamma.min=0.05,compute.evals='est2.my.ev2',acc=1e-04,verbose=FALSE){

  n1 <- nrow(x1)
  n2 <- nrow(x2)
  
  res.multisplit <- lapply(seq(b.splits),
                          function(i){
                            split1 <- sample(1:n1,round(n1*frac.split),replace=FALSE)
                            split2 <- sample(1:n2,round(n2*frac.split),replace=FALSE)
                            res.singlesplit <- diffregr_singlesplit(y1,y2,x1,x2,split1,split2,screen.meth,
                                                                    compute.evals,acc,verbose)                      
                          })
  pval.onesided <- sapply(res.multisplit,function(x){x[['pval.onesided']]},simplify='array')
  pval.twosided <- sapply(res.multisplit,function(x){x[['pval.twosided']]},simplify='array')
  teststat <- sapply(res.multisplit,function(x){x[['teststat']]},simplify='array')
  weights.nulldistr <- sapply(res.multisplit,function(x){x[['weights.nulldistr']]},simplify='array')
  aggpval.onesided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.onesided[!is.na(pval.onesided)])$objective)
  aggpval.twosided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.twosided[!is.na(pval.twosided)])$objective)
    
  return(list(pval.onesided=pval.onesided,pval.twosided=pval.twosided,
              sspval.onesided=pval.onesided[1],sspval.twosided=pval.twosided[1],
              medpval.onesided=median(pval.onesided,na.rm=TRUE),medpval.twosided=median(pval.twosided,na.rm=TRUE),
              aggpval.onesided=aggpval.onesided,aggpval.twosided=aggpval.twosided,
              teststat=teststat,weights.nulldistr=weights.nulldistr,
              active.last=res.multisplit[[b.splits]]$active,beta.last=res.multisplit[[b.splits]]$beta))
}

twosample_single_regr <- function(y1,y2,x1,x2,n.screen.pop1=100,n.screen.pop2=100,screen.meth=lasso.cvmin,compute.evals='est2.my.ev2'){

  ##Single-Split Pvalues
  ##
  ##Input:
  ##
  ## -Data:y1,y2,x1,x2
  ## -n.screen.pop1, n.screen.pop2: no samples for screening
  ## -screen.meth={lasso.cvmin,lasso.cv1se}
  ## -compute.evals: 'est2.my.ev2'

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
      
    teststat <- logratio(y1[-split1],y2[-split2],c(y1[-split1],y2[-split2]),
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


## twosample_regr_oracle.eval<- function(y1,y2,x1,x2,b.splits=50,frac.split=1/2,lambda.cv='lambda.min',gamma.min=0.05,Sig){

##   ##Multisplit Pvalues
##   ##
##   ##Input:
##   ##
##   ## -Data:y1,y2,x1,x2
##   ## -b.splits: number of splits
##   ## -frac.split: fraction train-data (for model selection) / test-data (for pval calculations)
##   ## -lambda.cv={lambda.min,lambda.1se}
##   ## -gamma.min: see Meinshausen&Meier&Buehlmann
##   ## -Sig: true (under H0) information matrix

##   n1 <- nrow(x1)
##   n2 <- nrow(x2)
##   pval.onesided <- pval.twosided<- rep(NA,length=b.splits)

##   for (i in 1:b.splits){
##     split1 <- sample(1:n1,round(n1*frac.split),replace=FALSE)
##     split2 <- sample(1:n2,round(n2*frac.split),replace=FALSE)

##     est.beta <- active <- list()

##     ##model joint:
##     xx.train <- rbind(x1[split1,],x2[split2,])
##     yy.train <- c(y1[split1],y2[split2])
##     xx.valid <- rbind(x1[-split1,],x2[-split2,])
##     yy.valid <- c(y1[-split1],y2[-split2])
##     fit.cv <- cv.glmnet(xx.train,yy.train)
##     fit <- glmnet(xx.train,yy.train,lambda=fit.cv[[lambda.cv]])
##     active[['modJ']] <- which(as.numeric(coef(fit)[-1])!=0)
##     if(length(active[['modJ']])!=0){
##       est.beta[['modJ']] <- as.numeric(coef(lm(yy.valid~xx.valid[,active[['modJ']]]-1)))
##     }else{
##       est.beta[['modJ']]<-numeric(0)
##     }

##     ##model individual:
##     for (j in c('1','2')){
##       split.train <- eval(as.name(paste('split',j,sep='')))
##       xx.train <- eval(as.name(paste('x',j,sep='')))[split.train,]
##       yy.train <- eval(as.name(paste('y',j,sep='')))[split.train]
##       xx.valid <- eval(as.name(paste('x',j,sep='')))[-split.train,]
##       yy.valid <- eval(as.name(paste('y',j,sep='')))[-split.train]
##       fit.cv <- cv.glmnet(xx.train,yy.train)
##       fit <- glmnet(xx.train,yy.train,lambda=fit.cv[[lambda.cv]])
##       active[[paste('modIpop',j,sep='')]] <- which(as.numeric(coef(fit)[-1])!=0)

##       if(length(active[[paste('modIpop',j,sep='')]])!=0){
##         est.beta[[paste('modIpop',j,sep='')]] <- as.numeric(coef(lm(yy.valid~xx.valid[,active[[paste('modIpop',j,sep='')]]]-1)))
##       }else{
##         est.beta[[paste('modIpop',j,sep='')]]<-numeric(0)
##       }
##     }
      
##     l.act <- lapply(active,length)
##     n1.valid <- nrow(x1[-split1,])
##     n2.valid <- nrow(x2[-split2,])
##     if (any(l.act==0)){ warning('at least one active-set is empty')}
    
##     if (all(l.act< c(n1.valid+n2.valid,n1.valid,n2.valid))){
      
##       teststat <- logratio(y1[-split1],y2[-split2],c(y1[-split1],y2[-split2]),
##                            x1[-split1,active[['modIpop1']],drop=FALSE],x2[-split2,active[['modIpop2']],drop=FALSE],
##                            rbind(x1[-split1,active[['modJ']],drop=FALSE],x2[-split2,active[['modJ']],drop=FALSE]),
##                            est.beta[['modIpop1']],est.beta[['modIpop2']],est.beta[['modJ']])
    

##       ev.nulldistr <- my.ev2(Sig=Sig,act=active[['modJ']],act1=active[['modIpop1']],act2=active[['modIpop2']])$eval

##       if (any(is.na(ev.nulldistr))){
##         warning('Eigenval is NA: pval=NA')
##       }else{
##         pval_onesided <- davies(teststat,lambda=ev.nulldistr)$Qq
##         pval.onesided[i] <- pval_onesided
##         pval.twosided[i] <- 2*min(pval_onesided,1-pval_onesided)
##       }
##     }else{warning('dim(model) > n-1: pval=NA')}
##   }
##   sspval.onesided<-pval.onesided[1]
##   sspval.twosided<-pval.twosided[1]
##   aggpval.onesided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.onesided[!is.na(pval.onesided)])$objective)
##   aggpval.twosided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.twosided[!is.na(pval.twosided)])$objective)
  
##   return(list(pval.onesided=pval.onesided,pval.twosided=pval.twosided,sspval.onesided=sspval.onesided,sspval.twosided=sspval.twosided,
##               aggpval.onesided=aggpval.onesided,aggpval.twosided=aggpval.twosided,
##               LR.last=teststat,active.last=active,beta.last=est.beta))
## }

## q.matrix2 <- function(y,x,beta.a,beta.b,act.a,act.b,ss){
  
##  ##Estimate Q 
##     sc.a<-est.score.onlyregrcoef(y,x,beta.a,act.a)
##     sc.b<-est.score.onlyregrcoef(y,x,beta.b,act.b)
##     var.ab<-var(sc.a,sc.b)
##     aa<-seq(1,length(act.a))[!(act.a%in%ss)]
##     bb<-seq(1,length(act.b))[!(act.b%in%ss)]
##     s.a<-seq(1,length(act.a))[(act.a%in%ss)]
##     s.b<-seq(1,length(act.b))[(act.b%in%ss)]
##     if(length(ss)==0){
##         return(var.ab[aa,bb,drop=FALSE])
##     }else{
##         return(var.ab[aa,bb,drop=FALSE]-(var.ab[aa,s.b,drop=FALSE]%*%solve(var.ab[s.a,s.b,drop=FALSE])%*%var.ab[s.a,bb,drop=FALSE]))
##     }
## }

## est.my.ev2 <- function(y1,y2,x1,x2,beta1,beta2,beta,act1,act2,act){

##    ##Estimate Evals ('2nd order' simplification of W)
##    ##
##    ##Input:
##    ##  data: y1,y2,x1,x2
##    ##  beta: estimate for joint model
##    ##  beta1,beta2: estimates for individual model
##    ##  act: active beta's for joint model
##    ##  act1,act2: active beta's for individual models

##   #dimension of models
##   dimf1 <- length(act1)+1
##   dimf2 <- length(act2)+1
##   dimf <- dimf1+dimf2
##   dimg <- length(act)+1
  
##   #intersection of models
##   ss <- intersect(act,intersect(act1,act2))
  
##   if(length(ss)==0){warning('no intersection between models')}
  
##   aa <- setdiff(act1,ss)
##   bb <- setdiff(act2,ss)
##   cc <- setdiff(act,ss)

##   ev.aux <- ev.aux.complex<-numeric(0)
##   if (dimf>=dimg){
##     if (length(cc)!=0){
##       qcc <- q.matrix2(c(y1,y2),rbind(x1,x2),beta,beta,act,act,ss)
##       aux.mat <- matrix(0,length(cc),length(cc))
##       if(length(aa)!=0){
##         qac <- q.matrix2(y1,x1,beta1,beta,act1,act,ss)
##         qaa <- q.matrix2(y1,x1,beta1,beta1,act1,act1,ss)
##         aux.mat <- aux.mat+(t(qac)%*%solve(qaa)%*%qac)%*%solve(qcc)
##       }
##       if(length(bb)!=0){
##         qbc <- q.matrix2(y2,x2,beta2,beta,act2,act,ss)
##         qbb <- q.matrix2(y2,x2,beta2,beta2,act2,act2,ss)
##         aux.mat <- aux.mat+(t(qbc)%*%solve(qbb)%*%qbc)%*%solve(qcc)
##       }
##       ev.aux.complex <- eigen(aux.mat)$values
##       ev.aux <- as.real(ev.aux.complex)
##       ev.aux <-  sqrt(1-ev.aux/2)
##     }
##     eval<-rep(1,dimf-dimg)
##     eval <- c(eval,rep(0,2*(length(ss)+1)),ev.aux,-ev.aux)
##   }# end if (dimf>=dimg){

##   if (dimf<dimg){
##     if (length(cc)!=0){
##       qcc <- q.matrix2(c(y1,y2),rbind(x1,x2),beta,beta,act,act,ss)
##       if(length(aa)!=0){
##         qac <- q.matrix2(y1,x1,beta1,beta,act1,act,ss)
##         qaa <- q.matrix2(y1,x1,beta1,beta1,act1,act1,ss)
##         oaa <- qac%*%solve(qcc)%*%t(qac)
##         aux.mat.aa<- oaa%*%solve(qaa)
##         aux.mat <- aux.mat.aa
##       }
##       if(length(bb)!=0){
##         qbc <- q.matrix2(y2,x2,beta2,beta,act2,act,ss)
##         qbb <- q.matrix2(y2,x2,beta2,beta2,act2,act2,ss)
##         obb <- qbc%*%solve(qcc)%*%t(qbc)
##         aux.mat.bb<- obb%*%solve(qbb)
##         aux.mat <- aux.mat.bb
##       }
##       if((length(aa)!=0)&(length(bb)!=0)){
##         oab <- qac%*%solve(qcc)%*%t(qbc)
##         aux.mat<- rbind(cbind(aux.mat.aa,oab%*%solve(qbb)),cbind(t(oab)%*%solve(qaa),aux.mat.bb))
##       }
##     }
##     if ((length(aa)!=0)|(length(bb)!=0)){ ##if (length(aa)==0)&(length(bb)==0) 'MI included in MJ; therefore -X^2(dim(MJ)-dim(MI)) distributed'
##       ev.aux.complex <- eigen(aux.mat)$values
##       ev.aux <- as.real(ev.aux.complex)
##       ev.aux <-  sqrt(1-ev.aux/2)
##     }
##     eval<-rep(-1,dimg-dimf)
##     eval <- c(eval,rep(-1,(length(ss)+1)),rep(1,(length(ss)+1)),rep(0,2*(length(ss)+1)),ev.aux,-ev.aux)
##   }
##   return(list(eval=eval,ev.aux.complex=ev.aux.complex))
## }

## est.score.onlyregrcoef<-function(y,x,beta,act){
##  ##Estimate Score-Function
 
##     n <- length(y)
##     res <- as.numeric(y-x[,act,drop=FALSE]%*%beta)
##     sig <- sum(res^2)/n
##     res*x[,act,drop=FALSE]/sig

## }

## est.score<-function(y,x,beta,act){
##  ##Estimate Score-Function
##   n <- length(y)
##   if(length(act)!=0){
##     res <- as.numeric(y-x[,act,drop=FALSE]%*%beta)
##     sig <- sum(res^2)/n
##     est.score <- cbind(res*x[,act,drop=FALSE]/sig,((res^2/sig)-1)/(2*sig))
##   }else{
##     res <- y
##     sig <- sum(res^2)/n
##     est.score <- ((res^2/sig)-1)/(2*sig)
##   }
##   est.score
## }

## est.ww.mat <- function(y1,y2,x1,x2,beta1,beta2,beta,act1,act2,act){

##    ##Estimate W and Eval(W) (without simplification of W)
##    ##
##    ##Input:
##    ##  data: y1,y2,x1,x2
##    ##  beta: estimate for joint model
##    ##  beta1,beta2: estimates for individual model
##    ##  act: active variables for joint model
##    ##  act1,act2: active variables for individual models
  
##     sc.i.1<- est.score(y1,x1,beta1,act1)
##     sc.i.2<- est.score(y2,x2,beta2,act2)
##     sc.j.1<- est.score(y1,x1,beta,act)
##     sc.j.2<- est.score(y2,x2,beta,act)

##     bfg <- rbind(var(sc.i.1,sc.j.1),var(sc.i.2,sc.j.2))
##     bf <- matrix(0,length(act1)+length(act2)+2,length(act1)+length(act2)+2)
##     bf[1:(length(act1)+1),1:(length(act1)+1)] <- var(sc.i.1)
##     bf[(length(act1)+1)+(1:(length(act2)+1)),(length(act1)+1)+(1:(length(act2)+1))] <- var(sc.i.2)
##     bg <- var(sc.j.1)+var(sc.j.2)
##     mat <- rbind(cbind(diag(1,length(act1)+length(act2)+2),bfg%*%solve(bg)),cbind(-t(bfg)%*%solve(bf),diag(-1,length(act)+1)))

##     eval <- as.real(eigen(mat)$values)
##     eval[abs(eval)<10^{-6}] <- 0

##     return(list(ww.mat=mat,eval=eval))
## }

## est.ww.mat2 <- function(y1,y2,x1,x2,beta1,beta2,beta,act1,act2,act){

##    ##Estimate W and Eval(W) ('1st order' simplification of W)
##    ##
##    ##Input:
##    ##  data: y1,y2,x1,x2
##    ##  beta: estimate for joint model
##    ##  beta1,beta2: estimates for individual model
##    ##  act: active beta's for joint model
##    ##  act1,act2: active beta's for individual models
##     sc.i.1<- est.score(y1,x1,beta1,act1)
##     sc.i.2<- est.score(y2,x2,beta2,act2)
##     sc.j.1<- est.score(y1,x1,beta,act)
##     sc.j.2<- est.score(y2,x2,beta,act)

##     dimf1 <- length(act1)+1
##     dimf2 <- length(act2)+1
##     dimf <- dimf1+dimf2
##     dimg <- length(act)+1

##     bfg <- rbind(var(sc.i.1,sc.j.1),var(sc.i.2,sc.j.2))
##     bgf <- t(bfg)
##     bf <- matrix(0,dimf,dimf)
##     bf[1:dimf1,1:dimf1] <- var(sc.i.1)
##     bf[dimf1+(1:dimf2),dimf1+(1:dimf2)] <- var(sc.i.2)
##     bg <- var(sc.j.1)+var(sc.j.2)

##     if (dimf>=dimg){
##         mat <- bgf%*%solve(bf)%*%bfg%*%solve(bg)
##         eval<-rep(1,dimf-dimg)
##     }
##     if (dimf<dimg){
##         mat <- bfg%*%solve(bg)%*%bgf%*%solve(bf)
##         eval<-rep(-1,dimg-dimf)
##     }
##     eval.mu.complex<-eigen(mat)$values
##     eval.mu <- as.real(eval.mu.complex)
##     eval <- c(eval,sqrt(1-eval.mu),-sqrt(1-eval.mu))

##     return(list(ww.mat=mat,eval=eval,eval.mu.complex=eval.mu.complex))
## }

## my.ev <- function(Sig,act,act1,act2){

##    ##Compute Eval (with simplification of W)
##    ##Sig: E_0[s(X)s(X)']
##    ##act: active variables for joint model
##    ##act1,act2: active variables for individual models
##     ss <- intersect(act,intersect(act1,act2))
##     aa <- setdiff(act1,ss)
##     bb <- setdiff(act2,ss)
##     cc <- setdiff(act,ss)

##     ##initialization
##     rank <- 0
##     ev.aux <- numeric(0)

##     if (length(cc)!=0){
##         qcc <- q.matrix(Sig,cc,cc,ss)
##         aux.mat <- matrix(0,length(cc),length(cc))
##         if(length(aa)!=0){
##             qac <- q.matrix(Sig,aa,cc,ss)
##             qaa <- q.matrix(Sig,aa,aa,ss)
##             aux.mat <- aux.mat+(t(qac)%*%solve(qaa)%*%qac)%*%solve(qcc)
##         }
##         if(length(bb)!=0){
##             qbc <- q.matrix(Sig,bb,cc,ss)
##             qbb <- q.matrix(Sig,bb,bb,ss)
##             aux.mat <- aux.mat+(t(qbc)%*%solve(qbb)%*%qbc)%*%solve(qcc)
##         }
##         ev.aux <- as.real(eigen(aux.mat)$values)
##         ev.aux[abs(ev.aux)<10^{-6}] <-0
##         rank <- sum(ev.aux!=0)
##         ev.aux <-  sqrt(1-ev.aux[ev.aux!=0]/2)
##     }
##     eigenval <- c(rep(0,2*length(ss)),ev.aux,-ev.aux,rep(1,length=length(act1)+length(act2)-length(ss)-rank),rep(-1,length=length(act)-length(ss)-rank))

##     return(list(eval=c(1,0,0,eigenval)))#1,0,0 for sig1,sig2 and sig

## }

## est2.ww.mat2 <- function(y1,y2,x1,x2,beta1,beta2,beta,act1,act2,act){

##   ##Estimate W and Eval(W) ('1st order' simplification of W)
##   ##
##   ##Input:
##   ##  data: y1,y2,x1,x2
##   ##  beta: estimate for joint model
##   ##  beta1,beta2: estimates for individual model
##   ##  act: active beta's for joint model
##   ##  act1,act2: active beta's for individual models
  
##   n1 <- length(y1)
##   n2 <- length(y2)
##   n <- length(y)
##   mu1<-mu2<-mu<-0
##   sig1 <- (sum(y1^2)/n1);sig2 <- (sum(y2^2)/n2);sig <- (sum(y^2)/n)
##   if(length(beta1)!=0){
##     mu1<-x1[,act1,drop=FALSE]%*%beta1
##     sig1 <- sum((y1-mu1)^2)/n1
##   }
##   if(length(beta2)!=0){
##     mu2<-x2[,act2,drop=FALSE]%*%beta2
##     sig2 <- sum((y2-mu2)^2)/n2
##   }
##   if(length(beta)!=0){
##     mu<-rbind(x1,x2)[,act,drop=FALSE]%*%beta
##     sig <- sum((y-mu)^2)/n
##   }
##   bact.pop1<- beta.mat(beta,beta,beta1,sig,sig,sig1,var(x1))
##   bact.pop2<- beta.mat(beta,beta,beta2,sig,sig,sig2,var(x2))
##   bact1<- beta.mat(beta1,beta1,beta1,sig1,sig1,sig1,var(x1))
##   bact2<- beta.mat(beta2,beta2,beta2,sig2,sig2,sig2,var(x2))
##   bact1.act<- beta.mat(beta1,beta,beta1,sig1,sig,sig1,var(x1))
##   bact2.act<- beta.mat(beta2,beta,beta2,sig2,sig,sig2,var(x2))
    
##   dimf1 <- length(act1)+1
##   dimf2 <- length(act2)+1
##   dimf <- dimf1+dimf2
##   dimg <- length(act)+1
  
##   bfg <- rbind(bact1.act,bact2.act)
##   bgf <- t(bfg)
##   bf <- matrix(0,dimf,dimf)
##   bf[1:dimf1,1:dimf1] <- bact1
##   bf[dimf1+(1:dimf2),dimf1+(1:dimf2)] <- bact2
##   bg <- bact.pop1+bact.pop2
  
##   if (dimf>=dimg){
##     mat <- bgf%*%solve(bf)%*%bfg%*%solve(bg)
##     eval<-rep(1,dimf-dimg)
##   }
##   if (dimf<dimg){
##     mat <- bfg%*%solve(bg)%*%bgf%*%solve(bf)
##     eval<-rep(-1,dimg-dimf)
##   }
##   eval.mu.complex<-eigen(mat)$values
##   eval.mu <- as.real(eval.mu.complex)
##   eval <- c(eval,sqrt(1-eval.mu),-sqrt(1-eval.mu))
  
##   return(list(ww.mat=mat,eval=eval,eval.mu.complex=eval.mu.complex))
## }
