###TwoSampleTest for Gaussian Graphical Model: assume mu1=mu2=mu0=0 
###Date: 20/07/2012
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


##Packages
library(mvtnorm)
library(glasso)
library(CompQuadForm)

#############
##Screening##
#############

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

glasso.parcor.launi <- function(x,maxiter=1000,term=10^{-3}){
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
  list(w=gl$w,wi=gl$wi,iter=iter)
}

glasso.launi <- function(x){
  s <- var(x)
  rho.uni <- sqrt(2*log(ncol(x))/nrow(x))
  gl <- glasso(s,rho=rho.uni,penalize.diagonal=FALSE)
  list(w=gl$w,wi=gl$wi)
}


############
##P-VALUES##
############

logratio <- function(x1,x2,x,sig1,sig2,sig,mu1,mu2,mu){
  ##Compute 2*log-likelihood ratio
  ##Input:
  ##-x1,x2,x 
  ##-sig1,sig2, sig: mle
  twiceLR <- 2*(sum(dmvnorm(x1,mean=mu1,sigma=sig1,log=TRUE))+sum(dmvnorm(x2,mean=mu2,sigma=sig2,log=TRUE))-sum(dmvnorm(x,mean=mu,sigma=sig,log=TRUE)))
  list(twiceLR=twiceLR,sig1=sig1,sig2=sig2,sig=sig)
}

inf.mat<-function(Sig,ind){
  ##Computes information matrix of Gaussian Graphical model with Sig=solve(InvCov)
  p <- ncol(Sig)
    uptri.rownr <- row(Sig)[upper.tri(Sig,diag=TRUE)]
    uptri.colnr <- col(Sig)[upper.tri(Sig,diag=TRUE)]

    infmat <- matrix(NA,p*(p+1)/2,p*(p+1)/2)
    for (j in ind){
        colnr.j <- uptri.colnr[j]
        rownr.j <- uptri.rownr[j]
        for (k in ind){
             colnr.k <- uptri.colnr[k]
             rownr.k <- uptri.rownr[k]
             infmat[j,k] <- Sig[colnr.j,rownr.k]*Sig[rownr.j,colnr.k]+Sig[colnr.j,colnr.k]*Sig[rownr.j,rownr.k]
         }
    }
    return(infmat)
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
  ##  act: active variables for joint model (including diag)
  ##  act1,act2: active variables for individual models (including diag)

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

q.matrix <- function(imat,a,b,s){
  ##Compute Q: based on true information matrix (imat)
  if(length(s)==0){
        return(imat[a,b,drop=FALSE])
    }else{
        return(imat[a,b,drop=FALSE]-imat[a,s,drop=FALSE]%*%solve(imat[s,s,drop=FALSE])%*%imat[s,b,drop=FALSE])
    }
}

my.ev2 <- function(Sig,act,act1,act2){
  ##Compute Eval (with 2nd order simplification of W)
  ##
  ##Input:
  ##  Sig: E_0[s(X)s(X)´]
  ##  act: active variables for joint model (including diag)
  ##  act1,act2: active variables for individual models (including diag)

  dimf1 <- length(act1)
  dimf2 <- length(act2)
  dimf <- dimf1+dimf2
  dimg <- length(act)

  ss <- intersect(act,intersect(act1,act2))

  if(length(ss)==0){warning('no intersection between models')}
  aa <- setdiff(act1,ss)
  bb <- setdiff(act2,ss)
  cc <- setdiff(act,ss)

  ev.aux <- ev.aux.complex <- numeric(0)
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
    eval <- c(eval,rep(0,2*length(ss)),ev.aux,-ev.aux)
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
    if ((length(aa)!=0)|(length(bb)!=0)){
      ev.aux.complex <- eigen(aux.mat)$values
      ev.aux <- as.real(ev.aux.complex)
      ev.aux <-  sqrt(1-ev.aux/2)
    }
    eval<-rep(-1,dimg-dimf)
    eval <- c(eval,rep(-1,length(ss)),rep(1,length(ss)),rep(0,2*length(ss)),ev.aux,-ev.aux)
  }
   return(list(eval=eval,ev.aux.complex=ev.aux.complex))
}

beta.mat<-function(ind1,ind2,sig1,sig2,sig){
  ##Compute beta-matrix
  p <- ncol(sig)
  uptri.rownr <- row(sig)[upper.tri(sig,diag=TRUE)]
  uptri.colnr <- col(sig)[upper.tri(sig,diag=TRUE)]

  bmat <- matrix(NA,p*(p+1)/2,p*(p+1)/2)
  for (j in ind1){
    colnr.j <- uptri.colnr[j]
    rownr.j <- uptri.rownr[j]
    for (k in ind2){
      colnr.k <- uptri.colnr[k]
      rownr.k <- uptri.rownr[k]
      bmat[j,k] <- sig[colnr.j,rownr.j]*sig[colnr.k,rownr.k]-sig1[colnr.j,rownr.j]*sig[colnr.k,rownr.k]-sig[colnr.j,rownr.j]*sig2[colnr.k,rownr.k]+sig1[colnr.j,rownr.j]*sig2[colnr.k,rownr.k]+sig[colnr.j,rownr.k]*sig[rownr.j,colnr.k]+sig[colnr.j,colnr.k]*sig[rownr.j,rownr.k]
    }
  }
  return(bmat[ind1,ind2,drop=FALSE])
}

est2.ww.mat2 <- function(sig1,sig2,sig,act1,act2,act){

  ##Estimate W and Eval(W) ('1st order' simplification of W)
  ##
  ##Input:
  ##  data: xx1,xx2
  ##  sig: estimate for joint model
  ##  sig1,sig2: estimates for individual model
  ##  act: active variables for joint model
  ##  act1,act2: active variables for individual models
  
  bact.pop1<- beta.mat(act,act,sig,sig,sig1)
  bact.pop2<- beta.mat(act,act,sig,sig,sig2)
  bact1<- beta.mat(act1,act1,sig1,sig1,sig1)
  bact2<- beta.mat(act2,act2,sig2,sig2,sig2)
  bact1.act<- beta.mat(act1,act,sig1,sig,sig1)
  bact2.act<- beta.mat(act2,act,sig2,sig,sig2)

  dimf1 <- length(act1)
  dimf2 <- length(act2)
  dimf <- dimf1+dimf2
  dimg <- length(act)
  if (any(c(dimf1,dimf2,dimg)>min(nrow(xx1),nrow(xx2)))){
    warning('|active-set| > n')
    return(list(eval=NA))
  }else{
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
    onemineval.mu <- 1-as.real(eval.mu.complex)
    onemineval.mu[abs(onemineval.mu)<10^{-6}]<-0
    eval <- c(eval,sqrt(onemineval.mu),-sqrt(onemineval.mu))
    
    return(list(ww.mat=mat,eval=eval,eval.mu.complex=eval.mu.complex))
    #return(list(ww.mat=mat,eval=eval,eval.mu.complex=eval.mu.complex,invbg=solve(bg),bfg=bfg))
  }
}

q.matrix3 <- function(sig,sig.a,sig.b,act.a,act.b,ss){
 ##Estimate Q 
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

est2.my.ev2 <- function(sig1,sig2,sig,act1,act2,act){

  ##Estimate Evals ('2nd order' simplification of W)
  ##
  ##Input:
  ##  data: x1,x2
  ##  sig1,sig2: estimate for individual model
  ##  sig: estimates for joint model
  ##  act1,act2: active variables for individual model
  ##  act: active variables for joint model
  
  #dimension of models
  dimf1 <- length(act1)
  dimf2 <- length(act2)
  dimf <- dimf1+dimf2
  dimg <- length(act)
  
  #intersection of models
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
  }# end if (dimf>=dimg){
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
  }
  return(list(eval=eval,ev.aux.complex=ev.aux.complex))
}
   
agg.pval <- function(gamma,pval){
    min(quantile(pval/gamma,probs=gamma),1)
}

twosample_diffnet2 <- function(x1,x2,b.splits=50,frac.split=1/2,screen.lambda='lambda.cv',gamma.min=0.05,compute.evals='est2.my.ev2',diag.invcov=TRUE,acc=1e-04,...){

  ##Multisplit Pvalues
  ##
  ##Input:
  ##
  ## -Data:x1,x2 [scaled data (var=1,mean=0), in particular mu1=mu2=mu=0]
  ## -b.splits: number of splits
  ## -frac.split: fraction train-data (for model selection) / test-data (for pval calculations)
  ## -screen.lambda={lambda.cv,lambda.cv.1se,lambda.uni.parcor,lambda.uni.invcov}
  ## -gamma.min: see Meinshausen&Meier&Buehlmann
  ## -compute.evals: 'est2.my.ev2'/'est2.ww.mat2 '
  ## -diag.invcov: TRUE (parameter of interest is invcov; diag(invcov) also considered) / FALSE (param of interest is par.cor; no diagonal)
  ## -...: additional arguments for cv.glasso(x,...) or glasso.parcor.launi(x,...)

  ##????Important Notes: Pval can be NA, because...
  ##????
  ##????                 1) Eval=sqrt(neg. value)=NA
  ##????                 2) W cannot be estimated (|active-set| > n, problems matrix inversion) 

  n1 <- nrow(x1)
  n2 <- nrow(x2)
  p <- ncol(x1)
  pval.onesided <- pval.twosided <- rep(NA,length=b.splits)

  if(screen.lambda=='lambda.cv'){screen.meth <- 'cv.glasso'}
  if(screen.lambda=='lambda.cv.1se'){screen.meth <- 'cv.glasso.1se'}
  if(screen.lambda=='lambda.uni.parcor'){screen.meth <- 'glasso.parcor.launi'}
  if(screen.lambda=='lambda.uni.invcov'){screen.meth <- 'glasso.launi'}

  if(diag.invcov){df.param <- p*(p+1)/2}
  if(!diag.invcov){df.param <- p*(p-1)/2}
  
  for (i in 1:b.splits){

    cat('split:',i,'\n')
    
    split1 <- sample(1:n1,round(n1*frac.split),replace=FALSE)
    split2 <- sample(1:n2,round(n2*frac.split),replace=FALSE)

    est.mu <- est.sig <- est.wi <- active <- list()#save est.beta & active variables

    ##model joint:
    xx.train <- rbind(x1[split1,],x2[split2,])
    xx.valid <- rbind(x1[-split1,],x2[-split2,])
    fit.screen <- eval(as.name(screen.meth))(xx.train,...)
    act <- which(fit.screen$wi[upper.tri(diag(1,p),diag=TRUE)]!=0)
    if(!diag.invcov){
      ind.diag <- which(diag(1,p)[upper.tri(diag(1,p),diag=TRUE)]!=0)
      act <- setdiff(act,ind.diag)
    }
    active[['modJ']] <- act
    est.mu[['modJ']] <- rep(0,p)
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
    
    ##model individual:
    for (j in c('1','2')){
      split.train <- eval(as.name(paste('split',j,sep='')))
      xx.train <- eval(as.name(paste('x',j,sep='')))[split.train,]
      xx.valid <- eval(as.name(paste('x',j,sep='')))[-split.train,]
      fit.screen <- eval(as.name(screen.meth))(xx.train,...)
      act <- which(fit.screen$wi[upper.tri(diag(1,p),diag=TRUE)]!=0)
      if(!diag.invcov){
        ind.diag <- which(diag(1,p)[upper.tri(diag(1,p),diag=TRUE)]!=0)
        act <- setdiff(act,ind.diag)
      }
      active[[paste('modIpop',j,sep='')]] <- act
      est.mu[[paste('modIpop',j,sep='')]] <- rep(0,p)
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
    l.act <- lapply(active,length)
    n1.valid <- nrow(x1[-split1,])
    n2.valid <- nrow(x2[-split2,])
    if (any(l.act==0)){ warning('at least one active-set is empty')}
    if (all(l.act>= c(n1.valid+n2.valid,n1.valid,n2.valid))){warning('dim(model) > n-1')}

    teststat <- logratio(x1=x1[-split1,,drop=FALSE],
                         x2=x2[-split2,,drop=FALSE],
                         x=rbind(x1[-split1,,drop=FALSE],x2[-split2,,drop=FALSE]),
                         sig1=est.sig[['modIpop1']],sig2=est.sig[['modIpop2']],sig=est.sig[['modJ']],
                         mu1=est.mu[['modIpop1']],mu2=est.mu[['modIpop2']],mu=est.mu[['modJ']])
    
    ev.nulldistr <- eval(as.name(compute.evals))(est.sig[['modIpop1']],est.sig[['modIpop2']],est.sig[['modJ']],
                                                 active[['modIpop1']],active[['modIpop2']],active[['modJ']])$eval
    
    if (any(is.na(ev.nulldistr))){
      warning('Eigenval is NA: pval=NA')
    }else{
      pval_onesided <- davies(teststat$twiceLR,lambda=ev.nulldistr,acc=acc)$Qq;cat('ifault(davies):',davies(teststat$twiceLR,lambda=ev.nulldistr,acc=acc)$ifault,'\n')
      pval.onesided[i] <- pval_onesided
      pval.twosided[i] <- 2*min(pval_onesided,1-pval_onesided)
    }
  }
  sspval.onesided<-pval.onesided[1]
  sspval.twosided<-pval.twosided[1]
  aggpval.onesided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.onesided[!is.na(pval.onesided)])$objective)
  aggpval.twosided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.twosided[!is.na(pval.twosided)])$objective)
    
  return(list(pval.onesided=pval.onesided,pval.twosided=pval.twosided,sspval.onesided=sspval.onesided,
              sspval.twosided=sspval.twosided,
              aggpval.onesided=aggpval.onesided,aggpval.twosided=aggpval.twosided,
              LR.last=teststat$twiceLR,active.last=active,sig.last=est.sig,wi.last=est.wi,ev.last=ev.nulldistr))
}


twosample_single_diffnet <- function(x1,x2,n.screen.pop1=150,n.screen.pop2=150,screen.lambda='lambda.cv',
                                     compute.evals='est2.my.ev2',diag.invcov=TRUE,acc=1e-04,...){

  ##Single-Split Pvalues
  ##
  ##Input:
  ##
  ## -Data:x1,x2 [scaled data (var=1,mean=0), in particular mu1=mu2=mu=0]
  ## -n.screen.pop1, n.screen.pop2: no samples for screening
  ## -screen.lambda={lambda.cv,lambda.cv.1se,lambda.uni.parcor,lambda.uni.invcov}
  ## -compute.evals: 'est2.my.ev2'/'est2.ww.mat2 '
  ## -diag.invcov: TRUE (parameter of interest is invcov; diag(invcov) also considered) / FALSE (param of interest is par.cor; no diagonal)
  ## -...: additional arguments for cv.glasso(x,...) or glasso.parcor.launi(x,...)

  n1 <- nrow(x1)
  n2 <- nrow(x2)
  p <- ncol(x1)

  if(screen.lambda=='lambda.cv'){screen.meth <- 'cv.glasso'}
  if(screen.lambda=='lambda.cv.1se'){screen.meth <- 'cv.glasso.1se'}
  if(screen.lambda=='lambda.uni.parcor'){screen.meth <- 'glasso.parcor.launi'}
  if(screen.lambda=='lambda.uni.invcov'){screen.meth <- 'glasso.launi'}

  if(diag.invcov){df.param <- p*(p+1)/2}
  if(!diag.invcov){df.param <- p*(p-1)/2}

  ##split data
  split1 <- sample(1:n1,n.screen.pop1,replace=FALSE)
  split2 <- sample(1:n2,n.screen.pop2,replace=FALSE)
  est.mu <- est.sig <- est.wi <- active <- list()#save est.beta & active variables

  ##model joint:
  xx.train <- rbind(x1[split1,],x2[split2,])
  xx.valid <- rbind(x1[-split1,],x2[-split2,])
  fit.screen <- eval(as.name(screen.meth))(xx.train,...)
  act <- which(fit.screen$wi[upper.tri(diag(1,p),diag=TRUE)]!=0)
  if(!diag.invcov){
    ind.diag <- which(diag(1,p)[upper.tri(diag(1,p),diag=TRUE)]!=0)
    act <- setdiff(act,ind.diag)
  }
  active[['modJ']] <- act
  est.mu[['modJ']] <- rep(0,p)
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
    
  ##model individual:
  for (j in c('1','2')){
    split.train <- eval(as.name(paste('split',j,sep='')))
    xx.train <- eval(as.name(paste('x',j,sep='')))[split.train,]
    xx.valid <- eval(as.name(paste('x',j,sep='')))[-split.train,]
    fit.screen <- eval(as.name(screen.meth))(xx.train,...)
    act <- which(fit.screen$wi[upper.tri(diag(1,p),diag=TRUE)]!=0)
    if(!diag.invcov){
      ind.diag <- which(diag(1,p)[upper.tri(diag(1,p),diag=TRUE)]!=0)
      act <- setdiff(act,ind.diag)
    }
    active[[paste('modIpop',j,sep='')]] <- act
    est.mu[[paste('modIpop',j,sep='')]] <- rep(0,p)
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
  l.act <- lapply(active,length)
  n1.valid <- nrow(x1[-split1,])
  n2.valid <- nrow(x2[-split2,])
  if (any(l.act==0)){ warning('at least one active-set is empty')}
  if (all(l.act>= c(n1.valid+n2.valid,n1.valid,n2.valid))){warning('dim(model) > n-1')}

  teststat <- logratio(x1=x1[-split1,,drop=FALSE],
                       x2=x2[-split2,,drop=FALSE],
                       x=rbind(x1[-split1,,drop=FALSE],x2[-split2,,drop=FALSE]),
                       sig1=est.sig[['modIpop1']],sig2=est.sig[['modIpop2']],sig=est.sig[['modJ']],
                       mu1=est.mu[['modIpop1']],mu2=est.mu[['modIpop2']],mu=est.mu[['modJ']])
  
  ev.nulldistr <- eval(as.name(compute.evals))(est.sig[['modIpop1']],est.sig[['modIpop2']],est.sig[['modJ']],
                                               active[['modIpop1']],active[['modIpop2']],active[['modJ']])$eval
  if (any(is.na(ev.nulldistr))){
    warning('Eigenval is NA: pval=NA')
  }else{
    sspval.onesided <- davies(teststat$twiceLR,lambda=ev.nulldistr,acc=acc)$Qq;cat('ifault(davies):',davies(teststat$twiceLR,lambda=ev.nulldistr,acc=acc)$ifault,'\n')
    sspval.twosided <- 2*min(sspval.onesided,1-sspval.onesided)
  }

  return(list(sspval.onesided=sspval.onesided,
              sspval.twosided=sspval.twosided,
              LR=teststat$twiceLR,active=active,sig=est.sig,wi=est.wi,ev=ev.nulldistr))
}


#################
##Old Functions##
#################

## est.score<-function(x,Sig,ind){
##    ##Estimate Score-Function
##    ##
##    ##x: data
##    ##Sig: estimated covariance Sig=solve(SigInv)
##    ##ind: active variables
##     p <- ncol(Sig)
##     n <- nrow(x)
##     uptri.rownr <- row(Sig)[upper.tri(Sig,diag=TRUE)]
##     uptri.colnr <- col(Sig)[upper.tri(Sig,diag=TRUE)]

##     sc <- matrix(NA,n,p*(p+1)/2)
##     for (j in ind){
##         colnr.j <- uptri.colnr[j]
##         rownr.j <- uptri.rownr[j]
##         sc[,j] <- x[,colnr.j]*x[,rownr.j]-Sig[colnr.j,rownr.j]
##     }
## return(sc[,ind,drop=FALSE])

## }

## q.matrix2 <- function(x,sig.a,sig.b,act.a,act.b,ss){
##  ##Estimate Q 
##     sc.a<-est.score(x,sig.a,act.a)
##     sc.b<-est.score(x,sig.b,act.b)
##     var.ab<-var(sc.a,sc.b)
##     aa<-seq(1,length(act.a))[!(act.a%in%ss)]
##     bb<-seq(1,length(act.b))[!(act.b%in%ss)]
##     s.a<-seq(1,length(act.a))[(act.a%in%ss)]
##     s.b<-seq(1,length(act.b))[(act.b%in%ss)]
    
##     if(length(ss)==0){
##         return(var.ab[aa,bb])
##     }else{
##         return(var.ab[aa,bb,drop=FALSE]-(var.ab[aa,s.b,drop=FALSE]%*%solve(var.ab[s.a,s.b,drop=FALSE])%*%var.ab[s.a,bb,drop=FALSE]))
##     }
## }

## est.ww.mat <- function(xx1,xx2,sig1,sig2,sig,act1,act2,act){

##   ##Estimate W and Eval(W) (without simplification of W)
##    ##
##    ##Input:
##    ##  data: xx1,xx2
##    ##  sig1,sig2: estimates for individual model
##    ##  sig: estimate for joint model
##    ##  act: active variables for joint model
##    ##  act1,act2: active variables for individual models
  
##     sc.i.1<- est.score(xx1,sig1,act1)
##     sc.i.2<- est.score(xx2,sig2,act2)
##     sc.j.1<- est.score(xx1,sig,act)
##     sc.j.2<- est.score(xx2,sig,act)

##     dimf1 <- length(act1)
##     dimf2 <- length(act2)
##     dimf <- dimf1+dimf2
##     dimg <- length(act)

##     if (any(c(dimf1,dimf2,dimg)>min(nrow(xx1),nrow(xx2)))){
##       warning('|active-set| > n')
##       return(list(eval=NA))
##     }else{

##       bfg <- rbind(var(sc.i.1,sc.j.1),var(sc.i.2,sc.j.2))
##       bf <- matrix(0,dimf,dimf)
##       bf[1:dimf1,1:dimf1] <- var(sc.i.1)
##       bf[dimf1+(1:dimf2),dimf1+(1:dimf2)] <- var(sc.i.2)
##       bg <- var(sc.j.1)+var(sc.j.2)
##       mat <- rbind(cbind(diag(1,dimf),bfg%*%solve(bg)),cbind(-t(bfg)%*%solve(bf),diag(-1,dimg)))

##       eval <- as.real(eigen(mat)$values)
##       eval[abs(eval)<10^{-6}] <- 0


##       return(list(ww.mat=mat,eval=eval))
##     }
## }

## est.ww.mat2 <- function(xx1,xx2,sig1,sig2,sig,act1,act2,act){

##    ##Estimate W and Eval(W) ('1st order' simplification of W)
##    ##
##    ##Input:
##    ##  data: xx1,xx2
##     ##  sig: estimate for joint model
##    ##  sig1,sig2: estimates for individual model
##    ##  act: active variables for joint model
##    ##  act1,act2: active variables for individual models
  
##     sc.i.1<- est.score(xx1,sig1,act1)
##     sc.i.2<- est.score(xx2,sig2,act2)
##     sc.j.1<- est.score(xx1,sig,act)
##     sc.j.2<- est.score(xx2,sig,act)

##     dimf1 <- length(act1)
##     dimf2 <- length(act2)
##     dimf <- dimf1+dimf2
##     dimg <- length(act)
##     if (any(c(dimf1,dimf2,dimg)>min(nrow(xx1),nrow(xx2)))){
##       warning('|active-set| > n')
##       return(list(eval=NA))
##     }else{
##       bfg <- rbind(var(sc.i.1,sc.j.1),var(sc.i.2,sc.j.2))
##       bgf <- t(bfg)
##       bf <- matrix(0,dimf,dimf)
##       bf[1:dimf1,1:dimf1] <- var(sc.i.1)
##       bf[dimf1+(1:dimf2),dimf1+(1:dimf2)] <- var(sc.i.2)
##       bg <- var(sc.j.1)+var(sc.j.2)
      
##       if (dimf>=dimg){
##         mat <- bgf%*%solve(bf)%*%bfg%*%solve(bg)
##         eval<-rep(1,dimf-dimg)
##       }
##       if (dimf<dimg){
##         mat <- bfg%*%solve(bg)%*%bgf%*%solve(bf)
##         eval<-rep(-1,dimg-dimf)
##       }
##       eval.mu.complex<-eigen(mat)$values
##       onemineval.mu <- 1-as.real(eval.mu.complex)
##       onemineval.mu[abs(onemineval.mu)<10^{-6}]<-0
##       eval <- c(eval,sqrt(onemineval.mu),-sqrt(onemineval.mu))
      
##       return(list(ww.mat=mat,eval=eval,eval.mu.complex=eval.mu.complex))
##                                         #return(list(ww.mat=mat,eval=eval,eval.mu.complex=eval.mu.complex,invbg=solve(bg),bfg=bfg))
##     }
##   }

## est.my.ev2 <- function(x1,x2,sig1,sig2,sig,act1,act2,act){

##    ##Estimate Evals ('2nd order' simplification of W)
##    ##
##    ##Input:
##    ##  data: x1,x2
##    ##  sig1,sig2: estimate for individual model
##    ##  sig: estimates for joint model
##    ##  act1,act2: active variables for individual model
##    ##  act: active variables for joint model

##    #dimension of models
##   dimf1 <- length(act1)
##   dimf2 <- length(act2)
##   dimf <- dimf1+dimf2
##   dimg <- length(act)
  
##   #intersection of models
##   ss <- intersect(act,intersect(act1,act2))
  
##   if(length(ss)==0){warning('no intersection between models')}
  
##   aa <- setdiff(act1,ss)
##   bb <- setdiff(act2,ss)
##   cc <- setdiff(act,ss)
##  #print(length(aa));print(length(bb));print(length(cc));print(length(ss))
##   if (any(c(length(aa),length(bb),length(cc),length(ss))>min(nrow(x1),nrow(x2)))){
##     warning('|active-set| > n')
##     return(list(eval=NA))
##     print(length(aa));print(length(bb));print(length(cc));print(length(ss))
##   }else{
##     ev.aux <- ev.aux.complex <- numeric(0)
##     if (dimf>=dimg){
##       if (length(cc)!=0){
##         qcc <- q.matrix2(rbind(x1,x2),sig,sig,act,act,ss)
##         aux.mat <- matrix(0,length(cc),length(cc))
##         if(length(aa)!=0){
##           qac <- q.matrix2(x1,sig1,sig,act1,act,ss)
##           qaa <- q.matrix2(x1,sig1,sig1,act1,act1,ss)
##           aux.mat <- aux.mat+(t(qac)%*%solve(qaa)%*%qac)%*%solve(qcc)
##         }
##         if(length(bb)!=0){
##           qbc <- q.matrix2(x2,sig2,sig,act2,act,ss)
##           qbb <- q.matrix2(x2,sig2,sig2,act2,act2,ss)
##           aux.mat <- aux.mat+(t(qbc)%*%solve(qbb)%*%qbc)%*%solve(qcc)
##         }
##         ev.aux.complex <- eigen(aux.mat)$values
##         ev.aux <- as.real(ev.aux.complex)
##         ev.aux <-  sqrt(1-ev.aux/2)
##       }
##       eval<-rep(1,dimf-dimg)
##       eval <- c(eval,rep(0,2*length(ss)),ev.aux,-ev.aux)
##     }# end if (dimf>=dimg){
  
##     if (dimf<dimg){
##       if (length(cc)!=0){
##         qcc <- q.matrix2(rbind(x1,x2),sig,sig,act,act,ss)
##         if(length(aa)!=0){
##           qac <- q.matrix2(x1,sig1,sig,act1,act,ss)
##           qaa <- q.matrix2(x1,sig1,sig1,act1,act1,ss)
##           oaa <- qac%*%solve(qcc)%*%t(qac)
##           aux.mat.aa<- oaa%*%solve(qaa)
##           aux.mat <- aux.mat.aa
##         }
##         if(length(bb)!=0){
##           qbc <- q.matrix2(x2,sig2,sig,act2,act,ss)
##           qbb <- q.matrix2(x2,sig2,sig2,act2,act2,ss)
##           obb <- qbc%*%solve(qcc)%*%t(qbc)
##           aux.mat.bb<- obb%*%solve(qbb)
##           aux.mat <- aux.mat.bb
##         }
##         if((length(aa)!=0)&(length(bb)!=0)){
##           oab <- qac%*%solve(qcc)%*%t(qbc)
##           aux.mat<- rbind(cbind(aux.mat.aa,oab%*%solve(qbb)),cbind(t(oab)%*%solve(qaa),aux.mat.bb))
##         }
##       }
##       if ((length(aa)!=0)|(length(bb)!=0)){ ##if (length(aa)==0)&(length(bb)==0) 'MI included in MJ; therefore -X^2(dim(MJ)-dim(MI)) distributed'
##         ev.aux.complex <- eigen(aux.mat)$values
##         ev.aux <- as.real(ev.aux.complex)
##         ev.aux <-  sqrt(1-ev.aux/2)
##       }
##       eval<-rep(-1,dimg-dimf)
##       eval <- c(eval,rep(-1,length(ss)),rep(1,length(ss)),rep(0,2*length(ss)),ev.aux,-ev.aux)
##     }
##     return(list(eval=eval,ev.aux.complex=ev.aux.complex))
##   }
   
## }

## twosample_diffnet <- function(x1,x2,b.splits=50,frac.split=1/2,screen.lambda='lambda.cv',gamma.min=0.05,compute.evals='est.my.ev2',diag.invcov=TRUE,...){

##   ##Multisplit Pvalues
##   ##
##   ##Input:
##   ##
##   ## -Data:x1,x2 [scaled data (var=1,mean=0), in particular mu1=mu2=mu=0]
##   ## -b.splits: number of splits
##   ## -frac.split: fraction train-data (for model selection) / test-data (for pval calculations)
##   ## -screen.lambda={lambda.cv,lambda.uni}
##   ## -gamma.min: see Meinshausen&Meier&Buehlmann
##   ## -compute.evals: 'est.my.ev2'/'est.ww.mat2'
##   ## -diag.invcov: TRUE (parameter of interest is invcov; diag(invcov) also considered) / FALSE (param of interest is par.cor; no diagonal)
##   ## -...: additional arguments for cv.glasso(x,...) or glasso.parcor.launi(x,...)

##   ##Important Notes: Pval can be NA, because...
##   ##
##   ##                 1) Eval=sqrt(neg. value)=NA
##   ##                 2) W cannot be estimated (|active-set| > n, problems matrix inversion) 

##   n1 <- nrow(x1)
##   n2 <- nrow(x2)
##   p <- ncol(x1)

##   pval.onesided <- pval.twosided <- rep(NA,length=b.splits)

##   if(screen.lambda=='lambda.cv'){screen.meth <- 'cv.glasso'}
##   if(screen.lambda=='lambda.uni'){screen.meth <- 'glasso.parcor.launi'}

##   if(diag.invcov){df.param <- p*(p+1)/2}
##   if(!diag.invcov){df.param <- p*(p-1)/2}
  
##   for (i in 1:b.splits){
##     split1 <- sample(1:n1,round(n1*frac.split),replace=FALSE)
##     split2 <- sample(1:n2,round(n2*frac.split),replace=FALSE)

##     est.mu <- est.sig <- est.wi <- active <- list()#save est.beta & active variables

##     ##model joint:
##     xx.train <- rbind(x1[split1,],x2[split2,])
##     xx.valid <- rbind(x1[-split1,],x2[-split2,])
##     fit.screen <- eval(as.name(screen.meth))(xx.train,...)
##     act <- which(fit.screen$wi[upper.tri(diag(1,p),diag=TRUE)]!=0)
##     if(!diag.invcov){
##       ind.diag <- which(diag(1,p)[upper.tri(diag(1,p),diag=TRUE)]!=0)
##       act <- setdiff(act,ind.diag)
##     }
##     active[['modJ']] <- act
##     est.mu[['modJ']] <- rep(0,p)
    
##     w <- var(xx.valid)
##     if(!diag.invcov){
##       w <-w*sqrt(tcrossprod(diag(solve(w))))
##     }
##     est.sig[['modJ']] <- w
##     est.wi[['modJ']] <- solve(w)
##     if (length(active[['modJ']])!=df.param){
##       fit.mle <- glasso(var(xx.valid),rho=10^{-10},zero=which(fit.screen$wi==0,arr.in=TRUE))
##       w <- fit.mle$w
##       if(!diag.invcov){
##          w <-w*sqrt(tcrossprod(diag(solve(fit.mle$w))))
##       }
##       est.sig[['modJ']] <- w
##       est.wi[['modJ']] <- fit.mle$wi
##     }
    
##     ##model individual:
##     for (j in c('1','2')){
##       split.train <- eval(as.name(paste('split',j,sep='')))
##       xx.train <- eval(as.name(paste('x',j,sep='')))[split.train,]
##       xx.valid <- eval(as.name(paste('x',j,sep='')))[-split.train,]
##       fit.screen <- eval(as.name(screen.meth))(xx.train,...)
##       act <- which(fit.screen$wi[upper.tri(diag(1,p),diag=TRUE)]!=0)
##       if(!diag.invcov){
##         ind.diag <- which(diag(1,p)[upper.tri(diag(1,p),diag=TRUE)]!=0)
##         act <- setdiff(act,ind.diag)
##       }
##       active[[paste('modIpop',j,sep='')]] <- act
##       est.mu[[paste('modIpop',j,sep='')]] <- rep(0,p)
      
##       w <- var(xx.valid)
##       if(!diag.invcov){
##         w <-w*sqrt(tcrossprod(diag(solve(w))))
##       }
##       est.sig[[paste('modIpop',j,sep='')]] <- w
##       est.wi[[paste('modIpop',j,sep='')]] <- solve(w)
##       if (length(active[[paste('modIpop',j,sep='')]])!=df.param){
##         fit.mle <- glasso(var(xx.valid),rho=10^{-10},zero=which(fit.screen$wi==0,arr.in=TRUE))
##         w <- fit.mle$w
##         if(!diag.invcov){
##           w <-w*sqrt(tcrossprod(diag(solve(fit.mle$w))))
##         }
##         est.sig[[paste('modIpop',j,sep='')]] <- w
##         est.wi[[paste('modIpop',j,sep='')]] <- fit.mle$wi
##       }
      
##     }
    
##     l.act <- lapply(active,length)
##     n1.valid <- nrow(x1[-split1,])
##     n2.valid <- nrow(x2[-split2,])
    
##     if (any(l.act==0)){ warning('at least one active-set is empty')}
##     if (all(l.act>= c(n1.valid+n2.valid,n1.valid,n2.valid))){warning('dim(model) > n-1')}
      
##     teststat <- logratio(x1=x1[-split1,,drop=FALSE],
##                          x2=x2[-split2,,drop=FALSE],
##                          x=rbind(x1[-split1,,drop=FALSE],x2[-split2,,drop=FALSE]),
##                          sig1=est.sig[['modIpop1']],sig2=est.sig[['modIpop2']],sig=est.sig[['modJ']],
##                          mu1=est.mu[['modIpop1']],mu2=est.mu[['modIpop2']],mu=est.mu[['modJ']])
    
      
##     ev.nulldistr <- eval(as.name(compute.evals))(x1[-split1,],x2[-split2,],
##                                                  est.sig[['modIpop1']],est.sig[['modIpop2']],est.sig[['modJ']],
##                                                  active[['modIpop1']],active[['modIpop2']],active[['modJ']])$eval
        
##     if (any(is.na(ev.nulldistr))){
##       warning('Eigenval is NA: pval=NA')
##     }else{
##       pval_onesided <- davies(teststat$twiceLR,lambda=ev.nulldistr)$Qq
##       pval.onesided[i] <- pval_onesided
##       pval.twosided[i] <- 2*min(pval_onesided,1-pval_onesided)
##     }
    
##   }
##   sspval.onesided<-pval.onesided[1]
##   sspval.twosided<-pval.twosided[1]
##   aggpval.onesided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.onesided[!is.na(pval.onesided)])$objective)
##   aggpval.twosided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.twosided[!is.na(pval.twosided)])$objective)
    
##   return(list(pval.onesided=pval.onesided,pval.twosided=pval.twosided,sspval.onesided=sspval.onesided,
##               sspval.twosided=sspval.twosided,
##               aggpval.onesided=aggpval.onesided,aggpval.twosided=aggpval.twosided,
##               LR.last=teststat$twiceLR,active.last=active,sig.last=est.sig,wi.last=est.wi))
## }

## twosample_diffnet_oracle.eval<- function(x1,x2,b.splits=50,frac.split=1/2,screen.lambda='lambda.cv',gamma.min=0.05,infmat,diag.invcov=TRUE,...){

##   ##Multisplit Pvalues
##   ##
##   ##Input:
##   ##
##   ## -Data:x1,x2 [scaled data (var=1,mean=0), in particular mu1=mu2=mu=0]
##   ## -b.splits: number of splits
##   ## -frac.split: fraction train-data (for model selection) / test-data (for pval calculations)
##   ## -screen.lambda={lambda.cv,lambda.uni}
##   ## -gamma.min: see Meinshausen&Meier&Buehlmann
##   ## -infmat: true information-matrix
##   ## -diag.invcov: TRUE (parameter of interest is invcov; diag(invcov) also considered) / FALSE (param of interest is par.cor; no diagonal)
##   ## -...: additional arguments for cv.glasso(x,...) or glasso.parcor.launi(x,...)

##   n1 <- nrow(x1)
##   n2 <- nrow(x2)
##   p <- ncol(x1)

##   pval.onesided <- pval.twosided <- rep(NA,length=b.splits)

##   if(screen.lambda=='lambda.cv'){screen.meth <- 'cv.glasso'}
##   if(screen.lambda=='lambda.uni'){screen.meth <- 'glasso.parcor.launi'}

##   if(diag.invcov){df.param <- p*(p+1)/2}
##   if(!diag.invcov){df.param <- p*(p-1)/2}
  
##   for (i in 1:b.splits){
##     split1 <- sample(1:n1,round(n1*frac.split),replace=FALSE)
##     split2 <- sample(1:n2,round(n2*frac.split),replace=FALSE)

##     est.mu <- est.sig <- est.wi <- active <- list()#save est.beta & active variables

##     ##model joint:
##     xx.train <- rbind(x1[split1,],x2[split2,])
##     xx.valid <- rbind(x1[-split1,],x2[-split2,])
##     fit.screen <- eval(as.name(screen.meth))(xx.train,...)
##     act <- which(fit.screen$wi[upper.tri(diag(1,p),diag=TRUE)]!=0)
##     if(!diag.invcov){
##       ind.diag <- which(diag(1,p)[upper.tri(diag(1,p),diag=TRUE)]!=0)
##       act <- setdiff(act,ind.diag)
##     }
##     active[['modJ']] <- act
##     est.mu[['modJ']] <- rep(0,p)
    
##     w <- var(xx.valid)
##     if(!diag.invcov){
##       w <-w*sqrt(tcrossprod(diag(solve(w))))
##     }
##     est.sig[['modJ']] <- w
##     est.wi[['modJ']] <- solve(w)
##     if (length(active[['modJ']])!=df.param){
##       fit.mle <- glasso(var(xx.valid),rho=10^{-10},zero=which(fit.screen$wi==0,arr.in=TRUE))
##       w <- fit.mle$w
##       if(!diag.invcov){
##         w <-w*sqrt(tcrossprod(diag(solve(fit.mle$w))))#w <-w*sqrt(tcrossprod(diag(fit.mle$wi)))
##       }
##       est.sig[['modJ']] <- w
##       est.wi[['modJ']] <- fit.mle$wi
##     }
    
##     ##model individual:
##     for (j in c('1','2')){
##       split.train <- eval(as.name(paste('split',j,sep='')))
##       xx.train <- eval(as.name(paste('x',j,sep='')))[split.train,]
##       xx.valid <- eval(as.name(paste('x',j,sep='')))[-split.train,]
##       fit.screen <- eval(as.name(screen.meth))(xx.train,...)
##       act <- which(fit.screen$wi[upper.tri(diag(1,p),diag=TRUE)]!=0)
##       if(!diag.invcov){
##         ind.diag <- which(diag(1,p)[upper.tri(diag(1,p),diag=TRUE)]!=0)
##         act <- setdiff(act,ind.diag)
##       }
##       active[[paste('modIpop',j,sep='')]] <- act
##       est.mu[[paste('modIpop',j,sep='')]] <- rep(0,p)
      
##       w <- var(xx.valid)
##       if(!diag.invcov){
##         w <-w*sqrt(tcrossprod(diag(solve(w))))
##       }
##       est.sig[[paste('modIpop',j,sep='')]] <- w
##       est.wi[[paste('modIpop',j,sep='')]] <- solve(w)
##       if (length(active[[paste('modIpop',j,sep='')]])!=df.param){
##         fit.mle <- glasso(var(xx.valid),rho=10^{-10},zero=which(fit.screen$wi==0,arr.in=TRUE))
##         w <- fit.mle$w
##         if(!diag.invcov){
##           w <-w*sqrt(tcrossprod(diag(solve(fit.mle$w))))
##         }
##         est.sig[[paste('modIpop',j,sep='')]] <- w
##         est.wi[[paste('modIpop',j,sep='')]] <- fit.mle$wi
##       }
      
##     }
    
##     l.act <- lapply(active,length)
##     n1.valid <- nrow(x1[-split1,])
##     n2.valid <- nrow(x2[-split2,])
    
##     if (any(l.act==0)){ warning('at least one active-set is empty')}
##     if (all(l.act>= c(n1.valid+n2.valid,n1.valid,n2.valid))){warning('dim(model) > n-1')}
      
##     teststat <- logratio(x1=x1[-split1,,drop=FALSE],
##                          x2=x2[-split2,,drop=FALSE],
##                          x=rbind(x1[-split1,,drop=FALSE],x2[-split2,,drop=FALSE]),
##                          sig1=est.sig[['modIpop1']],sig2=est.sig[['modIpop2']],sig=est.sig[['modJ']],
##                          mu1=est.mu[['modIpop1']],mu2=est.mu[['modIpop2']],mu=est.mu[['modJ']])

    
##     ev.nulldistr <- my.ev2(infmat,act=active[['modJ']],act1=active[['modIpop1']],act2=active[['modIpop2']])$eval 
    
        
##     if (any(is.na(ev.nulldistr))){
##       warning('Eigenval is NA: pval=NA')
##     }else{
##       pval_onesided <- davies(teststat$twiceLR,lambda=ev.nulldistr)$Qq
##       pval.onesided[i] <- pval_onesided
##       pval.twosided[i] <- 2*min(pval_onesided,1-pval_onesided)
##     }
    
##   }
##   sspval.onesided<-pval.onesided[1]
##   sspval.twosided<-pval.twosided[1]
##   aggpval.onesided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.onesided[!is.na(pval.onesided)])$objective)
##   aggpval.twosided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.twosided[!is.na(pval.twosided)])$objective)
    
##   return(list(pval.onesided=pval.onesided,pval.twosided=pval.twosided,sspval.onesided=sspval.onesided,
##               sspval.twosided=sspval.twosided,
##               aggpval.onesided=aggpval.onesided,aggpval.twosided=aggpval.twosided,
##               LR.last=teststat$twiceLR,active.last=active,sig.last=est.sig,wi.last=est.wi))
## }

## twosample_diffnet_oracle.act<- function(x1,x2,b.splits=50,frac.split=1/2,adj.oracle,gamma.min=0.05,compute.evals='est.my.ev2',diag.invcov=TRUE,...){

##   ##Multisplit Pvalues
##   ##
##   ##Input:
##   ##
##   ## -Data:x1,x2 [scaled data (var=1,mean=0), in particular mu1=mu2=mu=0]
##   ## -b.splits: number of splits
##   ## -frac.split: fraction train-data (for model selection) / test-data (for pval calculations)
##   ## -adj.oracle: 
##   ## -gamma.min: see Meinshausen&Meier&Buehlmann
##   ## -compute.evals: 'est.my.ev2'/'est.ww.mat2'
##   ## -diag.invcov: TRUE (parameter of interest is invcov; diag(invcov) also considered) / FALSE (param of interest is par.cor; no diagonal)
##   ## -...: additional arguments for cv.glasso(x,...) or glasso.parcor.launi(x,...)

##   n1 <- nrow(x1)
##   n2 <- nrow(x2)
##   p <- ncol(x1)

##   pval.onesided <- pval.twosided <- rep(NA,length=b.splits)

##   if(screen.lambda=='lambda.cv'){screen.meth <- 'cv.glasso'}
##   if(screen.lambda=='lambda.uni'){screen.meth <- 'glasso.parcor.launi'}

##   if(diag.invcov){df.param <- p*(p+1)/2}
##   if(!diag.invcov){df.param <- p*(p-1)/2}
  
##   for (i in 1:b.splits){
##     split1 <- sample(1:n1,round(n1*frac.split),replace=FALSE)
##     split2 <- sample(1:n2,round(n2*frac.split),replace=FALSE)

##     est.mu <- est.sig <- est.wi <- active <- list()#save est.beta & active variables

##     ##model joint:
##     xx.train <- rbind(x1[split1,],x2[split2,])
##     xx.valid <- rbind(x1[-split1,],x2[-split2,])
##     act <- which(adj.oracle[['modJ']]!=0)
##     if(!diag.invcov){
##       ind.diag <- which(diag(1,p)[upper.tri(diag(1,p),diag=TRUE)]!=0)
##       act <- setdiff(act,ind.diag)
##     }
##     active[['modJ']] <- act
##     est.mu[['modJ']] <- rep(0,p)
    
##     w <- var(xx.valid)
##     if(!diag.invcov){
##       w <-w*sqrt(tcrossprod(diag(solve(w))))
##     }
##     est.sig[['modJ']] <- w
##     est.wi[['modJ']] <- solve(w)
##     if (length(active[['modJ']])!=df.param){
##       fit.mle <- glasso(var(xx.valid),rho=10^{-10},zero=which(adj.oracle[['modJ']]==0,arr.in=TRUE))
##       w <- fit.mle$w
##       if(!diag.invcov){
##         w <-w*sqrt(tcrossprod(diag(fit.mle$wi)))
##       }
##       est.sig[['modJ']] <- w
##       est.wi[['modJ']] <- fit.mle$wi
##     }
    
##     ##model individual:
##     for (j in c('1','2')){
##       split.train <- eval(as.name(paste('split',j,sep='')))
##       xx.train <- eval(as.name(paste('x',j,sep='')))[split.train,]
##       xx.valid <- eval(as.name(paste('x',j,sep='')))[-split.train,]
##       act <- which(adj.oracle[[paste('modIpop',j,sep='')]]!=0)
##       if(!diag.invcov){
##         ind.diag <- which(diag(1,p)[upper.tri(diag(1,p),diag=TRUE)]!=0)
##         act <- setdiff(act,ind.diag)
##       }
##       active[[paste('modIpop',j,sep='')]] <- act
##       est.mu[[paste('modIpop',j,sep='')]] <- rep(0,p)
      
##       w <- var(xx.valid)
##       if(!diag.invcov){
##         w <-w*sqrt(tcrossprod(diag(solve(w))))
##       }
##       est.sig[[paste('modIpop',j,sep='')]] <- w
##       est.wi[[paste('modIpop',j,sep='')]] <- solve(w)
##       if (length(active[[paste('modIpop',j,sep='')]])!=df.param){
##         fit.mle <- glasso(var(xx.valid),rho=10^{-10},zero=which(adj.oracle[[paste('modIpop',j,sep='')]]==0,arr.in=TRUE))
##         w <- fit.mle$w
##         if(!diag.invcov){
##           w <-w*sqrt(tcrossprod(diag(fit.mle$wi)))
##         }
##         est.sig[[paste('modIpop',j,sep='')]] <- w
##         est.wi[[paste('modIpop',j,sep='')]] <- fit.mle$wi
##       }
      
##     }
    
##     l.act <- lapply(active,length)
##     n1.valid <- nrow(x1[-split1,])
##     n2.valid <- nrow(x2[-split2,])
    
##     if (any(l.act==0)){ warning('at least one active-set is empty')}
##     if (all(l.act>= c(n1.valid+n2.valid,n1.valid,n2.valid))){warning('dim(model) > n-1')}
      
##     teststat <- logratio(x1=x1[-split1,,drop=FALSE],
##                          x2=x2[-split2,,drop=FALSE],
##                          x=rbind(x1[-split1,,drop=FALSE],x2[-split2,,drop=FALSE]),
##                          sig1=est.sig[['modIpop1']],sig2=est.sig[['modIpop2']],sig=est.sig[['modJ']],
##                          mu1=est.mu[['modIpop1']],mu2=est.mu[['modIpop2']],mu=est.mu[['modJ']])
    
      
##     ev.nulldistr <- eval(as.name(compute.evals))(x1[-split1,],x2[-split2,],
##                                                  est.sig[['modIpop1']],est.sig[['modIpop2']],est.sig[['modJ']],
##                                                  active[['modIpop1']],active[['modIpop2']],active[['modJ']])$eval
        
##     if (any(is.na(ev.nulldistr))){
##       warning('Eigenval is NA: pval=NA')
##     }else{
##       pval_onesided <- davies(teststat$twiceLR,lambda=ev.nulldistr)$Qq
##       pval.onesided[i] <- pval_onesided
##       pval.twosided[i] <- 2*min(pval_onesided,1-pval_onesided)
##     }
    
##   }
##   sspval.onesided<-pval.onesided[1]
##   sspval.twosided<-pval.twosided[1]
##   aggpval.onesided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.onesided[!is.na(pval.onesided)])$objective)
##   aggpval.twosided <- min(1,(1-log(gamma.min))*optimize(f=agg.pval,interval=c(gamma.min,1),maximum=FALSE,pval=pval.twosided[!is.na(pval.twosided)])$objective)
    
##   return(list(pval.onesided=pval.onesided,pval.twosided=pval.twosided,sspval.onesided=sspval.onesided,
##               sspval.twosided=sspval.twosided,
##               aggpval.onesided=aggpval.onesided,aggpval.twosided=aggpval.twosided,
##               LR.last=teststat$twiceLR,active.last=active,sig.last=est.sig,wi.last=est.wi))
## }

