#####################################################################################
#####################################################################################
##HMMGLasso; Greedy Backward Pruning
## +different chromosomes, with possible different prob.trans
## +type={gamma}
## +backward pruning; merge/delete 'similar' components
##Revised:  Date:10/10/2012
##Changes: * in -vers2.R the SigInv are wrongly truncated !
##         * agglo_hmmglasso: don't save fit.hmmgl[[i]] for all number of states i
##         * new merge & delete operation
##         * take prob.init=rep(1/K,K) for computing test-loglikelihood using loglik_chr
##         * new name: bwprun_hmmglasso
##         * pen.wi, pen.ci, pen.pc
##         * merge-step was wrong !
##         * if (any(statesize<=min.statesize)){ KEEP 'OLD' PARAMETERS }
##         * do 'full' optimization when deleting/merging
##         * like func_hmmgl_bwpruning-24/04/2012 but: extra function for backward-pruning with only deleting states/backward-pruning with approx merge/del (see func_hmmgl_bwpruning-31032012)
##         * -19062012.r: miniter option
##         * -25062012: f.init (function for initialisation); small change in loglik_chr (if Phi[i,k]=0 for all k then set Phi[i,k]=1 for all k)
##         * -09082012: nstart.kmeans option, term-option in glasso.invcov, glasso.invcor, glasso.parcor
##         * -09082012: re-initialize bwprun if 'state too small'
##         * -14082012: cleaning-up;
##                      change in glasso.parcor [ww <- diag(s);param <- as.vector(diag(s))];
##                      change in loglik_chr: Phi[which(rowSums(Phi==0)==nr.states),] <- 1, instead of Phi[which(rowSums(Phi==0)==2),] <- 1;
##                      reinit.in/reinit.out-option in bwprun
##         * -15082012: (bwprun) output which state(s) are merged/deleted; plot_bwprun function; mer/del-option
##         * -10102012: mstep: if \lambda=0 then use cov(x.compl) and solve(cov(x.compl)) (and do not run glasso)
##                      add: mixglasso functions
##         * -11102012: mstep: if \lambda=Inf then use diag(cov(x.compl)) and solve(diag(cov(x.compl)))
##         * -12102012: mstep: if we consider consider chromosomes independent we have to update prob.init=apply(u[!duplicated(chromosome),,drop=FALSE],2,mean)
##         * -21102012: estep: if p is large then emission probablities get very small (underflow !);
##                             EXPStep.mix, EXPStep.hmm take log(emission) as input and perform estep by "scaling"
##                      bwprun_hmmglasso: take always prob.init.init=rep(1/(k-1),k-1) (do not take prob.init from previous hmmgl-solution)
##         * -22102012: bwprun_hmmglasso: include option equal.prob.trans=FALSE
##                      func.uinit: new options
##         * -04112012: hmmgl_mixgl-22102012.r originates from func_hmmgl_mixgl_diagcov-22102012.r
##         * -05112012: if (equal.prob.trans==FALSE) then it can happen that prob.trans[k,]=='nearly 0' (for the kth row). this happens if the per chromosome statesize (statesize.chr) is zero. this gives problem because of log in forward-backward equation. therefore: if norm.v.chr_sum_t==Inf, (for one k), then perform estep only on other chromosomes.
##
##
##

#####################
##Required Packages##
#####################
library('glasso')
library(mvtnorm)
#library(HiddenMarkov)
#library(abind)
library(mclust)
library(multicore)

simMIX <- function(n,nr.states,mix.prob,Mu,Sig){
  # mix.prob: mixing probabilities
  # Mu: means
  # Sig: covariance matrix
  # n: nr. of observations
  K <- nr.states
  p <- dim(Mu)[1]
  x <- matrix(0,nrow=n,ncol=p)
  s <- sample(1:K,n,replace=TRUE,p=mix.prob)
  for (i in 1:n){
    x[i,] <- rmvnorm(1,mean=Mu[,s[i]],sigma=Sig[,,s[i]])
  }
  list(S=s,X=x)
}
###########################################
##Sum of non-diagonal entries of a matrix##
###########################################
sumoffdiag <- function(m){
  sum(m)-sum(diag(m))
}
##################
##Initialization##
##################
func.uinit <- function(x,nr.states,init='kmeans',my.cl=NULL,nstart.kmeans=1,iter.max.kmeans=10,modelname.hc="EII"){
  ##Initializing EM (provides uinit)
  ##
  ##x: data-matrix
  ##nr.states:
  ##init: method of initialization {cl.init,r.means,random,kmeans,kmeans.hc,hc}
  ##cl: initial clustering (only necessary if init='cl.init')
  ##nstart.kmeans: no random start in kmeans-initialization (see init='kmeans')
  ##iter.max.kmeans: see documentation of kmeans
  ##modelname.hc: model used in hc-initialization (see init='hc')

  nonvirtual <- apply(!is.na(x),1,all)
  x.complete <- x[nonvirtual,]
  n <- nrow(x.complete)

  if (init=='cl.init'){
    cl <- my.cl
    u <- matrix(0.1,n,nr.states)
    u[col(u)==cl] <- 0.9
    u <- u/rowSums(u)
  }
  if (init=='r.means'){
    cl <- rep(NA,n)
    sub.sample <- sample(1:n,size=nr.states,replace=FALSE)
    cl[sub.sample] <- 1:nr.states
    off.sub.sample <- t(x.complete[setdiff(1:n,sub.sample),])
    init.mean <- x.complete[sub.sample,]
    ss_score <- matrix(NA,nrow=n-nr.states,nr.states)
    for (k in 1:nr.states){
      ss_score[,k]<- colSums((off.sub.sample-init.mean[k,])^2)
    }
    off_cl<- apply(ss_score,1,which.min)
    cl[setdiff(1:n,sub.sample)] <- off_cl
    u <- matrix(0.1,n,nr.states)
    u[col(u)==cl] <- 0.9
    u <- u/rowSums(u)
  }
  if(init=='random'){
    cl <- sample(1:nr.states,size=n,replace=TRUE)
    u <- matrix(0.1,n,nr.states)
    u[col(u)==cl] <- 0.9
    u <- u/rowSums(u)
  }
  if (init=='kmeans'){
    fit.kmeans <- kmeans(x.complete,nr.states,nstart=nstart.kmeans,iter.max=iter.max.kmeans)
    cl <- fit.kmeans$cluster
    u <- matrix(0.1,n,nr.states)
    u[col(u)==cl] <- 0.9
    u <- u/rowSums(u)
  }
  if (init=='kmeans.hc'){
    fit.hc <- hc(modelName = 'EII', data = x.complete)
    cl.hc <- as.vector(hclass(fit.hc,G=nr.states))
    init.center <- by(x.complete,cl.hc,colMeans)
    init.center <- do.call(rbind,init.center)
    fit.kmeans <- kmeans(x.complete,centers=init.center,iter.max=iter.max.kmeans)
    cl <- fit.kmeans$cluster
    u <- matrix(0.1,n,nr.states)
    u[col(u)==cl] <- 0.9
    u <- u/rowSums(u)
  }
  if (init=='hc'){
    fit.hc <- hc(modelName = modelname.hc, data = x.complete)
    cl <- as.vector(hclass(fit.hc,G=nr.states))
    u <- matrix(0.1,n,nr.states)
    u[col(u)==cl] <- 0.9
    u <- u/rowSums(u)
  }
  u.all <- matrix(NA,nrow(x),nr.states)
  u.all[nonvirtual,] <- u
  return(list(u=u.all))
}
###############
##KL-Distance##
###############
tr <- function(m){sum(diag(m))}

symmkldist <- function(mu1,mu2,sig1,sig2){
    symmkl <- 0.5*tr((sig1-sig2)%*%(solve(sig2)-solve(sig1)))+0.5*t(mu1-mu2)%*%(solve(sig1)+solve(sig2))%*%(mu1-mu2)
    return(symmkl)
}

w.kldist <- function(Mu,Sig){
    nr.states <- ncol(Mu)
    res<- matrix(NA,nr.states,nr.states)
    for (k in 1:nr.states){
        for (kk in 1:nr.states){
            res[k,kk] <- symmkldist(Mu[,k],Mu[,kk],Sig[,,k],Sig[,,kk])
        }
    }

  list(state.kldist=res,min.state.kldist=min(res[upper.tri(res)]))
}

EXPStep.mix <- function(logphi,mix.prob){
    u <- logphi+matrix(log(mix.prob),nrow(logphi),ncol(logphi),byrow=TRUE)
    max.u <- apply(u,1,max)
    u.scale <- exp(u-max.u)
    u <- u.scale/rowSums(u.scale)
    loglik <- sum(log(rowSums(u.scale))+max.u)
    return(list(u=u,loglik=loglik))
}
  
##################
##GLasso and Co.##
##################
glasso.parcor <- function(s,rho,penalize.diagonal,maxiter=1000,term=10^{-3}){
    ww <- diag(s)
    iter <- 0
    err <- Inf #convergence of parameters
    param <- as.vector(diag(s))

    while((err>term)&(iter<maxiter)){
        gl <- glasso(s,rho=rho*ww,penalize.diagonal=penalize.diagonal)
        ww <- 1/(diag(gl$wi))
        param.old <- param
        param <- as.vector(gl$w)
        err <- max(abs(param-param.old)/(1+abs(param)))
        iter <- iter+1
    }
    list(w=gl$w,wi=gl$wi,iter=iter)
}
glasso.invcor <- function(s,rho,penalize.diagonal,term=10^{-3}){
  if(penalize.diagonal==FALSE){
    ww <- diag(s)
    gl <- glasso(s,rho=rho*ww,penalize.diagonal=penalize.diagonal)
    return(list(w=gl$w,wi=gl$wi))
  }
  if(penalize.diagonal==TRUE){
    ww <- diag(s)/(1-rho)
    gl <- glasso(s,rho=rho*ww,penalize.diagonal=penalize.diagonal)
    return(list(w=gl$w,wi=gl$wi))
  }
}
glasso.invcov <- function(s,rho,penalize.diagonal,term=10^{-3}){
  gl <- glasso(s,rho=rho,penalize.diagonal=penalize.diagonal)
  list(w=gl$w,wi=gl$wi)
}
##########
##M-Step##
##########
MStepGlasso <- function(x,chromosome=NULL,u,v=NULL,lambda,gamma,pen,penalize.diagonal,equal.prob.trans=NULL,term,model='hmm'){
    ##Mstep for HMMGLasso/MixGLasso (optimizes -completeloglik+lambda*pen)
    ##different chromosomes
    ##x: nxp-data
    ##chromosome: a sequence of the form 1,1,1,1,1,...,2,2,2,2,2,...,3,3,3,3,3,... ; set to NULL for model='mixture'
    ##equal.prob.trans: TRUE (=same prob.trans for (independent) chromosomes); FALSE (=different prob.trans for (independent) chromosomes) ; set to NULL for model='mixture'
    ##u[t,j]=P[S_t=j|obs]
    ##v[t,j,j']=P[S_t-1=j,S_t=j'|obs] ; set to NULL for model='mixture'
    ##model={'hmm','mixture'}

  if(model=='hmm'){
    nr.states <- ncol(u)
    nonvirtual <- apply(!is.na(x),1,all)
    x.compl <- x[nonvirtual,]
    n <- dim(x.compl)[1]
    p <- dim(x.compl)[2]

    Mu <- matrix(NA,ncol=nr.states,nrow=p)
    SigInv <- Sig <- array(NA,dim=c(p,p,nr.states))
    prob.init <- matrix(NA,nr.states,length(levels(chromosome)))
    prob.trans <- array(NA,dim=c(nr.states,nr.states,length(levels(chromosome))))
    #exp.cloglik <- rep(NA,nr.states)

    for (l in 1:nr.states){
      pi.states <- mean(u[nonvirtual,l])
      obj <- cov.wt(x.compl, wt = u[nonvirtual,l], cor = FALSE, center =TRUE, method = c("ML"))
      Mu[,l] <- obj$center
      samplecov <- obj$cov
      if ((lambda!=0)&(lambda!=Inf)){
        fit.glasso <- eval(as.name(pen))(samplecov,rho=2*(pi.states^{gamma})*lambda/(sum(u[nonvirtual,l])),
                                         penalize.diagonal=penalize.diagonal,term=term)
        SigInv[,,l] <- fit.glasso$wi
        Sig[,,l] <- fit.glasso$w
      }
      if(lambda==0){
        SigInv[,,l] <- solve(samplecov)
        Sig[,,l] <- samplecov
      }
      if(lambda==Inf){
        SigInv[,,l] <- diag(1/diag(samplecov))
        Sig[,,l] <- diag(diag(samplecov))
      }
      #exp.cloglik[l] <- -(sum(u[nonvirtual,l])/2)*log(det(fit.glasso$w))-0.5*sum(diag(samplecov%*%fit.glasso$wi))#-lambda*(pi.states^{gamma})*sumoffdiag(abs(fit.glasso$wi))
    }
    if (equal.prob.trans==TRUE){
        prob.init[,1:length(levels(chromosome))] <- apply(u[!duplicated(chromosome),,drop=FALSE],2,mean)
        v_sum_t <- colSums(v,dims=1)
        norm.v_sum_t <- 1/rowSums(v_sum_t)
        prob.trans[,,1:length(levels(chromosome))] <- diag(norm.v_sum_t)%*%v_sum_t
    }
    if (equal.prob.trans==FALSE){
        for (chr in 1:length(levels(chromosome))){
            index.chr <- which(chromosome==chr)
            n.chr <- length(index.chr)
            u.chr <- u[chromosome==chr,]
            prob.init[,chr] <- u.chr[1,]
            v.chr <- v[chromosome==chr,,,drop=FALSE]
            v.chr_sum_t <- colSums(v.chr,dims=1)
            norm.v.chr_sum_t <- 1/rowSums(v.chr_sum_t)
            prob.trans[,,chr] <- diag(norm.v.chr_sum_t)%*%v.chr_sum_t
            prob.trans[norm.v.chr_sum_t==Inf,,chr] <- 0 #set prob.trans[k,k']=0 if number of obs. in state k is zero (otherwise we get Inf)
        }
    }
    return(list(prob.init=prob.init,prob.trans=prob.trans,Mu=Mu,Sig=Sig,SigInv=SigInv))#,exp.cloglik=exp.cloglik)
  }
  if(model=='mixture'){
    nr.states <- ncol(u)
    n <- dim(x)[1]
    p <- dim(x)[2]
    Mu <- matrix(0,ncol=nr.states,nrow=p)
    SigInv <- Sig <- array(0,dim=c(p,p,nr.states))
    mix.prob <- rep(NA,nr.states)
    for (l in 1:nr.states){
      pi.states <- mean(u[,l])
      mix.prob[l] <- pi.states
      obj <- cov.wt(x, wt = u[,l], cor = FALSE, center =TRUE, method = c("ML"))
      Mu[,l] <- obj$center
      samplecov <- obj$cov
      if((lambda!=0)&(lambda!=Inf)){
        fit.glasso <- eval(as.name(pen))(samplecov,rho=2*(pi.states^{gamma})*lambda/(sum(u[,l])),
                                         penalize.diagonal=penalize.diagonal)
        SigInv[,,l] <- fit.glasso$wi
        Sig[,,l] <- fit.glasso$w
      }
      if(lambda==0){
        SigInv[,,l] <- solve(samplecov)
        Sig[,,l] <- samplecov
      }
      if(lambda==Inf){
        SigInv[,,l] <- diag(1/diag(samplecov))
        Sig[,,l] <- diag(diag(samplecov))
      }
    }
    return(list(mix.prob=mix.prob,Mu=Mu,Sig=Sig,SigInv=SigInv))
  }
}

mixglasso_init<- function(x,nr.states,lambda,
                          u.init,mix.prob.init,
                          gamma=0.5,pen='glasso.parcor',penalize.diagonal=FALSE,term=10^{-3},miniter=5,maxiter=1000,min.statesize=5,
                          show.trace=FALSE){

  ##MixGLasso (optimizes -loglik+lambda*pen using EM); Requires u.init as initialization
  ##x: nxp-data
  ##nr.states
  ##lambda
  ##u.init, mix.prob.init: initialization
  ##gamma: see penalty
  ##pen: 'glasso','glasso.parcor','glasso.corinv'
  ##penalize.diagonal: see ?glasso !!!! in general: use PENALIZE.DIAGONAL=FALSE;
  ##                                    more detailed: for gamma=0 use PENALIZE.DIAGONAL=FALSE (PENALIZE.DIAGONAL=TRUE does NOT work !);
  ##                                                   for gamma=1 use PENALIZE.DIAGONAL=TRUE (PENALIZE.DIAGONAL=FALSE does also work)
  ##term: see termination of EM
  ##miniter: minimal number of EM iteration before 'stop EM if any(statesize)<min.statesize' applies
  ##maxiter: maximal number of EM iteration
  ##min.statesize: stop EM if any(statesize)<min.statesize; default value: min.statesize=5; minimial statesize possible is 5

  n <- nrow(x)
  p <- ncol(x)
  Mu <- matrix(NA,p,nr.states)
  Sig <- SigInv <- array(NA,dim=c(p,p,nr.states))

  ##If nr.states=1 do Glasso
  if (nr.states==1){
    obj <- cov.wt(x, cor = FALSE, center =TRUE, method = c("ML"))
    Mu <- obj$center
    samplecov <- obj$cov
    if((lambda!=0)&(lambda!=Inf)){
      fit.glasso <- eval(as.name(pen))(samplecov,rho=2*lambda/n,
                                       penalize.diagonal=penalize.diagonal,term=term)
      Sig <- fit.glasso$w
      SigInv <- fit.glasso$wi
    }
    if(lambda==0){
      Sig <- samplecov
      SigInv <- solve(samplecov)
    }
    if(lambda==Inf){
      SigInv[,,l] <- diag(1/diag(samplecov))
      Sig[,,l] <- diag(diag(samplecov))
    }
    loglik <- sum(dmvnorm(x,Mu,Sig,log=TRUE))
    SigInv[abs(SigInv)<10^{-3}] <- 0
    n.zero <- sum(SigInv==0)
    DF <- p+(p*(p+1)/2)-n.zero/2#df.mean+df.inv.covariance
    list(Mu=matrix(Mu,p,1),Sig=array(Sig,dim=c(p,p,1)),SigInv=array(SigInv,dim=c(p,p,1)),
         loglik=loglik,bic=-loglik+log(n)*DF/2,mmdl=-loglik+log(n)*DF/2,warn='NONE')
  }#end if (nr.states==1){
  else{
    ##Initialisation of parameters
    if(any(colSums(u.init)<=min.statesize)){cat('         -mixglasso: n.init_k <= min.statesize','\n')}
    u <- u.init
    for (l in 1:nr.states){
      obj <- cov.wt(x, wt = u[,l], cor = FALSE, center =TRUE, method = c("ML"))
      Mu[,l] <- obj$center
      samplecov <- obj$cov
      pi.states <- mean(u[,l])
      if((lambda!=0)&(lambda!=Inf)){
        fit.glasso <- eval(as.name(pen))(samplecov,rho=2*(pi.states^{gamma})*lambda/(sum(u[,l])),
                                         penalize.diagonal=penalize.diagonal,term=term)
        Sig[,,l] <- fit.glasso$w
        SigInv[,,l] <- fit.glasso$wi
      }
      if(lambda==0){
        Sig[,,l] <- samplecov
        SigInv[,,l] <- solve(samplecov)
      }
      if(lambda==Inf){
        SigInv[,,l] <- diag(1/diag(samplecov))
        Sig[,,l] <- diag(diag(samplecov))
      }
    }
    mix.prob <- mix.prob.init
    
    ##Start EM iteration
    iter <- 0
    err1 <- Inf #convergence of parameters
    param <- as.vector(Sig)
    warn <- 'NONE'

    while((err1>term)&(iter<maxiter)){
      if (show.trace==TRUE){
        cat('         -mixglasso: iter',iter,'\n')
      }
      ##Estep
      logphi <- matrix(NA,n,nr.states)
      for (l in 1:nr.states){
        logphi[,l] <- dmvnorm(x,Mu[,l],Sig[,,l],log=TRUE)
      }
      fit.E <- EXPStep.mix(logphi,mix.prob)
      unew <- fit.E$u
      loglik <- fit.E$loglik
      if (((any(colSums(unew)<=min.statesize))&(iter>miniter))|any(colSums(unew)<=5)){#min.statesize is a tuning-param for stoping EM; EM stops always if statesize<5
        statesize <- colSums(u)
        cat("         -mixglasso: state too small; min(n_k)=",min(colSums(unew)),'\n')
        warn <- "state too small"
        ##compute Bic&Mmdl
        SigInv[abs(SigInv)<10^{-3}] <- 0
        n.zero <- sum(SigInv==0)
        DF <- p*nr.states+(nr.states-1)+nr.states*(p*(p+1)/2)-n.zero/2#df.mean+df.mix.prob+df.inv.covariance
        n.zero.perstate <- apply(SigInv==0,3,sum)
        penmmdl <- p*sum(log(statesize))+log(n)*(nr.states-1)+sum(log(statesize)*((p*(p+1)/2)-n.zero.perstate/2))
  
        return(list(mix.prob=mix.prob,Mu=Mu,Sig=Sig,SigInv=SigInv,iter=iter,
                    loglik=loglik,bic=-loglik+log(n)*DF/2,mmdl=-loglik+penmmdl/2,
                    u=u,state=apply(u,1,which.max),
                    statesize=statesize,pi.states=colMeans(u),
                    warn=warn))
        break
      }#if (((any(colSums(unew)<=min.statesize))&(iter>miniter))|any(colSums(unew)<=5)){
      
      u<-unew#if statesize>min.statesize then continue EM
      ##Mstep
      fit.M <- MStepGlasso(x=x,u=u,lambda=lambda,gamma=gamma,pen=pen,penalize.diagonal=penalize.diagonal,term=term,model='mixture')
      mix.prob <- fit.M$mix.prob
      Mu <- fit.M$Mu
      Sig <- fit.M$Sig
      SigInv <- fit.M$SigInv

      param.old <- param
      param <- as.vector(Sig)
      err1 <- max(abs(param-param.old)/(1+abs(param)))
      iter <- iter+1
    }#end while((err1>term)&(iter<maxiter)){
    if(maxiter==0){
      logphi <- matrix(NA,n,nr.states)
      for (l in 1:nr.states){
        logphi[,l] <- dmvnorm(x,Mu[,l],Sig[,,l],log=TRUE)
      }
      loglik <- EXPStep.mix(logphi,mix.prob)$loglik
    }#end if(maxiter==0){

    ##compute Bic/Mmdl
    SigInv[abs(SigInv)<10^{-3}] <- 0
    n.zero <- sum(SigInv==0)
    statesize <- colSums(u)
    DF <- p*nr.states+(nr.states-1)+nr.states*(p*(p+1)/2)-n.zero/2#df.mean+df.mix.prob+df.inv.covariance
    n.zero.perstate <- apply(SigInv==0,3,sum)
    penmmdl <- p*sum(log(statesize))+log(n)*(nr.states-1)+sum(log(statesize)*((p*(p+1)/2)-n.zero.perstate/2))
    
    list(mix.prob=mix.prob,Mu=Mu,Sig=Sig,SigInv=SigInv,iter=iter,
         loglik=loglik,bic=-loglik+log(n)*DF/2,mmdl=-loglik+penmmdl/2,
         u=u,state=apply(u,1,which.max),
         statesize=statesize,pi.states=colMeans(u),
         warn=warn)
  }#end else { [(nr.states!=1)]
}


#######################
##hmmglasso/mixglasso##
#######################


mixglasso<- function(x,nr.states,lambda=sqrt(2*nrow(x[apply(!is.na(x),1,all),])*log(ncol(x[apply(!is.na(x),1,all),])))/2,pen='glasso.parcor',
                     init='kmeans.hc',my.cl=NULL,modelname.hc="VVV",nstart.kmeans=1,iter.max.kmeans=10,
                     term=10^{-3},min.statesize=5,...){

  ##MixGLasso (optimizes -loglik+lambda*pen using EM)
  
  ##x: nxp-data
  ##nr.states: = no of mixture components
  ##lambda:
  ##init: initialization
  ##modelname.hc:
  ##nstart.kmeans:
  ##pen: 'glasso.invcov','glasso.parcor','glasso.corinv'
  ##term: see termination of EM
  ##min.statesize: stop EM if any(statesize)<min.statesize
  
  nonvirtual <- apply(!is.na(x),1,all)
  x <- x[nonvirtual,]
  n <- nrow(x)
  p <- ncol(x)
  
  u <- mix.prob <- NULL
  if (nr.states>1){
    fit.init.u <- func.uinit(x,nr.states,init=init,my.cl=my.cl,nstart.kmeans=nstart.kmeans,iter.max.kmeans=iter.max.kmeans,modelname.hc=modelname.hc)
    u <- fit.init.u$u
    mix.prob <- rep(1/nr.states,nr.states)
  }
  fit.mixgl <-  mixglasso_init(x=x,nr.states=nr.states,lambda=lambda,
                               pen=pen,u.init=u,mix.prob.init=mix.prob,
                               min.statesize=min.statesize,term=term,...)
  return(fit.mixgl)
}

Mixglasso <- function(x,nr.states,
                      lambda=sqrt(2*nrow(x)*log(ncol(x)))/2,
                      pen='glasso.parcor',
                      init='kmeans.hc',my.cl=NULL,modelname.hc="VVV",nstart.kmeans=1,iter.max.kmeans=10,
                      term=10^{-3},min.statesize=5,
                      save.allfits=TRUE,filename=NULL,
                      mc.set.seed=FALSE, mc.preschedule = FALSE,...){
                  
  res <- mclapply(nr.states,
                  FUN=function(k){
                    fit.mixgl <-mixglasso(x,k,lambda=lambda,pen=pen,
                                          init=init,my.cl=my.cl,modelname.hc=modelname.hc,
                                          nstart.kmeans=nstart.kmeans,iter.max.kmeans=iter.max.kmeans,
                                          term=term,min.statesize=min.statesize,...)
                    if (save.allfits){
                      save(fit.mixgl,file=paste(filename,'_','fit.mixgl_k',k,'.rda',sep=''))
                    }
                    return(list(bic=fit.mixgl$bic,state=fit.mixgl$state,iter=fit.mixgl$iter,warn=fit.mixgl$warn))},
                  mc.set.seed=mc.set.seed, mc.preschedule = mc.preschedule)
  
  res.bic <- sapply(res,function(x){x[['bic']]})
  res.state <- sapply(res,function(x){x[['state']]})
  res.iter <- sapply(res,function(x){x[['iter']]})
  res.warn <- sapply(res,function(x){x[['warn']]})
  return(list(bic=res.bic,state=res.state,iter=res.iter,warn=res.warn))
}



bwprun_mixglasso <- function(x,
                             nr.states.min=1,nr.states.max,
                             lambda=sqrt(nrow(x[apply(!is.na(x),1,all),])*log(ncol(x[apply(!is.na(x),1,all),])))/2,
                             pen='glasso.parcor',selection.crit='mmdl',
                             term=10^{-3},min.statesize=5,
                             init='kmeans.hc',my.cl=NULL,modelname.hc="VVV",nstart.kmeans=1,iter.max.kmeans=10,
                             reinit.out=FALSE,reinit.in=FALSE,mer=TRUE,del=TRUE,...){
  ##Backward Pruning MixGLasso 8!!!!! so far only working for equal.prob.trans=TRUE
  ##
  ##x: nxp-data
  ##nr.states.min:  should be smaller than the optimal number of states, e.g. nr.states.min=1
  ##nr.states.max:  should be considerably larger than the optimal number of states
  ##lambda:
  ##pen: default 'glasso.parcor'
  ##selection.crit: 'mmdl'/'bic'
  ##term: termination of algorithm
  ##min.statesize: default =5
  ##init: initialisation at Kmax
  ##reinit.out: do re-initialization, if statesize<min.statesize, at the start of algorithm ?
  ##reinit.in: do re-initialization, if statesize<min.statesize, at the bwprun-loop level of algorithm ?
  ##mer:
  ##del:

  ##initialization
  nonvirtual <- apply(!is.na(x),1,all)
  x <- x[nonvirtual,]
  n <- nrow(x)
  p <- ncol(x)
  fit.init.u <- func.uinit(x,nr.states.max,init=init,my.cl=my.cl,nstart.kmeans=nstart.kmeans,iter.max.kmeans=iter.max.kmeans,modelname.hc=modelname.hc)
  u <- fit.init.u$u
  mix.prob <- rep(1/nr.states.max,nr.states.max)
  fit.mixgl <-  mixglasso_init(x=x,nr.states=nr.states.max,lambda=lambda,
                               gamma=0.5,pen=pen,
                               u.init=u,mix.prob.init=mix.prob,
                               min.statesize=min.statesize,term=term,...)

  if(reinit.out==TRUE){
    while(fit.mixgl$warn=='state too small'){#if 'statesize<min.statesize' then set Kmax<-Kmax-1 and re-initialize with 'init'
      cat('bwprun at Kmax: state too small, re-initialize & Kmax<-Kmax-1','\n')
      nr.states.max <- nr.states.max-1
      fit.init.u <- func.uinit(x,nr.states.max,init=init,my.cl=my.cl,nstart.kmeans=nstart.kmeans,iter.max.kmeans=iter.max.kmeans,modelname.hc=modelname.hc)
      u <- fit.init.u$u
      mix.prob <- rep(1/nr.states.max,nr.states.max)
      fit.mixgl <-  mixglasso_init(x=x,nr.states=nr.states.max,lambda=lambda,
                                   gamma=0.5,pen=pen,
                                   u.init=u,mix.prob.init=mix.prob,
                                   min.statesize=min.statesize,term=term,...)
    }
  }
  cat('Kmax=',nr.states.max,'\n')

  ##Start with Bwprun
  ##save initialization,state.name,selcrit.mixgl,re_init_in for k=1,...,K_max
  res.init<-list()
  state.name <- list()
  state.name[[paste('K=',nr.states.max,sep='')]] <- as.character(1:nr.states.max)
  selcrit.mixgl <- rep(NA,nr.states.max)
  re_init_in <- rep(FALSE,nr.states.max)
  
  fit.mixgl.del <- fit.mixgl.mer <- list()
  fit.mixgl.del$warn <- fit.mixgl.mer$warn <- FALSE
  selcrit.del <- selcrit.mer <- Inf

  k <- nr.states.max
  while (k>=nr.states.min){
    cat('bwprun: nr states',k,'\n')
    u <- fit.mixgl$u
    mix.prob <- fit.mixgl$mix.prob
    statesize <- fit.mixgl$statesize
    pi.states <- fit.mixgl$pi.states
    selcrit.mixgl[k] <- fit.mixgl[[selection.crit]]
    res.init[[paste('K=',k,sep='')]]<-list(u=u,mix.prob=fit.mixgl$mix.prob)
    
    if(k!=nr.states.min){

      if(del==TRUE){
        ##delete smallest state
        sm.state <- which.min(statesize)
        u.del <- u[,-sm.state,drop=FALSE]
        u.del[which(rowSums(u.del)==0),]<-1/(k-1)
        u.del <- u.del/rowSums(u.del)
        mix.prob.del <- mix.prob[-sm.state]
        mix.prob.del[sum(mix.prob.del)==0]<-1/(k-1)
        fit.mixgl.del <-  mixglasso_init(x,nr.states=k-1,lambda=lambda,
                                         gamma=0.5,pen=pen,
                                         u.init=u.del,mix.prob.init=mix.prob.del,
                                         min.statesize=min.statesize,term=term,...)
        selcrit.del<-fit.mixgl.del[[selection.crit]]
      }
      if(mer==TRUE){
        ##merge closest states
        bwstate.kldist <- w.kldist(fit.mixgl$Mu,fit.mixgl$Sig)
        sel.k <- which(bwstate.kldist$state.kldist==bwstate.kldist$min.state.kldist,arr.ind=TRUE)[1,]#row(matrix(,k,k))[bwstate.kldist$state.kldist==bwstate.kldist$min.state.kldist]
        u.mer <- cbind(u[,-sel.k],rowSums(u[,sel.k]))
        mix.prob.mer <- c(mix.prob[-sel.k],sum(mix.prob[sel.k]))
        fit.mixgl.mer <-  mixglasso_init(x=x,nr.states=k-1,lambda=lambda,
                                         gamma=0.5,pen=pen,
                                         u.init=u.mer,mix.prob.init=mix.prob.mer,
                                         min.statesize=min.statesize,term=term,...)
        selcrit.mer<-fit.mixgl.mer[[selection.crit]]
      }
      
      if((reinit.in)&((fit.mixgl.del$warn=="state too small")|(fit.mixgl.mer$warn=="state too small"))){#if 'statesize<min.statesize' then set K<-K-1 and re-initialize with 'init'
        cat('        state too small, re-initialize at K-1:','\n')
        cat('        warn(mixglasso.del) ',fit.mixgl.del$warn,'\n')
        cat('        warn(mixglasso.mer) ',fit.mixgl.mer$warn,'\n')
        re_init_in[k-1] <- TRUE
        state.name[[paste('K=',k-1,sep='')]] <- as.character(1:(k-1))
        fit.init.u <- func.uinit(x,k-1,init=init,my.cl=my.cl,nstart.kmeans=nstart.kmeans,iter.max.kmeans=iter.max.kmeans,modelname.hc=modelname.hc)
        fit.mixgl <-  mixglasso_init(x=x,nr.states=k-1,lambda=lambda,
                                     gamma=0.5,pen=pen,
                                     u.init=fit.init.u$u,
                                     mix.prob.init=rep(1/(k-1),(k-1)),
                                     min.statesize=min.statesize,term=term,...)
      }else{#if 'no re-initialization'
        if (selcrit.del<selcrit.mer){
          cat('        merge or delete?','delete state ',state.name[[paste('K=',k,sep='')]][sm.state],'\n\n')
          fit.mixgl<-fit.mixgl.del
          state.name[[paste('K=',k-1,sep='')]] <- state.name[[paste('K=',k,sep='')]][-sm.state]
        }else{
          cat('        merge or delete?','merge states ',state.name[[paste('K=',k,sep='')]][sel.k],'\n\n')
          fit.mixgl<-fit.mixgl.mer
          state.name[[paste('K=',k-1,sep='')]] <- c(state.name[[paste('K=',k,sep='')]][-sel.k],paste(state.name[[paste('K=',k,sep='')]][sel.k],collapse='+'))
        }
      }
    }#end if(k!=nr.states.min){
    k <- k-1
  }
  #optimal selcrit-solution
  fit.mixgl.selcrit <-  mixglasso_init(x,nr.states=which.min(selcrit.mixgl),lambda=lambda,
                                       gamma=0.5,pen=pen,
                                       u.init=res.init[[paste('K=',which.min(selcrit.mixgl),sep='')]]$u,
                                       mix.prob.init=res.init[[paste('K=',which.min(selcrit.mixgl),sep='')]]$mix.prob,
                                       min.statesize=min.statesize,term=term,...)
  
  return(list(selcrit=selcrit.mixgl,res.init=res.init,state.name=state.name,re.init.in=re_init_in,fit.mixgl.selcrit=fit.mixgl.selcrit))
}

################################
##Log-likelihood (hmm/mixture)##
################################

loglik_mix<- function(x,mix.prob,Mu,Sig){
  ##computes loglikelihood of model
  nr.states <- ncol(Mu)
  n <- nrow(x)
  if(nr.states==1){
    return(sum(dmvnorm(x,Mu[,1],Sig[,,1],log=TRUE)))
  }
  else{
    logphi <- matrix(NA,n,nr.states)
    for (l in 1:nr.states){
      logphi[,l] <- dmvnorm(x,Mu[,l],Sig[,,l],log=TRUE)
    }
    loglik <- EXPStep.mix(logphi,mix.prob)$loglik
    return(loglik)
  }
}

