library(huge)

####################################
#####Generate networks (invcov)#####
####################################

##' Generate networks (invcov) for 2sample testing
##'
##' 
##' @title Generate networks (invcov) for 2sample testing
##' @param p number of nodes
##' @param graph 'random' or 'hub'
##' @param K K=1 (null hypothesis) or K=2 (alternative hypothesis)
##' @param g number of hubs (only for graph='hub')
##' @param n.nz number of edges per graph (only for graph='random')
##' @param n.nz.common number of edges incommon between graphs (only for graph='random')
##' @param magn.nz.diff 
##' @param magn.nz.common 
##' @param magn.diag 
##' @param emin default=0.1 (see ?huge.generator)
##' @return 
##' @author n.stadler
generate.2networks<- function(p,K=2,graph='random',n.nz=rep(p,2),n.nz.common=p,n.hub=2,magn.nz.diff=0.9,magn.nz.common=0.9,magn.diag=0,emin=0.1){

  if(graph=='random'){
    if(min(n.nz)<n.nz.common){stop('min(n.nz) < n.nz.common')}
    ind.upper.tri <- which(upper.tri(matrix(NA,p,p)))
    B.list <- list()
    B1 <- matrix(0,p,p)
    nz1 <- sample(ind.upper.tri,size=n.nz[1],replace=FALSE)
    nz12 <- sample(nz1,size=n.nz.common,replace=FALSE)
    remain.zero <- setdiff(ind.upper.tri,nz1)
    B1[nz12] <- magn.nz.common
    B1[setdiff(nz1,nz12)] <- magn.nz.diff
    B.list[[1]] <- B1+t(B1)

    if(K==2){
      B2 <- matrix(0,p,p)
      nz2 <- c(nz12,sample(remain.zero,size=n.nz[2]-n.nz.common,replace=FALSE))
      B2[nz2] <- magn.nz.common
      B2[setdiff(nz2,nz12)] <- magn.nz.diff
      B.list[[2]] <- B2+t(B2)
    }
  }
  
  if(graph=='hub'){
    B.list <- list()
    theta.hub <- as.matrix(huge.generator(d = p, graph = "hub",g=n.hub)$theta)##mk hub-network
    B1 <- B2 <- matrix(0,p,p)
    indhub.diff <- 1:floor(p/n.hub)##identify variables involved in 1st hub
    if(length(indhub.diff)==1){stop('1st hub has only 1 variable; choose smaller number of hubs')}
    B1[-indhub.diff,-indhub.diff] <- theta.hub[-indhub.diff,-indhub.diff]*magn.nz.common
    B.list[[1]] <- B1
    if(K==2){
      B2[-indhub.diff,-indhub.diff] <- theta.hub[-indhub.diff,-indhub.diff]*magn.nz.common
      B2[indhub.diff,indhub.diff] <- theta.hub[indhub.diff,indhub.diff]*magn.nz.diff
      B.list[[2]] <- B2
    }
  }
  
  ####compute (positive definite) concentration matrices
  SigInv <- list()
  for (k in 1:K){
    siginv <- B.list[[k]]
    e.min <- min(eigen(siginv)$values)
    siginv <- siginv+diag(abs(e.min)+emin+magn.diag,p)
    SigInv[[k]] <- siginv
    cat('ev.min: ',min(eigen(siginv)$values),'\n')
    cat('condnum: ',max(abs(eigen(siginv)$values))/min(abs(eigen(siginv)$values)),'\n')
  }
  return(SigInv)
}
    
##' Generate random concentration matrices with overlap
##'
##' Similar to Rothman 2008; fixed condition number
##' @title Generate random concentration matrices with overlap
##' @param p number of nodes
##' @param K number of concentration matrices
##' @param s number of edges
##' @param s.common number of edges in common 
##' @param magn.nz magnitude of non-zeros
##' @param scale.parcor return parcor ?
##' @param condnum value of concentration matrix 
##' @return SigInv (list)
##' @author n.stadler
sparse_conc <- function(p,K,s,s.common,magn.nz=0.5,scale.parcor=TRUE,condnum=p){
    ##Generate K different Sparse Inverse Covariance-Matrices of dimension p:
    ##
    ##-condition number=condnum (necessary when comparing different performance wrt various p's; if scale.parcor=TRUE, then sparse_conc does not depend on magn.nz)
    ##-for each SigInv there are s non-zero entries
    ##-s.common locations of non-zero entries are common among all SigInv;
    ## whereas s-s.common non-zero entries are at different locations
    
  if(s==0){
    SigInv <- list()
    for (k in 1:K){
      SigInv[[k]] <- diag(1,p)
    }
  }else{
    ind.upper.tri <- which(upper.tri(matrix(NA,p,p)))
    
    B.list <- list()
    B <- matrix(0,p,p)
    comp.nonzero <- tot.nonzero <- sample(ind.upper.tri,size=s,replace=FALSE)
    same.nonzero <- sample(comp.nonzero,size=s.common,replace=FALSE)
    remain.zero <- setdiff(ind.upper.tri,comp.nonzero)
    B[comp.nonzero] <- magn.nz
    B.list[[1]] <- B+t(B)
    if (K>1){
      for (k in 2:K){
        B <- matrix(0,p,p)
        comp.nonzero <- c(same.nonzero,sample(remain.zero,size=s-s.common,replace=FALSE))
        tot.nonzero <- union(tot.nonzero,comp.nonzero)
        B[comp.nonzero] <- magn.nz
        B.list[[k]] <- B+t(B)
        remain.zero <- setdiff(ind.upper.tri,tot.nonzero)
      }
    }

    SigInv <- list()
    for (k in 1:K){
      SigInv[[k]] <- B.list[[k]]
      ev <- eigen(SigInv[[k]])$values
      del <- (max(ev)-min(ev)*condnum)/(condnum-1)
      diag(SigInv[[k]]) <- del
      if(scale.parcor==TRUE){
        SigInv[[k]] <- SigInv[[k]]/del
      }
    }
  }
  return(SigInv)
}

##' Generate random concentration matrices with overlap
##'
##' similar to huge.generator
##' @title Generate random concentration matrices with overlap
##' @param p number of nodes
##' @param K number of siginv's
##' @param s number of edges
##' @param s.common number of common edges
##' @param magn.nz.diff magnitude of different nonzeros
##' @param magn.nz.common  magnitude of common nonzeros
##' @param magn.diag see ?huge.generator
##' @param emin see ?huge.generator
##' @return SigInv (list)
##' @author n.stadler
sparse_conc_v2<- function(p,K,s=p,s.common=p,magn.nz.diff=0.9,magn.nz.common=0.9,magn.diag=0,emin=0.01){

  #####generate adjacency matrices
  ind.upper.tri <- which(upper.tri(matrix(NA,p,p)))
  B.list <- list()
  B <- matrix(0,p,p)
  comp.nonzero <- tot.nonzero <- sample(ind.upper.tri,size=s,replace=FALSE)
  same.nonzero <- sample(comp.nonzero,size=s.common,replace=FALSE)
  remain.zero <- setdiff(ind.upper.tri,comp.nonzero)
  B[same.nonzero] <- magn.nz.common
  B[setdiff(comp.nonzero,same.nonzero)] <- magn.nz.diff
  B.list[[1]] <- B+t(B)
  if (K>1){
    for (k in 2:K){
      B <- matrix(0,p,p)
      comp.nonzero <- c(same.nonzero,sample(remain.zero,size=s-s.common,replace=FALSE))
      tot.nonzero <- union(tot.nonzero,comp.nonzero)
      B[same.nonzero] <- magn.nz.common
      B[setdiff(comp.nonzero,same.nonzero)] <- magn.nz.diff
      B.list[[k]] <- B+t(B)
      remain.zero <- setdiff(ind.upper.tri,tot.nonzero)
    }
  }
  ######compute concentration matrices
  SigInv <- list()
  for (k in 1:K){
    siginv <- B.list[[k]]
    e.min <- min(eigen(siginv)$values)
    siginv <- siginv+diag(abs(e.min)+emin+magn.diag,p)
    SigInv[[k]] <- siginv
    cat('ev.min: ',min(eigen(siginv)$values),'\n')
    cat('condnum: ',max(abs(eigen(siginv)$values))/min(abs(eigen(siginv)$values)),'\n')
  }
  return(SigInv)
}

####################################
#####Other HighdimT2 tests     #####
####################################

##' High-Dim Two-Sample Test based on multiple correction
##'
##' perform univariate t-tests and correct for multiple comparison
##' @title 
##' @param x1 
##' @param x2 
##' @param correction.meth 
##' @param sign.level 
##' @return 
##' @author n.stadler
test.multcor<- function(x1,x2,correction.meth='fdr',sign.level=0.05){
  p <- ncol(x1)
  pval <- sapply(1:p,function(k){t.test(x1[,k],x2[,k],var.equal=TRUE)$p.value})
  pval.adj <- p.adjust(pval, method = correction.meth)
  test.outcome <- as.numeric(!any(pval.adj<sign.level))##if at least one t-test is significant than reject H0 (=0); if no t-test is significant than do not reject (=1)
  return(test.outcome)
}

logratio.varfixed <- function(x1,x2,sig1,sig2,mu1,mu2,mu){
  twiceLR <- 2*(sum(dmvnorm(x1,mean=mu1,sigma=sig1,log=TRUE))+sum(dmvnorm(x2,mean=mu2,sigma=sig2,log=TRUE))-(sum(dmvnorm(x1,mean=mu,sigma=sig1,log=TRUE))+sum(dmvnorm(x2,mean=mu,sigma=sig2,log=TRUE))))
  list(twiceLR=twiceLR)
}

trace.mat <- function(m){
  sum(diag(m))
}

##' High-Dim Two-Sample Test (Srivastava, 2006)
##'
##' .. content for \details{} ..
##' @title 
##' @param x1 
##' @param x2 
##' @return 
##' @author n.stadler
test.sd <- function(x1,x2){
  k <- ncol(x1)
  N1 <- nrow(x1)
  N2 <- nrow(x2)
  n <- N1+N2-2
  mu1 <- colMeans(x1)
  mu2 <- colMeans(x2)
  pooled.cov <-(var(x1)*(N1-1)+var(x2)*(N2-1))/n
  d.pooled.cov <- diag(pooled.cov)
  cormat <- diag(1/sqrt(d.pooled.cov))%*%pooled.cov%*%diag(1/sqrt(d.pooled.cov))
  c.const <- 1+trace.mat(cormat%*%cormat)/(k^{3/2})
  teststat <- ((N1*N2)/(N1+N2))*crossprod(mu1-mu2,diag(1/d.pooled.cov))%*%(mu1-mu2)-(n*k/(n-2))
  teststat <- teststat/sqrt(2*(trace.mat(cormat%*%cormat)-(k^2)/n)*c.const)
  
  list(teststat=teststat,pval=1-pnorm(teststat))
}

############################
#####Plotting functions#####
############################

#####False positive -, True positive rate, ROC curve
my.fpr <- function(res,index=1:3,signi=0.05){
  return(mean(res[index]<signi))
}
my.tpr <- function(res,index=1:3,signi=0.05){
  return(mean(res[index]<signi))
}
my.fprtpr <- function(pval,labels){
  thresh <- sort(unique(c(0,(pval),1.1)))
  n <- ncol(pval)
  fp <- tp <- matrix(NA,length(thresh),n)
  for (i in 1:n){
    for (j in 1:length(thresh)){
      fp[j,i] <- my.fpr(pval[,i],which(labels==0),thresh[j])
      tp[j,i] <- my.tpr(pval[,i],which(labels==1),thresh[j])
    }
  }
  return(list(fpr=fp,tpr=tp))
}
plot.roc <- function(pval,labels,sd=FALSE,add=FALSE,...){
  res <- my.fprtpr(pval,labels)

  if(add==FALSE){
    plot(rowMeans(res$fpr),rowMeans(res$tpr),type='l',xlab='FPR',ylab='TPR',ylim=c(0,1),xlim=c(0,1),...)
  }
  if(add==TRUE){
    lines(rowMeans(res$fpr),rowMeans(res$tpr),type='l',xlab='FPR',ylab='TPR',ylim=c(0,1),xlim=c(0,1),...)
  }
}

