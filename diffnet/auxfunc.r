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

## sparse_mat<- function(p,K,s,s.common,magn.nz=0.5,myseed=1){
##     ##Generate K different Sparse Inverse Covariance-Matrices of dimension p:
##     ##
##     ##-
##     ##-for each SigInv there are s non-zero entries
##     ##-s.common locations of non-zero entries are common among all SigInv;
##     ## whereas s-s.common non-zero entries are at different locations

##   set.seed(myseed)
  
##   if(s==0){
##     SigInv <- list()
##     for (k in 1:K){
##       SigInv[[k]] <- diag(1,p)
##     }
##   }else{
    

##     ind.upper.tri <- which(upper.tri(matrix(NA,p,p)))
    
##     B.list <- list()
##     B <- matrix(0,p,p)
##     comp.nonzero <- tot.nonzero <- sample(ind.upper.tri,size=s,replace=FALSE)
##     same.nonzero <- sample(comp.nonzero,size=s.common,replace=FALSE)
##     remain.zero <- setdiff(ind.upper.tri,comp.nonzero)
##     B[comp.nonzero] <- magn.nz
##     B.list[[1]] <- B+t(B)
##     if (K>1){
##       for (k in 2:K){
##         B <- matrix(0,p,p)
##         comp.nonzero <- c(same.nonzero,sample(remain.zero,size=s-s.common,replace=FALSE))
##         tot.nonzero <- union(tot.nonzero,comp.nonzero)
##         B[comp.nonzero] <- magn.nz
##         B.list[[k]] <- B+t(B)
##         remain.zero <- setdiff(ind.upper.tri,tot.nonzero)
##       }
##     }

##     SigInv <- list()
##     for (k in 1:K){
##       SigInv[[k]] <- B.list[[k]]
##       diag(SigInv[[k]]) <- 1
##     }
##   }
    
##   return(SigInv)
## }

## sparse_conc_v2<- function(p,K,s,s.common,myseed=1){

##   my.magnnz <- 0
##   minev <- 1
##   while(minev>10^-3){
##     my.magnnz <- my.magnnz+0.01
##     spmat<- sparse_mat(p=p,K=K,s=s,s.common=s.common,magn.nz=my.magnnz,myseed=myseed)
##     ev <- sapply(spmat,function(x){min(eigen(x)$values)})
##     minev <- min(ev)
##   }
  
##   magnnz.max <- my.magnnz-0.01
##   set.seed(myseed)
##   spmat<- sparse_mat(p=p,K=K,s=s,s.common=s.common,magn.nz=magnnz.max)
##   ev <- sapply(spmat,function(x){min(eigen(x)$values)});cat('min eigenvalue:',min(ev),'\n');cat('magnnz.max:',magnnz.max,'\n')
##   return(spmat)
## }



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

######################
##Plotting functions##
######################

##False positive -, True positive rate, ROC curve
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

