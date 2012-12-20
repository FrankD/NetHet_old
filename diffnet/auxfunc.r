sparse_conc <- function(p,K,s,s.common){
    ##Generate K different Sparse Inverse Covariance-Matrices of dimension p:
    ##
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
    B[comp.nonzero] <- 0.9
    B.list[[1]] <- B+t(B)
    if (K>1){
      for (k in 2:K){
        B <- matrix(0,p,p)
        comp.nonzero <- c(same.nonzero,sample(remain.zero,size=s-s.common,replace=FALSE))
        tot.nonzero <- union(tot.nonzero,comp.nonzero)
        B[comp.nonzero] <- 0.5
        B.list[[k]] <- B+t(B)
        remain.zero <- setdiff(ind.upper.tri,tot.nonzero)
      }
    }

    SigInv <- list()
    for (k in 1:K){
      SigInv[[k]] <- B.list[[k]]
      ev <- eigen(SigInv[[k]])$values
      del <- (max(ev)-min(ev)*p)/(p-1)
      diag(SigInv[[k]]) <- del
      SigInv[[k]] <- SigInv[[k]]/del
    }
  }
    
    return(SigInv)
}

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


