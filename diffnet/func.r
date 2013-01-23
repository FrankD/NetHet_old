sparse_conc <- function(p,K,s,s.common){
    ##Generate K different Sparse Inverse Covariance-Matrices of dimension p:
    ##
    ##-for each SigInv there are s non-zero entries
    ##-s.common locations of non-zero entries are common among all SigInv;
    ## whereas s-s.common non-zero entries are at different locations

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

    return(SigInv)
}


