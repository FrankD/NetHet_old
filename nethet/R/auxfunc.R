
library(network)

##' Generate two sparse inverse covariance matrices with overlap
##'
##' 
##' @title Generate sparse invcov with overlap
##' @param p number of nodes
##' @param graph 'random' or 'hub'
##' @param K K=1 (null hypothesis) or K=2 (alternative hypothesis)
##' @param n.hub number of hubs (only for graph='hub')
##' @param n.nz number of edges per graph (only for graph='random')
##' @param n.nz.common number of edges incommon between graphs (only for graph='random')
##' @param magn.nz.diff default=0.9
##' @param magn.nz.common default=0.9
##' @param magn.diag default=0
##' @param emin default=0.1 (see ?huge.generator)
##' @export
##' @return 
generate.2networks<- function(p,graph='random',
                              n.nz=rep(p,2),n.nz.common=p,
                              n.hub=2,
                              magn.nz.diff=0.8,
                              magn.nz.common=0.9,magn.diag=0,emin=0.1){

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

        B2 <- matrix(0,p,p)
        nz2 <- c(nz12,sample(remain.zero,size=n.nz[2]-n.nz.common,replace=FALSE))
        B2[nz2] <- magn.nz.common
        B2[setdiff(nz2,nz12)] <- magn.nz.diff
        B.list[[2]] <- B2+t(B2)
    
    }
  
    if(graph=='hub'){
        B.list <- list()
        theta.hub <- as.matrix(huge.generator(d = p, graph = "hub",g=n.hub,verbose=FALSE)$theta)##mk hub-network
        B1 <- B2 <- matrix(0,p,p)
        indhub.diff <- 1:floor(p/n.hub)##identify variables involved in 1st hub
        if(length(indhub.diff)==1){stop('1st hub has only 1 variable; choose smaller number of hubs')}
        B1[-indhub.diff,-indhub.diff] <- theta.hub[-indhub.diff,-indhub.diff]*magn.nz.common
        B.list[[1]] <- B1
    
        B2[-indhub.diff,-indhub.diff] <- theta.hub[-indhub.diff,-indhub.diff]*magn.nz.common
        B2[indhub.diff,indhub.diff] <- theta.hub[indhub.diff,indhub.diff]*magn.nz.diff
        B.list[[2]] <- B2
    
    }
  
  ####compute (positive definite) concentration matrices
    SigInv <- list()
    for (k in 1:2){
        siginv <- B.list[[k]]
        e.min <- min(eigen(siginv)$values)
        siginv <- siginv+diag(abs(e.min)+emin+magn.diag,p)
        SigInv[[k]] <- siginv
        cat('ev.min: ',min(eigen(siginv)$values),'\n')
        cat('condnum: ',max(abs(eigen(siginv)$values))/min(abs(eigen(siginv)$values)),'\n')
    }
    return(list(invcov1=SigInv[[1]],invcov2=SigInv[[2]]))
}

plot.2networks <- function(invcov1,invcov2){
    par(mfrow=c(1,2))
    adj1 <- invcov1!=0
    adj2 <- invcov2!=0
    coord <- plot.network(network(adj1),
                          main="invcov1",
                          vertex.col="orange", 
                          vertex.border=grey(0.8),
                          displaylabels=TRUE,
                          boxed.labels=FALSE,
                          label=paste('X',1:nrow(invcov1),sep=''),
                          label.cex=0.7,
                          usearrows=FALSE, 
                          pad=0.02)
    plot.network(network(adj2),
                 coord=coord,
                 main="invcov2",
                 vertex.col="orange", 
                 vertex.border=grey(0.8),
                 displaylabels=TRUE,
                 boxed.labels=FALSE,
                 label=paste('X',1:nrow(invcov2),sep=''),
                 label.cex=0.7,
                 usearrows=FALSE, 
                 pad=0.02)
}

