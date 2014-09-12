
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
##' @param verbose If verbose=FALSE then tracing output is disabled.
##' @export
##' @return 
generate.2networks<- function(p,graph='random',
                              n.nz=rep(p,2),n.nz.common=p,
                              n.hub=2,n.hub.diff=1,
                              magn.nz.diff=0.8,
                              magn.nz.common=0.9,magn.diag=0,emin=0.1,
                              verbose=FALSE){

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

        if(n.hub<n.hub.diff){stop("n.hub less than n.hub.diff: choose smaller n.hub.diff")}
        
        ##generate hub network (see library(huge); huge.generator)
        theta.hub <- matrix(0,p,p)
        g.large = p%%n.hub#number of large hubs
        g.small = n.hub - g.large#number of small hubs
        n.small = floor(p/n.hub)#size small hub
        if(n.small<=1){
            stop('hub with less than 2 nodes: choose a smaller n.hub')
        }
        n.large = n.small + 1#size large hub
        g.list = c(rep(n.small, g.small), rep(n.large, g.large))
        g.ind = rep(c(1:n.hub), g.list)
         for (i in 1:n.hub) {
            tmp = which(g.ind == i)
            theta.hub[tmp[1], tmp] = 1
            theta.hub[tmp, tmp[1]] = 1
            rm(tmp)
            gc()
        }
        
        B.list <- list()
        if(n.hub.diff==0){
            B.list[[1]] <- B.list[[2]] <- theta.hub*magn.nz.common
        }else{
            B1 <- B2 <- matrix(0,p,p)
            tmp <- which(g.ind%in%1:n.hub.diff)
            
            B1[-tmp,-tmp] <- theta.hub[-tmp,-tmp]*magn.nz.common
            B.list[[1]] <- B1
    
            B2[-tmp,-tmp] <- theta.hub[-tmp,-tmp]*magn.nz.common
            B2[tmp,tmp] <- theta.hub[tmp,tmp]*magn.nz.diff
            B.list[[2]] <- B2
        }
    
    }
  
  ####compute (positive definite) concentration matrices
    SigInv <- list()
    for (k in 1:2){
        siginv <- B.list[[k]]
        e.min <- min(eigen(siginv)$values)
        siginv <- siginv+diag(abs(e.min)+emin+magn.diag,p)
        SigInv[[k]] <- siginv
        if(verbose){
            cat('ev.min: ',min(eigen(siginv)$values),'\n')
            cat('condnum: ',max(abs(eigen(siginv)$values))/min(abs(eigen(siginv)$values)),'\n')
        }
    }
    return(list(invcov1=SigInv[[1]],invcov2=SigInv[[2]]))
}

##' Plot two networks (GGMs)
##'
##' 
##' @title Plot two networks (GGMs)
##' @param invcov1 Inverse covariance matrix of GGM1.
##' @param invcov2 Inverse covariance matrix of GGM2.
##' @param node.label Names of nodes.
##' @param main Vector (two elements) with network names.
##' @param ... Other arguments (see plot.network).
##' @return Figure with two panels (for each network).
##' @author nicolas
plot.2networks <- function(invcov1,invcov2,
                           node.label=paste('X',1:nrow(invcov1),sep=''),
                           main=c('Network 1','Network 2'),...){
    par(mfrow=c(1,2),mar=c(1,1,1,1))
    adj1 <- invcov1!=0
    adj2 <- invcov2!=0
    adj <- list(adj1,adj2)
    nonzeros <- sapply(adj,sum)
    coord <- plot.network(network(adj[[which.max(nonzeros)]]),
                          main=main[which.max(nonzeros)],
                          displaylabels=TRUE,
                          label=node.label,
                          usearrows=FALSE,...)
    plot.network(network(adj[[which.min(nonzeros)]]),
                 coord=coord,
                 main=main[which.min(nonzeros)],
                 displaylabels=TRUE,
                 label=node.label,
                 usearrows=FALSE,...)
}
