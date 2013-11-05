
#' DiffNetGroup-package
#'
#' Performs Differential network (DiffNet) analyses for groups...
#' 
#'
#' DiffNet provides ....
#' 
#' 
#' @references St\"adler, N. and Mukherjee, S. (2013). Two-Sample Testing in High-Dimensional Models.
#' Preprint \url{http://arxiv.org/abs/1210.4584}.
#' @import multicore DiffNet
#' @docType package
#' @name DiffNetGroup-package
#' @useDynLib DiffNetGroup
NULL

#####################
##Required Packages##
#####################
library(multicore)
library(DiffNet)

######################
##P-value adjustment##
######################
my.p.adjust <- function(p,method='fdr'){
  if(method=='fdr'){
    return(p.adjust(p,method='fdr'))
  }
  if(method=='none'){
    return(p.adjust(p,method='none'))
  }
   if(method=='bonferroni'){
    return(p.adjust(p,method='bonferroni'))
  }
  if(method=='qvalue'){
    return(qvalue(p)$qvalues)
  }
}
#######################
##p-value aggregation##
#######################
aggpval <- function(pval,gamma.min=0.05){
  
  min(1,(1-log(gamma.min))*optimize(
                                    f=function(gamma){
                                      min(quantile(pval[!is.na(pval)]/gamma,probs=gamma),1)
                                    }
                                    ,interval=c(gamma.min,1),maximum=FALSE)$objective)
}
################
##single-split##
################
diffnet_groups_singlesplit<- function(x,groups,method.p.adjust="bonferroni",...){
  p <- ncol(x)
  nr.groups <- length(levels(groups))
  nr.comp <- 0.5*(nr.groups)*(nr.groups-1)
  splits.group <- list()
  for(i in levels(groups)){
    n.gr <- sum(groups==i)
    splits.group[[i]] <- sample(1:n.gr,round(n.gr*0.5),replace=FALSE)
  }

  uptri.rownr <- row(matrix(NA,nr.groups,nr.groups))[upper.tri(matrix(NA,nr.groups,nr.groups),diag=FALSE)]
  uptri.colnr <- col(matrix(NA,nr.groups,nr.groups))[upper.tri(matrix(NA,nr.groups,nr.groups),diag=FALSE)]
  
  res <- sapply(1:nr.comp,function(i){
    cl1 <- uptri.rownr[i]
    cl2 <- uptri.colnr[i]
    x1 <- x[which(groups==levels(groups)[cl1]),]
    x2 <- x[which(groups==levels(groups)[cl2]),]
    fit <- diffnet_singlesplit(scale(x1),scale(x2),splits.group[[levels(groups)[cl1]]],splits.group[[levels(groups)[cl2]]],...)
    dfu <- length(fit$active[['modIpop1']])
    dfv<- length(fit$active[['modIpop2']])
    dfuv<- length(fit$active[['modJ']])
    edge.intersect <- (length(intersect(fit$active[['modIpop1']],fit$active[['modIpop2']]))-p)/(length(union(fit$active[['modIpop1']],fit$active[['modIpop2']]))-p)
    teststat.aic <- fit$teststat-2*(dfu+dfv-dfuv)
    pval <- fit$pval.onesided
    ww <- fit$weights.nulldistr
    teststat.aic.sc <- teststat.aic/sqrt(sum(2*(ww)^2))
    cat(paste('comparison: ',levels(groups)[cl1],'-',levels(groups)[cl2],' pval=',pval,sep=''),'\n')
    c(pval,fit$teststat,teststat.aic,teststat.aic.sc,edge.intersect,dfu,dfv,dfuv)
  })
  rownames(res) <- c('pvals','teststat','teststat.aic','teststat.aic.sc','rel.edgeinter','dfu','dfv','dfuv')
  pvals.corrected <- my.p.adjust(res['pvals',],method=method.p.adjust)
  cat('pvals: ',pvals.corrected,'\n')
  return(list(pvals=pvals.corrected,teststat=res['teststat',],teststat.aic=res['teststat.aic',],teststat.aic.sc=res['teststat.aic.sc',],rel.edgeinter=res['rel.edgeinter',],dfu=res['dfu',],dfv=res['dfv',],dfuv=res['dfuv',]))
}
###############
##multi-split##
###############
par.diffnet_groups_multisplit<- function(x,groups,no.splits=50,method.p.adjust='bonferroni',order.adj.agg='adj-agg',...){

  res <- mclapply(seq(no.splits),
                  function(i){
                    cat('split:',i,'\n')
                    res <- diffnet_groups_singlesplit(x,groups,method.p.adjust='none',...)
                    mat <- cbind(res$pvals,res$teststat,res$teststat.aic,res$teststat.aic.sc,res$rel.edgeinter,res$dfu,res$dfv,res$dfuv)
                    colnames(mat) <- c('pvals','teststat','teststat.aic','teststat.aic.sc','rel.edgeinter','dfu','dfv','dfuv')
                    return(mat)
                  }, mc.set.seed=TRUE, mc.preschedule = TRUE)
  res.pval <- sapply(seq(no.splits),function(i){res[[i]][,'pvals']})
  res.teststat <- sapply(seq(no.splits),function(i){res[[i]][,'teststat']})
  res.teststat.aic <- sapply(seq(no.splits),function(i){res[[i]][,'teststat.aic']})
  res.teststat.aic.sc <- sapply(seq(no.splits),function(i){res[[i]][,'teststat.aic.sc']})
  res.rel.edgeinter <- sapply(seq(no.splits),function(i){res[[i]][,'rel.edgeinter']})
  res.dfu <- sapply(seq(no.splits),function(i){res[[i]][,'dfu']})
  res.dfv <- sapply(seq(no.splits),function(i){res[[i]][,'dfv']})
  res.dfuv <- sapply(seq(no.splits),function(i){res[[i]][,'dfuv']})

  if(order.adj.agg=='agg-adj'){
    pvalagg <- my.p.adjust(apply(res.pval,1,aggpval),method=method.p.adjust)
    pvalmed <- my.p.adjust(apply(res.pval,1,median,na.rm=TRUE),method=method.p.adjust)
  }
  if(order.adj.agg=='adj-agg'){
    pvalagg <- apply(apply(res.pval,2,my.p.adjust,method=method.p.adjust),1,aggpval)
    pvalmed <- apply(apply(res.pval,2,my.p.adjust,method=method.p.adjust),1,median)
  }

  return(list(pvalmed=pvalmed,
              pvalagg=pvalagg,
              pval=res.pval,
              teststatmed=apply(res.teststat,1,median,na.rm=TRUE),
              teststatmed.aic=apply(res.teststat.aic,1,median,na.rm=TRUE),
              teststatmed.aic.sc=apply(res.teststat.aic.sc,1,median,na.rm=TRUE),
              teststat=res.teststat,
              teststat.aic=res.teststat.aic,
              teststat.aic.sc=res.teststat.aic.sc,
              rel.edgeintermed=apply(res.rel.edgeinter,1,median,na.rm=TRUE),
              rel.edgeinter=res.rel.edgeinter,
              dfumed=apply(res.dfu,1,median,na.rm=TRUE),
              dfu=res.dfu,
              dfvmed=apply(res.dfv,1,median,na.rm=TRUE),
              dfv=res.dfv,
              dfuvmed=apply(res.dfuv,1,median,na.rm=TRUE),
              dfuv=res.dfuv))
}

convert2mat <- function(vec,groups){
  nr.groups <- length(levels(groups))
  nr.comp <- 0.5*(nr.groups)*(nr.groups-1)
  uptri.rownr <- row(matrix(NA,nr.groups,nr.groups))[upper.tri(matrix(NA,nr.groups,nr.groups),diag=FALSE)]
  uptri.colnr <- col(matrix(NA,nr.groups,nr.groups))[upper.tri(matrix(NA,nr.groups,nr.groups),diag=FALSE)]
  mat <- matrix(0,nr.groups,nr.groups)
  colnames(mat) <- rownames(mat) <- levels(groups)
  for (i in 1:nr.comp){
    cl1 <- uptri.rownr[i]
    cl2 <- uptri.colnr[i]
    mat[cl1,cl2] <- vec[i]
  }
  mat <- mat+t(mat)
  return(mat)
}

test.sd_groups <- function(x,groups,method.p.adjust='bonferroni'){
  nr.groups <- length(levels(groups))
  uptri.rownr <- row(matrix(NA,nr.groups,nr.groups))[upper.tri(matrix(NA,nr.groups,nr.groups),diag=FALSE)]
  uptri.colnr <- col(matrix(NA,nr.groups,nr.groups))[upper.tri(matrix(NA,nr.groups,nr.groups),diag=FALSE)]
  p <- ncol(x)
  res <- sapply(1:(0.5*(nr.groups)*(nr.groups-1)),function(i){
    cl1 <- uptri.rownr[i]
    cl2 <- uptri.colnr[i]
    cat(paste('comparison: ',levels(groups)[cl1],'-',levels(groups)[cl2],sep=''),'\n')
    x1 <- x[which(groups==levels(groups)[cl1]),]
    x2 <- x[which(groups==levels(groups)[cl2]),]
    l2norm <- sum((colMeans(x1)-colMeans(x2))^2)/p
    fit <- test.sd(x1,x2)
    return(c(fit$pval,fit$teststat,l2norm))
    
  })
  rownames(res) <- c('pval','teststat','l2norm')
  pval <- my.p.adjust(res['pval',],method=method.p.adjust)
  res <- list(pval=pval,teststat=res['teststat',],l2norm=res['l2norm',])
  return(res)
}

hotellingsT2_groups <- function(x,groups,method.p.adjust='bonferroni'){
  nr.groups <- length(levels(groups))
  uptri.rownr <- row(matrix(NA,nr.groups,nr.groups))[upper.tri(matrix(NA,nr.groups,nr.groups),diag=FALSE)]
  uptri.colnr <- col(matrix(NA,nr.groups,nr.groups))[upper.tri(matrix(NA,nr.groups,nr.groups),diag=FALSE)]
  p <- ncol(x)
  fit.pval <- lapply(1:(0.5*(nr.groups)*(nr.groups-1)),function(i){
    cl1 <- uptri.rownr[i]
    cl2 <- uptri.colnr[i]
    cat(paste('comparison: ',levels(groups)[cl1],'-',levels(groups)[cl2],sep=''),'\n')
    x1 <- x[which(groups==levels(groups)[cl1]),]
    x2 <- x[which(groups==levels(groups)[cl2]),]
    if((nrow(x1)+nrow(x2))>p){
      fit <- HotellingsT2(x1,x2,mu=NULL,test='f')
      pval <- fit$p.value
    }else{
      pval <- NA
    }
    return(pval)    
  })
  pval <- simplify2array(fit.pval)
  pval <- my.p.adjust(pval,method=method.p.adjust)
  res <- list(pval=pval)
  return(res)
}





## diffnet_cluster <- function(x,class,no.splits=1,screen.meth='screen_bic.glasso',trunc.k=5){
##   nr.class <- length(levels(class))
##   uptri.rownr <- row(matrix(NA,nr.class,nr.class))[upper.tri(matrix(NA,nr.class,nr.class),diag=FALSE)]
##   uptri.colnr <- col(matrix(NA,nr.class,nr.class))[upper.tri(matrix(NA,nr.class,nr.class),diag=FALSE)]
##   p <- ncol(x)
##   fit.pval <- mclapply(1:(0.5*(nr.class)*(nr.class-1)),function(i){
##     cl1 <- uptri.rownr[i]
##     cl2 <- uptri.colnr[i]
##     cat(paste('comparison: ',levels(class)[cl1],'-',levels(class)[cl2],sep=''),'\n')
##     x1 <- x[which(class==levels(class)[cl1]),]
##     x2 <- x[which(class==levels(class)[cl2]),]
##     fit <- diffnet_multisplit(scale(x1),scale(x2),b.splits=no.splits,
##                                    screen.meth=screen.meth,
##                                    trunc.k=trunc.k,verbose=FALSE)
##     pvalmed <- fit$medpval.onesided
##     return(pvalmed)
    
##   },mc.set.seed=TRUE, mc.preschedule = TRUE)
##   my.pval <- simplify2array(fit.pval)
##   mat.pval <- matrix(1,nr.class,nr.class)
##   colnames(mat.pval) <- rownames(mat.pval) <- levels(class)
##   for (i in 1:(0.5*(nr.class)*(nr.class-1))){
##     cl1 <- uptri.rownr[i]
##     cl2 <- uptri.colnr[i]
##     mat.pval[cl1,cl2] <- my.pval[i]
##   }
##   res <- list(pval=mat.pval)
##   return(res)
## }

## test.sd_cluster <- function(x,class){
##   nr.class <- length(levels(class))
##   uptri.rownr <- row(matrix(NA,nr.class,nr.class))[upper.tri(matrix(NA,nr.class,nr.class),diag=FALSE)]
##   uptri.colnr <- col(matrix(NA,nr.class,nr.class))[upper.tri(matrix(NA,nr.class,nr.class),diag=FALSE)]
##   p <- ncol(x)
##   fit.pval <- lapply(1:(0.5*(nr.class)*(nr.class-1)),function(i){
##     cl1 <- uptri.rownr[i]
##     cl2 <- uptri.colnr[i]
##     cat(paste('comparison: ',levels(class)[cl1],'-',levels(class)[cl2],sep=''),'\n')
##     x1 <- x[which(class==levels(class)[cl1]),]
##     x2 <- x[which(class==levels(class)[cl2]),]
##     fit <- test.sd(x1,x2)
##     pvalmed <- fit$pval
##     return(pvalmed)
    
##   })
##   my.pval <- simplify2array(fit.pval)
##   mat.pval <- matrix(1,nr.class,nr.class)
##   colnames(mat.pval) <- rownames(mat.pval) <- levels(class)
##   for (i in 1:(0.5*(nr.class)*(nr.class-1))){
##     cl1 <- uptri.rownr[i]
##     cl2 <- uptri.colnr[i]
##     mat.pval[cl1,cl2] <- my.pval[i]
##   }
##   res <- list(pval=mat.pval)
##   return(res)
## }
