
#' DiffNetGroup-package
#'
#' Performs Differential network (DiffNet) analysis for different disease comparsions
#' 
#'
#' 
#' 
#' 
#' @references St\"adler, N. and Mukherjee, S. (2013). Two-Sample Testing in High-Dimensional Models.
#' Preprint \url{http://arxiv.org/abs/1210.4584}.
#' @import multicore DiffNet ICSNP
#' @docType package
#' @name DiffNetGroup-package
#' @useDynLib DiffNetGroup
NULL

###########################
#####Required Packages#####
###########################
library(multicore)
library(DiffNet)
library(ICSNP)

############################
#####P-value adjustment#####
############################
##' pvalue correction
##'
##' 
##' @title pvalue correction
##' @param p vector of pvalues
##' @param method fdr, bonferroni, none
##' @return adjusted pvalues
##' @author n.stadler
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
#############################
#####p-value aggregation#####
#############################
##' pvalue aggregation (Meinshausen & Buehlmann)
##'
##' 
##' @title pvalue aggregation (Meinshausen & Buehlmann)
##' @param pval no decr
##' @param gamma.min no decr
##' @return no decr
##' @author n.stadler
aggpval <- function(pval,gamma.min=0.05){
  
  min(1,(1-log(gamma.min))*optimize(
                                    f=function(gamma){
                                      min(quantile(pval[!is.na(pval)]/gamma,probs=gamma),1)
                                    }
                                    ,interval=c(gamma.min,1),maximum=FALSE)$objective)
}
######################
#####single-split#####
######################
##' diffnet analysis for different disease comparisons for one a single data-split
##'
##'
##' @title diffnet analysis for different disease comparisons for one a single data-split
##' @param x data-matrix (rows: samples; cols: proteins)
##' @param groups vector specifying which sample belongs to which disease (factor)
##' @param method.p.adjust no descr
##' @param ... no descr
##' @return no descr
##' @author n.stadler
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
#####################
#####multi-split#####
#####################
##' diffnet analysis for different disease comparisons for many data-split
##'
##' this is the main function for performing the diffnet analysis of pancan methods paper
##' @title diffnet analysis for different disease comparisons for many data-split
##' @param x data-matrix (rows: samples; cols: proteins); take non-centered, non-scaled data-matrix
##' @param groups vector specifying which sample belongs to which disease (factor)
##' @param no.splits number of data-splits. default=50
##' @param method.p.adjust method for p-value correction (default='bonferroni). paper uses 'fdr'
##' @param order.adj.agg order of p-value adjustment and p-value aggregation (default='adj-agg'). paper uses 'agg-adj'
##' @param ... no descr
##' @return list of output
##' \item{pvalmed}{median aggregated p-values. this are the pvalues reported in the paper (with option 'fdr' and 'agg-adj') }
##' \item{pval}{matrix of pvalues for all sample-splits (not corrected and not aggregated)}
##' \item{teststatmed}{median aggregated test-statistics}
##' \item{teststat}{test-statistics for all sample-splits}
##' @author n.stadler
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
#################################################################################################
#####convert pvalue vector of disease comparisons into matrix of pvalues diseases by diseases####
#################################################################################################
##' convert pvalue vector of disease comparisons into matrix of pvalues (diseases by diseases)
##'
##' 
##' @title convert pvalue vector of disease comparisons into matrix of pvalues (diseases by diseases)
##' @param vec vector of p-values
##' @param groups vector with different diseases (8! this is not the same variable as above; it only
##' contains the levels of the factor groups from above)
##' @return matrix of pvalues diseases by diseases
##' @author n.stadler
convert2mat <- function(vec,groups){
  nr.groups <- length(groups)
  nr.comp <- 0.5*(nr.groups)*(nr.groups-1)
  uptri.rownr <- row(matrix(NA,nr.groups,nr.groups))[upper.tri(matrix(NA,nr.groups,nr.groups),diag=FALSE)]
  uptri.colnr <- col(matrix(NA,nr.groups,nr.groups))[upper.tri(matrix(NA,nr.groups,nr.groups),diag=FALSE)]
  mat <- matrix(0,nr.groups,nr.groups)
  colnames(mat) <- rownames(mat) <- groups
  for (i in 1:nr.comp){
    cl1 <- uptri.rownr[i]
    cl2 <- uptri.colnr[i]
    mat[cl1,cl2] <- vec[i]
  }
  mat <- mat+t(mat)
  return(mat)
}

###########################################################
#####Other functions not used for pancan methods paper#####
###########################################################
trace.mat <- function(m){
  sum(diag(m))
}
##' "SD-test": High-Dim Two-Sample Test for means (Srivastava, 2006)
##'
##' 
##' @title "SD-test": High-Dim Two-Sample Test for means (Srivastava, 2006)
##' @param x1 no descr
##' @param x2 no descr
##' @return no descr
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
##' "SD test" for disease comparison
##'
##' 
##' @title "SD test" for disease comparison
##' @param x  data-matrix (rows: samples; cols: proteins)
##' @param groups vector specifying which sample belongs to which disease (factor)
##' @param method.p.adjust method for p-value adjustment
##' @return list of output
##' @author n.stadler
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
##' HotellingT2-test for disease comparison
##'
##' 
##' @title HotellingT2-test for disease comparison
##' @param x data-matrix (rows: samples; cols: proteins)
##' @param groups vector specifying which sample belongs to which disease (factor)
##' @param method.p.adjust  method for p-value adjustment
##' @return list of output
##' @author n.stadler
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


