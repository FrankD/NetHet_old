######################################################################
# Multivariate gene-set testing based on graphical models (GGM-GSA)
# -------------------------------------------------------------------
#
#
######################################################################


#####################
##Required packages##
#####################

library(GSA)
library('limma')
library('multtest')
library(ICSNP)
library(parallel)
library(mvtnorm)

######################
##P-value adjustment##
######################

##' P-value adjustment
##'
##' 
##' @title P-value adjustment
##' @param p Vector of p-values.
##' @param method Method for p-value adjustment (default='fdr').
##' @return Vector of adjusted p-values.
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
  #if(method=='qvalue'){
  #  return(qvalue(p)$qvalues)
  #}
}

##################################################
##Classical gene-set testing (Irizarry approach)##
##################################################

##' T-test (equal variances)
##'
##' 
##' @title T-test
##' @param x1 no descr
##' @param x2 no descr
##' @return no descr
##' @author n.stadler
my.ttest <- function(x1,x2){
  ##equal variances
  s <- sqrt((var(x1)*(length(x1)-1)+var(x2)*(length(x2)-1))/(length(x1)+length(x2)-2))
  tstat <- (mean(x1)-mean(x2))/(s*sqrt(1/length(x1)+1/length(x2)))
  return(tstat)
}

##' T-test (unequal variances)
##'
##' 
##' @title T-test
##' @param x1 no descr
##' @param x2 no descr
##' @return no descr
##' @author n.stadler
my.ttest2 <- function(x1,x2){
  ##unequal variances
  s <- sqrt((var(x1)/length(x1))+(var(x2)/length(x2)))
  tstat <- (mean(x1)-mean(x2))/s
  return(tstat)
}

##' Irizarry aggregate score (shift)
##'
##' 
##' @title Irizarry aggregate score (shift)
##' @param ttstat no descr
##' @param geneset no descr
##' @param gene.name no descr
##' @return no descr
##' @author n.stadler
agg.score.iriz.shift <- function(ttstat,geneset,gene.name){
  genes.gs <- which(gene.name%in%geneset)
  ttstat <- ttstat[genes.gs]
  aggscore <- sqrt(length(genes.gs))*mean(ttstat)
  return(aggscore)
}

##' Irizarry aggregate score (scale)
##'
##' 
##' @title Irizarry aggregate score (scale)
##' @param ttstat no descr
##' @param geneset no descr
##' @param gene.name no descr
##' @return no descr
##' @author n.stadler
agg.score.iriz.scale <- function(ttstat,geneset,gene.name){
  genes.gs <- which(gene.name%in%geneset)
  ttstat <- ttstat[genes.gs]
  aggscore <- (sum((ttstat-mean(ttstat))^2)-(length(genes.gs)-1))/(2*(length(genes.gs)-1))
  return(aggscore)
}

##' Irizarry approach (shift only)
##'
##' 
##' @title Irizarry approach (shift only)
##' @param x1 no descr
##' @param x2 no descr
##' @param gene.sets no descr
##' @param gene.names no descr
##' @param gs.names no descr
##' @param method.p.adjust no descr
##' @param alternative no descr
##' @return no descr
##' @author n.stadler
gsea.iriz.shift <- function(x1,x2,gene.sets,gene.names,gs.names=NULL,method.p.adjust='fdr',alternative='two-sided'){
  
  no.genes <- ncol(x1)

  my.tstat <- sapply(1:no.genes,function(i){my.ttest2(x1[,i],x2[,i])})
  agg.scores <- sapply(gene.sets,function(x){agg.score.iriz.shift(my.tstat,x,gene.names)})
  if(alternative=='two-sided'){
    pvals <- my.p.adjust(2*(1-pnorm(abs(agg.scores))),method=method.p.adjust)
  }
  if(alternative=='one-sided'){
    pvals <- my.p.adjust(1-pnorm(agg.scores),method=method.p.adjust)
  }

  if(is.null(gs.names)){
    return(list(tstat=my.tstat,
                agg.scores.sort=agg.scores[order(pvals)],pvals.sort=pvals[order(pvals)],
                agg.scores=agg.scores,pvals=pvals))
  }else{
    return(list(tstat=my.tstat,
                agg.scores.sort=agg.scores[order(pvals)],pvals.sort=pvals[order(pvals)],
                gs.names.sort=gs.names[order(pvals)],
                agg.scores=agg.scores,pvals=pvals))
  }
}

##' Irizarry approach (scale only)
##'
##' 
##' @title Irizarry approach (scale only)
##' @param x1 no descr
##' @param x2 no descr
##' @param gene.sets no descr
##' @param gene.names no descr
##' @param gs.names no descr
##' @param method.p.adjust no descr
##' @param alternative no descr
##' @return no descr
##' @author n.stadler
gsea.iriz.scale <- function(x1,x2,gene.sets,gene.names,gs.names=NULL,method.p.adjust='fdr',alternative='two-sided'){
  
  no.genes <- ncol(x1)

  my.tstat <- sapply(1:no.genes,function(i){my.ttest2(x1[,i],x2[,i])})
  agg.scores <- sapply(gene.sets,function(x){agg.score.iriz.scale(my.tstat,x,gene.names)})
  if(alternative=='two-sided'){
    pvals <- my.p.adjust(2*(1-pnorm(abs(agg.scores))),method=method.p.adjust)
  }
  if(alternative=='one-sided'){
    pvals <- my.p.adjust(1-pnorm(agg.scores),method=method.p.adjust)
  }

  if(is.null(gs.names)){
    return(list(tstat=my.tstat,
                agg.scores.sort=agg.scores[order(pvals)],pvals.sort=pvals[order(pvals)],
                agg.scores=agg.scores,pvals=pvals))
  }else{
    return(list(tstat=my.tstat,
                agg.scores.sort=agg.scores[order(pvals)],pvals.sort=pvals[order(pvals)],
                gs.names.sort=gs.names[order(pvals)],
                agg.scores=agg.scores,pvals=pvals))
  }
}

##' Irizarry approach for gene-set testing
##'
##' Implements the approach described in
##' "Gene set enrichment analysis made simple" by Irizarry et al (2011).
##' It tests for shift and/or change in scale of the distribution.
##' 
##' @title Irizarry approach for gene-set testing
##' @param x1 Expression matrix (condition 1)
##' @param x2 Expression matrix (condition 2)
##' @param gene.sets List of gene-sets
##' @param gene.names Gene names
##' @param gs.names Gene-set names
##' @param method.p.adjust Method for p-value adjustment (default='fdr')
##' @param alternative Default='two-sided' (uses two-sided p-values).
##' @return List consisting of
##' \item{pval.shift}{p-values measuring shift}
##' \item{pval.scale}{p-values measuring scale}
##' \item{pval.combined}{combined p-values (minimum of pval.shift and pval.scale)}
##' @author n.stadler
##' @export
gsea.iriz <- function(x1,x2,gene.sets,gene.names,gs.names=NULL,method.p.adjust='fdr',alternative='two-sided'){

  fit.shift <- gsea.iriz.shift(x1,x2,gene.sets,gene.names,gs.names,method.p.adjust,alternative)
  fit.scale <- gsea.iriz.scale(x1,x2,gene.sets,gene.names,gs.names,method.p.adjust,alternative)
  pvals.shift <- fit.shift$pvals
  agg.scores.shift <- fit.shift$agg.scores
  pvals.scale <- fit.scale$pvals
  agg.scores.scale <- fit.scale$agg.scores

  return(list(pvals.shift=pvals.shift,agg.scores.shift=agg.scores.shift,
              pvals.scale=pvals.scale,agg.scores.scale=agg.scores.scale,
              pvals.combined=pmin(pvals.shift,pvals.scale)))
}


#######################
##Covariance T2-tests##
#######################

##' Classical likelihood-ratio test (equality of covariance matrices)
##'
##' 
##' @title Classical likelihood-ratio test
##' @param x1 no descr
##' @param x2 no descr
##' @param include.mean no descr
##' @return no descr
##' @author n.stadler
t2cov.lr <- function(x1,x2,include.mean=FALSE){
  k <- ncol(x1)
  if(include.mean){
    lr.mle <-logratio(x1,x2,rbind(x1,x2),var(x1),var(x2),var(rbind(x1,x2)),mean(x1),mean(x2),mean(rbind(x1,x2)))$twiceLR
    twosamp.df <- k+k*(k+1)/2
    pval.mle <- 1-pchisq(lr.mle, df=twosamp.df)
  }else{
    lr.mle <-logratio(x1,x2,rbind(x1,x2),var(x1),var(x2),var(rbind(x1,x2)),rep(0,k),rep(0,k),rep(0,k))$twiceLR
    twosamp.df <- k*(k+1)/2
    pval.mle <- 1-pchisq(lr.mle, df=twosamp.df)
  }
  list(pval=pval.mle)
}

##' Diagonal-restricted likelihood-ratio test
##'
##' 
##' @title Diagonal-restricted likelihood-ratio test
##' @param x1 no descr
##' @param x2 no descr
##' @param include.mean no descr
##' @return no descr
##' @author n.stadler
t2diagcov.lr <- function(x1,x2,include.mean=FALSE){
  k <- ncol(x1)
  if(include.mean){
    lr.mle <-logratio(x1,x2,rbind(x1,x2),diag(diag(var(x1))),diag(diag(var(x2))),diag(diag(var(rbind(x1,x2)))),
                      mean(x1),mean(x2),mean(rbind(x1,x2)))$twiceLR
    twosamp.df <- 2*k
    pval.mle <- 1-pchisq(lr.mle, df=twosamp.df)
  }else{
    lr.mle <-logratio(x1,x2,rbind(x1,x2),diag(diag(var(x1))),diag(diag(var(x2))),diag(diag(var(rbind(x1,x2)))),
                      rep(0,k),rep(0,k),rep(0,k))$twiceLR
    twosamp.df <- k
    pval.mle <- 1-pchisq(lr.mle, df=twosamp.df)
  }
  list(pval=pval.mle)
}

##' GSA using T2cov-test
##'
##' 
##' @title GSA using T2cov-test
##' @param x1 expression matrix (condition 1)
##' @param x2 expression matrix (condition 2)
##' @param gene.sets list of gene-sets
##' @param gene.names gene names
##' @param gs.names gene-set names
##' @param method method for testing equality of covariance matrices
##' @param method.p.adjust method for p-value adjustment (default: 'fdr')
##' @return list of results
##' @author n.stadler
gsea.t2cov <- function(x1,x2,gene.sets,gene.names,gs.names=NULL,method='t2cov.lr',method.p.adjust='fdr'){
  
  
  pvals<- sapply(gene.sets,
                 function(y){
                   ind.genes <- which(gene.names%in%y)
                   eval(as.name(method))(x1[,ind.genes],x2[,ind.genes])$pval
                 })
  pvals.corrected <- my.p.adjust(pvals,method=method.p.adjust)
  if(is.null(gs.names)){
    return(list(pvals.sort=pvals.corrected[order(pvals.corrected)],pvals=pvals.corrected))
  }else{
    return(list(pvals.sort=pvals.corrected[order(pvals.corrected)],
                gs.names.sort=gs.names[order(pvals.corrected)],
                pvals=pvals.corrected))
  }
  
}

#################
##Mean T2-tests##
#################

##' HotellingsT2
##'
##' 
##' @title HotellingsT2
##' @param x1 no descr
##' @param x2 no descr
##' @return no descr
##' @author n.stadler
test.t2 <- function(x1,x2){
  fit.t2 <- HotellingsT2(x1,x2, mu = NULL, test = "f")
  list(pval=fit.t2$p.value)
}

##' Trace of a matrix
##'
##' 
##' @title Trace of a matrix
##' @param m no descr
##' @return no descr
##' @author n.stadler
trace.mat <- function(m){
  sum(diag(m))
}

##' High-Dim Two-Sample Test (Srivastava, 2006)
##'
##' 
##' @title High-Dim Two-Sample Test (Srivastava, 2006)
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

##' GSA based on HighdimT2
##'
##' 
##' @title GSA based on HighdimT2
##' @param x1 no descr
##' @param x2 no descr
##' @param gene.sets no descr
##' @param gene.names no descr
##' @param gs.names no descr
##' @param method no descr
##' @param method.p.adjust no descr
##' @return no descr
##' @author n.stadler
gsea.highdimT2 <- function(x1,x2,gene.sets,gene.names,gs.names=NULL,method='test.sd',method.p.adjust='fdr'){
  
  
  pvals<- sapply(gene.sets,
                 function(y){
                   ind.genes <- which(gene.names%in%y)
                   eval(as.name(method))(x1[,ind.genes],x2[,ind.genes])$pval
                 })
  pvals.corrected <- my.p.adjust(pvals,method=method.p.adjust)
  if(is.null(gs.names)){
    return(list(pvals.sort=pvals.corrected[order(pvals.corrected)],pvals=pvals.corrected))
  }else{
    return(list(pvals.sort=pvals.corrected[order(pvals.corrected)],
                gs.names.sort=gs.names[order(pvals.corrected)],
                pvals=pvals.corrected))
  }
  
}

###########################################################
##------------GGMGSA---------------------------------------
###########################################################

##' Filter "non-normal" genes
##'
##' Discarding genes which have Shapiro-Wilk p-value (corrected for multiplicity)
##' smaller than sign.level in either of the two conditions. We used sign.level=0.001
##' in the GGMGSA paper.
##' 
##' @title Filter "non-normal" genes
##' @param x1 expression matrix (condition 1)
##' @param x2 expression matrix (condition 2)
##' @param sign.level sign.level in Shapiro-Wilk tests (default: sign.level=0.001)
##' @return list consisting of
##' \item{x1.filt}{expression matrix (condition 1) after filtering}
##' \item{x2.filt}{expression matrix (condition 2) after filtering}
##' @author n.stadler
shapiro_screen <- function(x1,x2,sign.level=0.001){
  p <- ncol(x1)
  pval1 <- sapply(1:p,function(j){shapiro.test(x1[,j])$p.value})
  pval1.adj <- p.adjust(pval1,method='bonferroni')
  pval2 <- sapply(1:p,function(j){shapiro.test(x2[,j])$p.value})
  pval2.adj <- p.adjust(pval2,method='bonferroni')
  return(list(x1.filt=x1[,-which(pmin(pval1.adj,pval2.adj)<sign.level)],x2.filt=x2[,-which(pmin(pval1.adj,pval2.adj)<sign.level)]))
}

##' Meinshausen p-value aggregation.
##'
##' Inf-quantile formula for p-value aggregation presented in Meinshausen et al 2009.
##' 
##' @title Meinshausen p-value aggregation
##' @param pval Vector of p-values.
##' @param gamma.min See inf-quantile formula of Meinshausen et al 2009 (default=0.05).
##' @return Aggregated p-value.
##' @author n.stadler
##' @export
aggpval <- function(pval,gamma.min=0.05){
  
  min(1,(1-log(gamma.min))*optimize(
                                    f=function(gamma){
                                      min(quantile(pval[!is.na(pval)]/gamma,probs=gamma),1)
                                    }
                                    ,interval=c(gamma.min,1),maximum=FALSE)$objective)
}

##' Single-split GGMGSA
##'
##' 
##' @title Single-split GGMGSA
##' @param x1 centered (scaled) data for condition 1
##' @param x2 centered (scaled) data for condition 2
##' @param gene.sets List of gene-sets.
##' @param gene.names Gene names. Each column in x1 (and x2) corresponds
##'                   to a gene.   
##' @param method.p.adjust Method for p-value adjustment (default='fdr').
##' @param verbose If TRUE, show output progess.
##' @param ... Other arguments (see diffnet_singlesplit).
##' @return List of results.
##' @author n.stadler
ggmgsa_singlesplit <- function(x1,x2,gene.sets,gene.names,method.p.adjust='fdr',verbose=TRUE,...){
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  split1 <- sample(1:n1,round(n1*0.5),replace=FALSE)
  split2 <- sample(1:n2,round(n2*0.5),replace=FALSE)
  res<- sapply(seq(length(gene.sets)),
                 function(i){
                   y <- gene.sets[[i]]
                   if(verbose){cat('  gene set:',i,'\n')}
                   ind.genes <- which(gene.names%in%y)
                   fit <- diffnet_singlesplit(x1[,ind.genes],x2[,ind.genes],split1,split2,...)
                   dfu <- length(fit$active[['modIpop1']])
                   dfv<- length(fit$active[['modIpop2']])
                   dfuv<- length(fit$active[['modJ']])
                   edge.intersect <- length(intersect(fit$active[['modIpop1']],fit$active[['modIpop2']]))/(length(union(fit$active[['modIpop1']],fit$active[['modIpop2']]))-length(ind.genes))
                   teststat.aic <- fit$teststat-2*(dfu+dfv-dfuv)
                   teststat.bic <- fit$teststat-log(0.5*(n1+n2))*(dfu+dfv-dfuv)
                   c(fit$pval.onesided,fit$teststat,teststat.bic,teststat.aic,edge.intersect,dfu,dfv,dfuv)
                 })
  rownames(res) <- c('pvals','teststat','teststat.bic','teststat.aic','rel.edgeinter','dfu','dfv','dfuv')
  pvals.corrected <- my.p.adjust(res['pvals',],method=method.p.adjust)
  return(list(pvals=pvals.corrected,teststat=res['teststat',],teststat.bic=res['teststat.bic',],teststat.aic=res['teststat.aic',],rel.edgeinter=res['rel.edgeinter',],dfu=res['dfu',],dfv=res['dfv',],dfuv=res['dfuv',]))
}

##' Multi-split GGMGSA (parallelized computation)
##'
##' Computation can be parallelized over many data splits.
##' 
##' @title Multi-split GGMGSA (parallelized computation)
##' @param x1 Expression matrix for condition 1 (mean zero is required).
##' @param x2 Expression matrix for condition 2 (mean zero is required).
##' @param no.splits Number of random data splits (default=50).
##' @param gene.sets List of gene-sets.
##' @param gene.names Gene names. Each column in x1 (and x2) corresponds
##'                   to a gene. 
##' @param gs.names Gene-set names (default=NULL).
##' @param method.p.adjust Method for p-value adjustment (default='fdr').
##' @param order.adj.agg Order of aggregation and adjustment of p-values.
##'                      Options: 'agg-adj' (default), 'adj-agg'.
##' @param mc.flag If \code{TRUE} use parallel execution for each no.splits via function 
##'                \code{mclapply} of package \code{parallel}.
##' @param mc.set.seed See mclapply. Default=TRUE
##' @param mc.preschedule See mclapply. Default=TRUE
##' @param verbose If TRUE, show output progess.
##' @param ... Other arguments (see diffnet_singlesplit).
##' @return List consisting of
##' \item{medagg.pval}{Median aggregated p-values}
##' \item{meinshagg.pval}{Meinshausen aggregated p-values}
##' \item{pval}{matrix of p-values before correction and adjustement, dim(pval)=(number of gene-sets)x(number of splits)}
##' \item{teststatmed}{median aggregated test-statistic}
##' \item{teststatmed.bic}{median aggregated bic-corrected test-statistic}
##' \item{teststatmed.aic}{median aggregated aic-corrected test-statistic}
##' \item{teststat}{matrix of test-statistics, dim(teststat)=(number of gene-sets)x(number of splits)}
##' \item{rel.edgeinter}{normalized intersection of edges in condition 1 and 2}
##' \item{df1}{degrees of freedom of GGM obtained from condition 1}
##' \item{df2}{degrees of freedom of GGM obtained from condition 2}
##' \item{df12}{degrees of freedom of GGM obtained from pooled data (condition 1 and 2)}
##' @author n.stadler
##' @export
##' @example ../ggmgsa_ex.R
ggmgsa_multisplit<- function(x1,x2,no.splits=50,gene.sets,gene.names,gs.names=NULL,
                                 method.p.adjust='fdr',order.adj.agg='agg-adj',
                                 mc.flag=FALSE,mc.set.seed=TRUE,mc.preschedule=TRUE,verbose=TRUE,...){

    if(is.null(gs.names)){
        gs.names <- paste('gs',1:length(gene.sets),sep='')
    }

    if(mc.flag==TRUE){
        res <- mclapply(seq(no.splits),
                        function(i){
                            if(verbose){cat('\n split: ',i,'\n\n')}
                            res <- ggmgsa_singlesplit(x1,x2,gene.sets=gene.sets,gene.names=gene.names,
                                                      method.p.adjust='none',verbose=verbose,...)

                            if(verbose){cat('  pvals: ',res$pvals,'\n')}
                            mat <- cbind(res$pvals,res$teststat,res$teststat.bic,res$teststat.aic,res$rel.edgeinter,res$dfu,res$dfv,res$dfuv)
                            colnames(mat) <- c('pvals','teststat','teststat.bic','teststat.aic','rel.edgeinter','dfu','dfv','dfuv')
                            return(mat)
                        }, mc.set.seed=mc.set.seed, mc.preschedule = mc.preschedule)
    }else{

        res <- lapply(seq(no.splits),
                      function(i){
                          if(verbose){cat('\n split: ',i,'\n\n')}
                          res <- ggmgsa_singlesplit(x1,x2,gene.sets=gene.sets,gene.names=gene.names,
                                                    method.p.adjust='none',verbose=verbose,...)
                          if(verbose){cat('  p-values: ',res$pvals,'\n\n')}
                          mat <- cbind(res$pvals,res$teststat,res$teststat.bic,res$teststat.aic,res$rel.edgeinter,res$dfu,res$dfv,res$dfuv)
                          colnames(mat) <- c('pvals','teststat','teststat.bic','teststat.aic','rel.edgeinter','dfu','dfv','dfuv')
                          return(mat)
                      }
                      )
    }
    
    res.pval <- sapply(seq(no.splits),function(i){res[[i]][,'pvals']})
    res.teststat <- sapply(seq(no.splits),function(i){res[[i]][,'teststat']})
    res.teststat.bic <- sapply(seq(no.splits),function(i){res[[i]][,'teststat.bic']})
    res.teststat.aic <- sapply(seq(no.splits),function(i){res[[i]][,'teststat.aic']})
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

    out <- list(medagg.pval=pvalmed,meinshagg.pval=pvalagg,
                pval=res.pval,
                teststatmed=apply(res.teststat,1,median,na.rm=TRUE),
                teststatmed.bic=apply(res.teststat.bic,1,median,na.rm=TRUE),
                teststatmed.aic=apply(res.teststat.aic,1,median,na.rm=TRUE),
                teststat=res.teststat,teststat.aic=res.teststat.aic,
                rel.edgeinter=res.rel.edgeinter,df1=res.dfu,df2=res.dfv,df12=res.dfuv,
                gs.names=gs.names)
    class(out) <- 'ggmgsa'
    return(out)
  
}





