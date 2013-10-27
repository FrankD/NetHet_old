#####################
##Required packages##
#####################

library(GSA)
library('limma')
library('multtest')
#library('genefilter')
library(ICSNP)
library(multicore)
library(DiffNet)
#library(qvalue)

#############################################################
##Some functions for data-format conversion (GSEA specific)##
#############################################################

GSEA.Gct2Frame <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in GCT format and converts it into an R data frame
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
   ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T)
   ds <- ds[-1]
   return(ds)
}

GSEA.Gct2Frame2 <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in GCT format and converts it into an R data frame
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
      content <- readLines(filename)
      content <- content[-1]
      content <- content[-1]
      col.names <- noquote(unlist(strsplit(content[1], "\t")))
      col.names <- col.names[c(-1, -2)]
      num.cols <- length(col.names)
      content <- content[-1]
      num.lines <- length(content)


      row.nam <- vector(length=num.lines, mode="character")
      row.des <- vector(length=num.lines, mode="character")
      m <- matrix(0, nrow=num.lines, ncol=num.cols)

      for (i in 1:num.lines) {
         line.list <- noquote(unlist(strsplit(content[i], "\t")))
         row.nam[i] <- noquote(line.list[1])
         row.des[i] <- noquote(line.list[2])
         line.list <- line.list[c(-1, -2)]
         for (j in 1:length(line.list)) {
            m[i, j] <- as.numeric(line.list[j])
         }
      }
      ds <- data.frame(m)
      names(ds) <- col.names
      row.names(ds) <- row.nam
      return(ds)
}

GSEA.Res2Frame <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in RES format and converts it into an R data frame
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

   header.cont <- readLines(filename, n = 1)
   temp <- unlist(strsplit(header.cont, "\t"))
   colst <- length(temp)
   header.labels <- temp[seq(3, colst, 2)]
   ds <- read.delim(filename, header=F, row.names = 2, sep="\t", skip=3, blank.lines.skip=T, comment.char="", as.is=T)
   colst <- length(ds[1,])
   cols <- (colst - 1)/2
   rows <- length(ds[,1])
   A <- matrix(nrow=rows - 1, ncol=cols)
   A <- ds[1:rows, seq(2, colst, 2)]
   table1 <- data.frame(A)
   names(table1) <- header.labels
   return(table1)
}

GSEA.ReadClsFile <- function(file = "NULL") { 
#
# Reads a class vector CLS file and defines phenotype and class labels vectors for the samples in a gene expression file (RES or GCT format)
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

      cls.cont <- readLines(file)
      num.lines <- length(cls.cont)
      class.list <- unlist(strsplit(cls.cont[[3]], " "))
      s <- length(class.list)
      t <- table(class.list)
      l <- length(t)
      phen <- vector(length=l, mode="character")
      phen.label <- vector(length=l, mode="numeric")
      class.v <- vector(length=s, mode="numeric")
      for (i in 1:l) {
         phen[i] <- noquote(names(t)[i])
         phen.label[i] <- i - 1
      }
      for (i in 1:s) {
         for (j in 1:l) {
             if (class.list[i] == phen[j]) {
                class.v[i] <- phen.label[j]
             }
         }
      }
      return(list(phen = phen, class.v = class.v))
}

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

##################################################
##Classical gene-set testing (Irizarry approach)##
##################################################

my.ttest <- function(x1,x2){
  ##equal variances
  s <- sqrt((var(x1)*(length(x1)-1)+var(x2)*(length(x2)-1))/(length(x1)+length(x2)-2))
  tstat <- (mean(x1)-mean(x2))/(s*sqrt(1/length(x1)+1/length(x2)))
  return(tstat)
}

my.ttest2 <- function(x1,x2){
  ##unequal variances
  s <- sqrt((var(x1)/length(x1))+(var(x2)/length(x2)))
  tstat <- (mean(x1)-mean(x2))/s
  return(tstat)
}

agg.score.iriz.shift <- function(ttstat,geneset,gene.name){
  genes.gs <- which(gene.name%in%geneset)
  ttstat <- ttstat[genes.gs]
  aggscore <- sqrt(length(genes.gs))*mean(ttstat)
  return(aggscore)
}
agg.score.iriz.scale <- function(ttstat,geneset,gene.name){
  genes.gs <- which(gene.name%in%geneset)
  ttstat <- ttstat[genes.gs]
  aggscore <- (sum((ttstat-mean(ttstat))^2)-(length(genes.gs)-1))/(2*(length(genes.gs)-1))
  return(aggscore)
}

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
##' @param x1 
##' @param x2 
##' @param include.mean 
##' @return 
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
##' @param x1 
##' @param x2 
##' @param include.mean 
##' @return 
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

test.t2 <- function(x1,x2){
  fit.t2 <- HotellingsT2(x1,x2, mu = NULL, test = "f")
  list(pval=fit.t2$p.value)
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

##############
##NetGeneSet##
##############
shapiro_screen <- function(x1,x2,sign.level=0.01){
  p <- ncol(x1)
  pval1 <- sapply(1:p,function(j){shapiro.test(x1[,j])$p.value})
  pval1.adj <- p.adjust(pval1,method='bonferroni')
  pval2 <- sapply(1:p,function(j){shapiro.test(x2[,j])$p.value})
  pval2.adj <- p.adjust(pval2,method='bonferroni')
  return(list(x1.filt=x1[,-which(pmin(pval1.adj,pval2.adj)<sign.level)],x2.filt=x2[,-which(pmin(pval1.adj,pval2.adj)<sign.level)]))
}

aggpval <- function(pval,gamma.min=0.05){
  
  min(1,(1-log(gamma.min))*optimize(
                                    f=function(gamma){
                                      min(quantile(pval[!is.na(pval)]/gamma,probs=gamma),1)
                                    }
                                    ,interval=c(gamma.min,1),maximum=FALSE)$objective)
}

gsea.diffnet.singlesplit <- function(x1,x2,gene.sets,gene.names,method.p.adjust='fdr',...){
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  split1 <- sample(1:n1,round(n1*0.5),replace=FALSE)
  split2 <- sample(1:n2,round(n2*0.5),replace=FALSE)
  res<- sapply(seq(length(gene.sets)),
                 function(i){
                   y <- gene.sets[[i]]
                   cat('gene set:',i,'\n')
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


gsea.diffnet.multisplit <- function(x1,x2,no.splits=50,gene.sets,gene.names,gs.names=NULL,method.p.adjust='fdr',order.adj.agg='adj-agg',...){

  res <- lapply(seq(no.splits),
                function(i){
                  cat('split:',i,'\n')
                  res <- gsea.diffnet.singlesplit(x1,x2,gene.sets=gene.sets,gene.names=gene.names,
                                                    method.p.adjust='none',...)
                  cat('pvals: ',res$pvals,'\n')
                  mat <- cbind(res$pvals,res$teststat,res$teststat.bic,res$teststat.aic,res$rel.edgeinter,res$dfu,res$dfv,res$dfuv)
                  colnames(mat) <- c('pvals','teststat','teststat.bic','teststat.aic','rel.edgeinter','dfu','dfv','dfuv')
                  return(mat)
                }
                )
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

  if(is.null(gs.names)){
    return(list(pvalmed=pvalmed,pvalagg=pvalagg,
                pval=res.pval,
                teststatmed=apply(res.teststat,1,median,na.rm=TRUE),
                teststatmed.bic=apply(res.teststat.bic,1,median,na.rm=TRUE),
                teststatmed.aic=apply(res.teststat.aic,1,median,na.rm=TRUE),
                teststat=res.teststat,teststat.aic=res.teststat.aic,rel.edgeinter=res.rel.edgeinter,dfu=res.dfu,dfv=res.dfv,dfuv=res.dfuv))
  }else{
    return(list(pvalmed=pvalmed,pvalagg=pvalagg,
                pval=res.pval,
                teststatmed=apply(res.teststat,1,median,na.rm=TRUE),
                teststatmed.bic=apply(res.teststat.bic,1,median,na.rm=TRUE),
                teststatmed.aic=apply(res.teststat.aic,1,median,na.rm=TRUE),
                teststat=res.teststat,teststat.aic=res.teststat.aic,rel.edgeinter=res.rel.edgeinter,dfu=res.dfu,dfv=res.dfv,dfuv=res.dfuv,
                gs.names=gs.names))
  }
}

par.gsea.diffnet.multisplit <- function(x1,x2,no.splits=50,gene.sets,gene.names,gs.names=NULL,method.p.adjust='fdr',order.adj.agg='adj-agg',...){

  res <- mclapply(seq(no.splits),
                  function(i){
                    cat('split:',i,'\n')
                    res <- gsea.diffnet.singlesplit(x1,x2,gene.sets=gene.sets,gene.names=gene.names,
                                                    method.p.adjust='none',...)

                    cat('pvals: ',res$pvals,'\n')
                    mat <- cbind(res$pvals,res$teststat,res$teststat.bic,res$teststat.aic,res$rel.edgeinter,res$dfu,res$dfv,res$dfuv)
                    colnames(mat) <- c('pvals','teststat','teststat.bic','teststat.aic','rel.edgeinter','dfu','dfv','dfuv')
                    return(mat)
                  }, mc.set.seed=TRUE, mc.preschedule = TRUE)
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

  if(is.null(gs.names)){
    return(list(pvalmed=pvalmed,pvalagg=pvalagg,
                pval=res.pval,
                teststatmed=apply(res.teststat,1,median,na.rm=TRUE),
                teststatmed.bic=apply(res.teststat.bic,1,median,na.rm=TRUE),
                teststatmed.aic=apply(res.teststat.aic,1,median,na.rm=TRUE),
                teststat=res.teststat,teststat.aic=res.teststat.aic,rel.edgeinter=res.rel.edgeinter,dfu=res.dfu,dfv=res.dfv,dfuv=res.dfuv))
  }else{
    return(list(pvalmed=pvalmed,pvalagg=pvalagg,
                pval=res.pval,
                teststatmed=apply(res.teststat,1,median,na.rm=TRUE),
                teststatmed.bic=apply(res.teststat.bic,1,median,na.rm=TRUE),
                teststatmed.aic=apply(res.teststat.aic,1,median,na.rm=TRUE),
                teststat=res.teststat,teststat.aic=res.teststat.aic,rel.edgeinter=res.rel.edgeinter,dfu=res.dfu,dfv=res.dfv,dfuv=res.dfuv,
                gs.names=gs.names))
  }
}

gsea.diffregr.singlesplit <- function(y1,y2,x1,x2,gene.sets,gene.names,method.p.adjust='fdr',screen.meth='lasso.cvtrunc'){
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  split1 <- sample(1:n1,round(n1*0.5),replace=FALSE)
  split2 <- sample(1:n2,round(n2*0.5),replace=FALSE)
  pvals<- sapply(seq(length(gene.sets)),
                 function(i){
                   y <- gene.sets[[i]]
                   cat('gene set:',i,'\n')
                   ind.genes <- which(gene.names%in%y)
                   diffregr_singlesplit(y1,y2,x1[,ind.genes],x2[,ind.genes],split1,split2,screen.meth=screen.meth)$pval.onesided
                 })
  pvals.corrected <- my.p.adjust(pvals,method=method.p.adjust)
  return(list(pvals=pvals.corrected))
}

#########################################################################################################################################
#########################################################################################################################################
#########################################################################################################################################
######################################################################################################################################### 

## gsea.iriz.new <- function(x,pheno,gene.sets,gs.names=NULL,gene.names,method.p.adjust='fdr',alternative='two-sided'){
  
##   no.genes <- ncol(x)
##   my.tstat <- rowttests(as.matrix(x),pheno)$stat
  
##   agg.scores <- sapply(gs.names,function(gs){
##     geneset <- gene.sets[which(gs.names==gs),]
##     geneset <- geneset[geneset!="null"]
##     Index2 <- match(geneset,gene.names)
##     x=my.tstat[Index2]
##     N <- length(x)
##     return(mean(x)*sqrt(N))
##   })
##   if(alternative=='two-sided'){
##     pvals <- my.p.adjust(2*(1-pnorm(abs(agg.scores))),method=method.p.adjust)
##   }
##   if(alternative=='one-sided'){
##     pvals <- my.p.adjust(1-pnorm(agg.scores),method=method.p.adjust)
##   }
##   if(is.null(gs.names)){
##     return(list(tstat=my.tstat,
##                 agg.scores.sort=agg.scores[order(pvals)],pvals.sort=pval[order(pvals)],
##                 agg.scores=agg.scores,pvals=pvals))
##   }else{
##     return(list(tstat=my.tstat,
##                 agg.scores.sort=agg.scores[order(pvals)],pvals.sort=pvals[order(pvals)],
##                 gs.names.sort=gs.names[order(pvals)],
##                 agg.scores=agg.scores,pvals=pvals))
##   }
## }

## rho.max <- function(s){
##   diag(s) <- 0
##   return(max(abs(s)))
## }

## bic.glasso.invcor <- function(x,lambda,penalize.diagonal=FALSE,plot.it=TRUE)
## {
##   ##glasso; lambda.opt with bic
##   aic.score <-rep(NA,length(lambda))
##   Mu <- colMeans(x)
##   samplecov <- var(x)

##   if(is.null(lambda)){
##     la <- lambda
##   }else{
##     la <- 2*lambda/nrow(x)
##   }

##   v.mat <- 1/sqrt(diag(samplecov))
##   v.tcross <- tcrossprod(v.mat)
##   vinv.tcross <- tcrossprod(1/v.mat)
##   fit.path <- glassopath(samplecov*v.tcross,rholist=la,penalize.diagonal=penalize.diagonal,trace=0)

##   loglik <- lapply(seq(length(fit.path$rholist)),
##                    function(i){
##                      return(sum(dmvnorm(x,mean=Mu,sigma=fit.path$w[,,i]*vinv.tcross,log=TRUE)))
##                    }
##                    )
##   degfree <- lapply(seq(length(fit.path$rholist)),
##                     function(i){
##                       wi <- fit.path$wi[,,i]
##                       wi[abs(wi)<10^{-3}]<-0
##                       wi <- wi*v.tcross
##                       p <- ncol(wi)
##                       n.zero<-sum(wi==0)
##                       return(p+(p*(p+1)/2)-n.zero/2)
##                     }
##                     )
##   loglik <- simplify2array(loglik,higher=TRUE)
##   degfree <- simplify2array(degfree,higher=TRUE)
##   myscore <- -loglik+log(nrow(x))*degfree/2

##   if(is.null(lambda)){
##     lambda <- 0.5*nrow(x)*fit.path$rholist
##   }
  
##   if(ncol(x)<nrow(x)){
##     loglik.la0 <- sum(dmvnorm(x,mean=Mu,sigma=var(x),log=TRUE))
##     degfree.la0 <- ncol(x)+(ncol(x)*(ncol(x)+1)/2)
##     myscore.la0 <- -loglik.la0+log(nrow(x))*degfree.la0/2
##     myscore <- c(myscore.la0,myscore)
##     index.opt <- which.min(myscore)
##     if(index.opt==1){
##       w <- var(x)
##       wi <- solve(var(x))
##     }else{
##       wi <- fit.path$wi[,,index.opt-1]
##       wi[abs(wi)<10^{-3}]<-0
##       wi <- wi*v.tcross
##       w <- fit.path$w[,,index.opt-1]*vinv.tcross
##     }
##     lambda <- c(0,lambda)
##   }else{
##     index.opt <- which.min(myscore)
##     wi <- fit.path$wi[,,index.opt]
##     wi[abs(wi)<10^{-3}]<-0
##     wi <- wi*v.tcross
##     w <- fit.path$w[,,index.opt]*vinv.tcross
##   }
##   if (plot.it){
##     plot(lambda,myscore,type='b',xlab='lambda')
##   }
  
##   list(lambda=lambda,bic.score=myscore,Mu=Mu,wi=wi,w=w)
## }

## screen_bic.glasso.invcor <- function(x,length.lambda=20,trunc.k=5,plot.it=TRUE,include.mean=TRUE){

##   gridmax <- nrow(x)*rho.max(cor(x))/2
##   gridmin <- gridmax/length.lambda
##   my.grid <- lambdagrid_mult(gridmin,gridmax,length.lambda)[length.lambda:1]
  
##   fit.bicgl <- bic.glasso.invcor(x,lambda=my.grid,penalize.diagonal=FALSE,plot.it=plot.it)
##   wi.trunc <- fit.bicgl$wi
##   diag(wi.trunc) <- 0
##   nonzero <- min(2*ceiling(ncol(x)*ceiling(nrow(x)/trunc.k)/2),sum(wi.trunc!=0))
##   wi.trunc[-order(abs(wi.trunc),decreasing=TRUE)[1:nonzero]] <- 0
##   diag(wi.trunc) <- diag(fit.bicgl$wi)

##   list(wi=wi.trunc)
## }

## glasso.invcor <- function(s,rho,penalize.diagonal,term=10^{-3}){
##   if(penalize.diagonal==FALSE){
##     ww <- diag(s)
##     gl <- glasso(s,rho=rho*ww,penalize.diagonal=penalize.diagonal)
##     return(list(w=gl$w,wi=gl$wi))
##   }
##   if(penalize.diagonal==TRUE){
##     ww <- diag(s)/(1-rho)
##     gl <- glasso(s,rho=rho*ww,penalize.diagonal=penalize.diagonal)
##     return(list(w=gl$w,wi=gl$wi))
##   }
## }
  
## rhomax <- function(x){
##   n <- nrow(x)
##   s.var <- var(x)
##   diag(s.var) <- 0
##   return(max(abs(s.var)))
## }

## ##Lambda-grid
## rhogrid_mult <- function(rho.min,rho.max,nr.gridpoints){
##     mult.const <- (rho.max/rho.min)^(1/(nr.gridpoints-1))
##     return(rho.min*mult.const^((nr.gridpoints-1):0))
## }

## ##Crossvalidation Plot
## error.bars <- function (x, upper, lower, width = 0.02, ...)
## {
##     xlim <- range(x)
##     barw <- diff(xlim) * width
##     segments(x, upper, x, lower, ...)
##     segments(x - barw, upper, x + barw, upper, ...)
##     segments(x - barw, lower, x + barw, lower, ...)
##     range(upper, lower)
## }
## plotCV <- function(lambda,cv,cv.error,se=TRUE,type='b',...){
##   if (se==TRUE){ylim <- range(cv,cv+cv.error,cv-cv.error)}
##   else {ylim <- range(cv)}
##   plot(lambda,cv,type=type,ylim=ylim,...)
##   if (se)
##     error.bars(lambda,cv+cv.error,cv-cv.error,width=1/length(lambda))
##   invisible()
## }
## cv.fold <- function(n,folds=10){
##   split(sample(1:n),rep(1:folds,length=n))
## }

## ##' Crossvalidation for GLasso 
## ##'
## ##' 8! lambda-grid has to be increasing (see glassopath)
## ##' @title Crossvalidation for GLasso
## ##' @param x 
## ##' @param folds 
## ##' @param lambda lambda-grid (increasing!)
## ##' @param penalize.diagonal 
## ##' @param plot.it 
## ##' @param se 
## ##' @param include.mean 
## ##' @return 
## ##' @author n.stadler
## cv.glasso <- function(x,folds=10,rho,penalize.diagonal=FALSE,plot.it=FALSE,se=TRUE,include.mean=TRUE)
## {
##   colnames(x)<-paste('x',1:ncol(x),sep='')  
##   all.folds <- cv.fold(nrow(x),folds)
##   residmat <- matrix(NA,folds,length(rho))
##   for (cvfold in 1:folds){
##     omit <- all.folds[[cvfold]]
##     s <- var(x[-omit,])
##     if (include.mean==TRUE){
##       mu <- colMeans(x[-omit,,drop=FALSE])
##     }else{
##       mu <- rep(0,ncol(x))
##     }
##     fit.path <- glassopath(s,rholist=rho,penalize.diagonal=penalize.diagonal,trace=0)
##     if(length(omit)==1){
##       residmat[cvfold,] <- -2*apply(fit.path$w,3,dmvnorm,log=TRUE,mean=mu,x=x[omit,,drop=FALSE])
##     }else{
##       residmat[cvfold,] <- apply(-2*apply(fit.path$w,3,dmvnorm,log=TRUE,mean=mu,x=x[omit,,drop=FALSE]),2,sum)
##     }
##   }
##   cv <- apply(residmat,2,mean)
##   cv.error <- sqrt(apply(residmat,2,var)/folds)
##   gl.opt<-glasso(var(x),rho=rho[which.min(cv)],penalize.diagonal=penalize.diagonal)
##   cat('rho.opt:',rho[which.min(cv)],'\n')
##   w<-gl.opt$w
##   wi<-gl.opt$wi
##   wi[abs(wi)<10^{-3}]<-0
##   colnames(w)<-rownames(w)<-colnames(wi)<-rownames(wi)<-colnames(x)  

##   object <- list(rho.opt=rho[which.min(cv)],rho=rho,residmat=residmat,cv=cv,cv.error=cv.error,w=w,wi=wi,mu=colMeans(x))
##   if (plot.it){
##     plotCV(rho,cv,cv.error,se=se)
##   }
##   invisible(object)
## }
  
## gsea.rhoopt <- function(x1,x2,gene.sets,gene.name,length.la=20){

##   rho.optimal <- sapply(seq(length(gene.sets)),
##                         function(i){
##                           y <- gene.sets[[i]]
##                           cat('gene set:',i,'\n')
##                           ind.genes <- gene.name%in%y
##                           y1 <- x1[,ind.genes]
##                           y2 <- x2[,ind.genes]
                   
##                           gridmax <- rhomax(y1)
##                           gridmin <- 0.02*gridmax
##                           length.grid <- 20
##                           my.grid <- rhogrid_mult(gridmin,gridmax,length.grid)[length.grid:1]
##                           fit.cv1<- cv.glasso(y1,rho=my.grid,plot.it=TRUE)
                          
##                           gridmax <- rhomax(y2)
##                           gridmin <- 0.02*gridmax
##                           length.grid <- 20
##                           my.grid <- rhogrid_mult(gridmin,gridmax,length.grid)[length.grid:1]
##                           fit.cv2<- cv.glasso(y2,rho=my.grid,plot.it=TRUE)

##                           gridmax <- rhomax(rbind(y1,y2))
##                           gridmin <- 0.02*gridmax
##                           length.grid <- 20
##                           my.grid <- rhogrid_mult(gridmin,gridmax,length.grid)[length.grid:1]
##                           fit.cv<- cv.glasso(rbind(y1,y2),rho=my.grid,plot.it=TRUE)

##                           return(c(fit.cv1$rho.opt,fit.cv2$rho.opt,fit.cv$rho.opt))
##                         })
  
##   return(list(rho.optimal=rho.optimal))
## }  



