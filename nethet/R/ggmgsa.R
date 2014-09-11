######################################################################
# Multivariate gene-set testing based on graphical models (GGM-GSA)
# -------------------------------------------------------------------
#
#
######################################################################


#####################
##Required packages##
#####################

#library(GSA)
#library('limma')
#library('multtest')
#library(ICSNP)
#library(multicore)
#library(DiffNet)
#library(mvtnorm)

#############################################################
##Some functions for data-format conversion (GSEA specific)##
#############################################################

##' Reads a gene expression dataset in GCT format and converts it into an R data frame
##'
##' 
##' @title Reads a gene expression dataset in GCT format and converts it into an R data frame
##' @param filename no descr
##' @return no descr
##' @author n.stadler
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

##' Reads a gene expression dataset in GCT format and converts it into an R data frame
##'
##' 
##' @title Reads a gene expression dataset in GCT format and converts it into an R data frame
##' @param filename no descr
##' @return no descr
##' @author n.stadler
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

##' Reads a gene expression dataset in RES format and converts it into an R data frame
##'
##' 
##' @title Reads a gene expression dataset in RES format and converts it into an R data frame
##' @param filename no descr
##' @return R data frame
##' @author n.stadler
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

##' Reads a class vector CLS file and defines phenotype
##' and class labels vectors for the samples in a gene expression file (RES or GCT format)
##'
##' 
##' @title Reads a class vector CLS file
##' @param file no descr
##' @return no descr
##' @author n.stadler
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

##' Pvalue adjustment
##'
##' 
##' @title Pvalue adjustment
##' @param p p-values
##' @param method method for p-value adjustment (default='fdr')
##' @return adjusted p-values
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
##' Tests for shift and change in scale of distribution.
##' 
##' @title Irizarry approach (shift and scale)
##' @param x1 expression matrix (condition 1)
##' @param x2 expression matrix (condition 2)
##' @param gene.sets list of gene-sets
##' @param gene.names gene names
##' @param gs.names gene-set names
##' @param method.p.adjust method for p-value adjustment (default: 'fdr')
##' @param alternative default='two-sided': two-sided p-values
##' @return list consisting of
##' \item{pval.shift}{p-values measuring shift}
##' \item{pval.scale}{p-values measuring scale}
##' \item{pval.combined}{combined p-values (minimum of pval.shift and pval.scale)}
##' @author n.stadler
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

##############
##GGM-GSA   ##
##############

##' Filter "non-normal" genes
##'
##' dicards genes with Shapiro-Wilk p-value (corrected for multiple comparision)
##' smaller than sign.level in either of the two conditions
##' (we used sign.level=0.001 in the GGMGSA paper)
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

##' P-value aggregation (inf-quantile formula of Meinshausen et al 2009)
##' 
##' 
##' @title P-value aggregation (inf-quantile formula of Meinshausen et al 2009)
##' @param pval p-values
##' @param gamma.min see inf-quantile formula of Meinshausen et al 2009 (default=0.05)
##' @return aggregated p-value
##' @author n.stadler
##' @export
aggpval <- function(pval,gamma.min=0.05){
  
  min(1,(1-log(gamma.min))*optimize(
                                    f=function(gamma){
                                      min(quantile(pval[!is.na(pval)]/gamma,probs=gamma),1)
                                    }
                                    ,interval=c(gamma.min,1),maximum=FALSE)$objective)
}

##' Single-split GGM-GSA
##'
##' 
##' @title ggmgsa_singlesplit
##' @param x1 centered (scaled) data for condition 1
##' @param x2 centered (scaled) data for condition 2
##' @param gene.sets list of gene-sets
##' @param gene.names gene names
##' @param method.p.adjust method for p-value adjustment (default: 'fdr')
##' @param ... other arguments (see diffnet_singlesplit)
##' @return list of results
##' @author n.stadler
ggmgsa_singlesplit <- function(x1,x2,gene.sets,gene.names,method.p.adjust='fdr',...){
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

##' Multi-split GGM-GSA
##'
##' 
##' @title ggmgsa_multisplit
##' @param x1 expression matrix for condition 1 (mean centered !)
##' @param x2 expression matrix for condition 2 (mean centered !)
##' @param no.splits number of random data splits (default: 50)
##' @param gene.sets list of gene-sets
##' @param gene.names gene names
##' @param gs.names gene-set names
##' @param method.p.adjust method for p-value adjustment (default: 'fdr')
##' @param order.adj.agg order of "pvalue-aggregation / -adjustment" (default='agg-adj')
##' @param ... other arguments (see diffnet_singlesplit)
##' @return list consisting of
##' \item{pvalmed}{median aggregated p-values}
##' \item{pvalagg}{inf-quantile aggregated p-values}
##' \item{pval}{p-value matrix (number of gene-sets)x(number of splits); before correction and adjustement}
##' @author n.stadler
##' @export
##' @example ../ggmgsa-pkg_test.r
ggmgsa_multisplit <- function(x1,x2,no.splits=50,gene.sets,gene.names,gs.names=NULL,method.p.adjust='fdr',order.adj.agg='agg-adj',...){

  res <- lapply(seq(no.splits),
                function(i){
                  cat('split:',i,'\n')
                  res <- ggmgsa_singlesplit(x1,x2,gene.sets=gene.sets,gene.names=gene.names,
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

##' Multi-split GGM-GSA
##'
##' computation parallelized over many data splits
##' 
##' @title par_ggmgsa_multisplit (parallelized computation)
##' @param x1 expression matrix for condition 1 (mean centered !)
##' @param x2 expression matrix for condition 2 (mean centered !)
##' @param no.splits number of random data splits (default: 50)
##' @param gene.sets list of gene-sets
##' @param gene.names gene names
##' @param gs.names gene-set names
##' @param method.p.adjust method for p-value adjustment (default: 'fdr')
##' @param order.adj.agg order of "pvalue-aggregation / -adjustment" (default='agg-adj')
##' @param ... other arguments (see diffnet_singlesplit)
##' @return list consisting of
##' \item{pvalmed}{median aggregated p-values}
##' \item{pvalagg}{inf-quantile aggregated p-values}
##' \item{pval}{p-value matrix (number of gene-sets)x(number of splits); before correction and adjustement}
##' @author n.stadler
##' @export
##' @example ../ggmgsa-pkg_test.r
par_ggmgsa_multisplit <- function(x1,x2,no.splits=50,gene.sets,gene.names,gs.names=NULL,method.p.adjust='fdr',order.adj.agg='agg-adj',...){

  res <- mclapply(seq(no.splits),
                  function(i){
                    cat('split:',i,'\n')
                    res <- ggmgsa_singlesplit(x1,x2,gene.sets=gene.sets,gene.names=gene.names,
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

