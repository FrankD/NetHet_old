library(GSA)

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

my.ttest <- function(x1,x2){
  s <- sqrt((var(x1)*(length(x1)-1)+var(x2)*(length(x2)-1))/(length(x1)+length(x2)-2))
  tstat <- (mean(x1)-mean(x2))/(s*sqrt(1/length(x1)+1/length(x2)))
  return(tstat)
}

agg.score.iriz <- function(ttstat,geneset,gene.name){
  ttstat <- ttstat[gene.name%in%geneset]
  aggscore <- sqrt(length(geneset))*mean(ttstat)
  return(aggscore)
}

gsea.iriz <- function(x1,x2,gene.sets,gene.name,genesets.name=NULL){
  
  no.genes <- ncol(x1)

  my.tstat <- sapply(1:no.genes,function(i){my.ttest(x1[,i],x2[,i])})
  agg.scores <- sapply(gene.sets,function(x){agg.score.iriz(my.tstat,x,gene.name)})
  pvals <- p.adjust(2*pnorm(-abs(agg.scores)),method='fdr')
  if(is.null(genesets.name)){
    return(list(agg.scores.sort=agg.scores[order(pvals)],pvals.sort=pval[order(pvals)],
                agg.scores=agg.scores,pvals=pvals))
  }else{
    return(list(agg.scores.sort=agg.scores[order(pvals)],pvals.sort=pvals[order(pvals)],
                genesets.name.sort=genesets.name[order(pvals)],
                agg.scores=agg.scores,pvals=pvals))
  }
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

gsea.highdimT2 <- function(x1,x2,gene.sets,gene.name,genesets.name=NULL,method='test.sd'){
  
  
  pvals<- sapply(gene.sets,
                 function(y){
                   ind.genes <- gene.name%in%y
                   eval(as.name(method))(x1[,ind.genes],x2[,ind.genes])$pval
                 })
  pvals.corrected <- p.adjust(pvals,method='fdr')
  if(is.null(genesets.name)){
    return(list(pvals.sort=pvals.corrected[order(pvals.corrected)],pvals=pvals.corrected))
  }else{
    return(list(pvals.sort=pvals.corrected[order(pvals.corrected)],
                genesets.name.sort=genesets.name[order(pvals.corrected)],
                pvals=pvals.corrected))
  }
  
}

gsea.diffnet.singlesplit <- function(x1,x2,gene.sets,gene.name,genesets.name=NULL,scale.mean=FALSE,scale.var=FALSE,...){

  pvals<- sapply(seq(length(gene.sets)),
                 function(i){
                   ## k <- length(y);print(k)
                   ## n1 <- nrow(x1)
                   ## n2 <- nrow(x2)
                   ## la <- sqrt(2*log(k)*sum(c(n1,n2))/2)/2
                   ## la <- seq(0.3,2,length=10)*la;print(la)
                   y <- gene.sets[[i]]
                   cat('gene set:',i,'\n')
                   ind.genes <- gene.name%in%y
                   diffnet_multisplit(scale(x1[,ind.genes],center=scale.mean,scale=scale.var),scale(x2[,ind.genes],center=scale.mean,scale=scale.var),b.splits=1,...)$pval.onesided
                 })
  pvals.corrected <- p.adjust(pvals,method='fdr')
  return(list(pvals=pvals.corrected))
}

aggpval <- function(pval,gamma.min=0.05){
  
  min(1,(1-log(gamma.min))*optimize(
                                    f=function(gamma){
                                      min(quantile(pval[!is.na(pval)]/gamma,probs=gamma),1)
                                    }
                                    ,interval=c(gamma.min,1),maximum=FALSE)$objective)
}

gsea.diffnet.multisplit <- function(x1,x2,no.splits=50,gene.sets,gene.name,genesets.name=NULL,scale.mean=FALSE,scale.var=FALSE,...){

  res <- lapply(seq(no.splits),
                function(i){
                  cat('split:',i,'\n')
                  pvals <- gsea.diffnet.singlesplit(x1,x2,gene.sets,gene.name,genesets.name=NULL,scale.mean=FALSE,scale.var=FALSE,...)$pvals
                  return(pvals)
                }
                )
  res <- matrix(simplify2array(res,higher=TRUE),nrow=length(gene.sets),ncol=no.splits)
  pval.agg <- apply(res,1,aggpval)

  if(is.null(genesets.name)){
    return(list(pvalagg.sort=pval.agg[order(pval.agg)],pvalagg=pval.agg,pval=res))
  }else{
    return(list(pvalagg.sort=pval.agg[order(pval.agg)],
                genesets.name.sort=genesets.name[order(pval.agg)],
                pvalagg=pval.agg,pval=res))
  }
}

par.gsea.diffnet.multisplit <- function(x1,x2,no.splits=50,gene.sets,gene.name,genesets.name=NULL,scale.mean=FALSE,scale.var=FALSE,...){

  res <- mclapply(seq(no.splits),
                function(i){
                  cat('split:',i,'\n')
                  pvals <- gsea.diffnet.singlesplit(x1,x2,gene.sets,gene.name,genesets.name=NULL,scale.mean=FALSE,scale.var=FALSE,...)$pvals
                  return(pvals)
                },
                  mc.set.seed=TRUE, mc.preschedule = TRUE)
  res <- matrix(simplify2array(res,higher=TRUE),nrow=length(gene.sets),ncol=no.splits)
  pval.agg <- apply(res,1,aggpval)

  if(is.null(genesets.name)){
    return(list(pvalagg.sort=pval.agg[order(pval.agg)],pvalagg=pval.agg,pval=res))
  }else{
    return(list(pvalagg.sort=pval.agg[order(pval.agg)],
                genesets.name.sort=genesets.name[order(pval.agg)],
                pvalagg=pval.agg,pval=res))
  }
}

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



