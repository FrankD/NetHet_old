##source code
source('../code/hmmgl_mixgl-22102012-backup-mixglasso-only.r')

##load packages
library(cluster)

##some functions
jaccard_coef <- function(set1,set2){
  return(length(intersect(set1,set2))/length(union(set1,set2)))
}

mixglasso_ss <-  function(x,k=3){
  launi <- sqrt(2*nrow(x)*log(ncol(x)))/2
  fit <- mixglasso(x,n.comp=k,lambda=launi,pen='glasso.parcor',
                   init='kmeans',nstart.kmeans=100,iter.max.kmeans=1000,
                   term=10^{-3},min.compsize=5,show.trace=TRUE)
  list(cl=fit$comp)
}

mixunpen_ss <-  function(x,k=3){
  la <- 0
  p <- ncol(x)
  fit <- mixglasso(x,n.comp=k,lambda=la,pen='glasso.parcor',
                   init='kmeans',nstart.kmeans=100,iter.max.kmeans=1000,
                   term=10^{-3},min.compsize=p,show.trace=TRUE,miniter=0)
  list(cl=fit$comp)
}

kmeans_ss <- function(x,k){
  list(cl=kmeans(x, centers=k, iter.max = 1000, nstart = 100)$cluster)
}

my_clustboot <- function(x,B,frac.sub=0.5,cl.orig=NULL,bootmethod='subset',method=mixglasso_ss,k=3){
  ##fit clustering on all data
  if(is.null(cl.orig)){
    cl.orig <- method(x,k)$cl
  }

  ##fit on subsampled data 
  n <- nrow(x)
  fitss <- mclapply(1:B,function(i){
    n <- nrow(x)
    if(bootmethod=='subset'){
      n.sub <- ceiling(n*frac.sub)
      my.ss <- sample(1:n,n.sub,replace=FALSE)
    }
    if(bootmethod=='boot'){
      my.ss <- sample(1:n,n,replace=TRUE)
    }
    if(bootmethod=='boot_nodup'){
      my.ss <- sample(1:n,n,replace=TRUE)
      my.ss <- my.ss[!duplicated(my.ss)]
    }
    x.ss <- x[my.ss,]
    my.cl <- rep(NA,n)
    my.cl[my.ss[!duplicated(my.ss)]] <- method(x.ss,k)$cl[!duplicated(my.ss)]
    return(list(cl=my.cl))
  },mc.preschedule = FALSE, mc.set.seed = TRUE)
  clss <- sapply(fitss,function(x){x[['cl']]})

  ##compute stability measures
  jaccss <- matrix(NA,length(unique(cl.orig)),B)
  for(i in 1:length(unique(cl.orig))){
    for (b in 1:B){
      jaccss[i,b] <- max(sapply(1:length(unique(clss[,b])),function(j){
        cl.orig.ss <- cl.orig
        cl.orig.ss[is.na(clss[,b])] <- NA
        jaccard_coef(which(clss[,b]==j),which(cl.orig.ss==i))
      }))
    }
  }
  return(list(jaccss=jaccss,clss=clss,cl.orig=cl.orig))
}



