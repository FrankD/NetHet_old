% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/diffnet.R
\name{diffnet_multisplit}
\alias{diffnet_multisplit}
\title{Differential Network}
\usage{
diffnet_multisplit(x1, x2, b.splits = 50, frac.split = 1/2,
  screen.meth = "screen_bic.glasso", include.mean = FALSE,
  gamma.min = 0.05, compute.evals = "est2.my.ev3",
  algorithm.mleggm = "glasso_rho0", method.compquadform = "imhof",
  acc = 1e-04, epsabs = 1e-10, epsrel = 1e-10, show.warn = FALSE,
  save.mle = FALSE, verbose = TRUE, mc.flag = FALSE, mc.set.seed = TRUE,
  mc.preschedule = TRUE, mc.cores = getOption("mc.cores", 2L), ...)
}
\arguments{
\item{x1}{Data-matrix sample 1.
You might need to center and scale your data-matrix.}

\item{x2}{Data-matrix sample 1.
You might need to center and scale your data-matrix.}

\item{b.splits}{Number of splits (default=50).}

\item{frac.split}{Fraction train-data (screening) / test-data (cleaning) (default=0.5).}

\item{screen.meth}{Screening procedure. Options: 'screen_bic.glasso' (default),
'screen_cv.glasso', 'screen_shrink' (not recommended), 'screen_mb'.}

\item{include.mean}{Should sample specific means be included in hypothesis?
Use include.mean=FALSE (default and recommended) which assumes mu1=mu2=0
and tests the hypothesis H0: Omega_1=Omega_2.}

\item{gamma.min}{Tuning parameter in p-value aggregation of Meinshausen et al (2009). (Default=0.05).}

\item{compute.evals}{Method to estimate the weights in the weighted-sum-of-chi2s distribution.
The default and (currently) the only available option
is the method 'est2.my.ev3'.}

\item{algorithm.mleggm}{Algorithm to compute MLE of GGM. The algorithm 'glasso_rho' is the
default and (currently) the only available option.}

\item{method.compquadform}{Method to compute distribution function of weighted-sum-of-chi2s
(default='imhof').}

\item{acc}{See ?davies (default 1e-04).}

\item{epsabs}{See ?imhof (default 1e-10).}

\item{epsrel}{See ?imhof (default 1e-10).}

\item{show.warn}{Should warnings be showed (default=FALSE)?}

\item{save.mle}{If TRUE, MLEs (inverse covariance matrices for samples 1 and 2)
are saved for all b.splits. The median aggregated inverse covariance matrix
is provided in the output as 'medwi'. The default is save.mle=FALSE.}

\item{verbose}{If TRUE, show output progress.}

\item{mc.flag}{If \code{TRUE} use parallel execution for each b.splits via function
\code{mclapply} of package \code{parallel}.}

\item{mc.set.seed}{See mclapply. Default=TRUE}

\item{mc.preschedule}{See mclapply. Default=TRUE}

\item{mc.cores}{Number of cores to use in parallel execution. Defaults to
mc.cores option if set, or 2 otherwise.}

\item{...}{Additional arguments for screen.meth.}
}
\value{
list consisting of
\item{ms.pval}{p-values for all b.splits}
\item{ss.pval}{single-split p-value}
\item{medagg.pval}{median aggregated p-value}
\item{meinshagg.pval}{meinshausen aggregated p-value (meinshausen et al 2009)}
\item{teststat}{test statistics for b.splits}
\item{weights.nulldistr}{estimated weights}
\item{active.last}{active-sets obtained in last screening-step}
\item{medwi}{median of inverse covariance matrices over b.splits}
\item{sig.last}{constrained mle (covariance matrix) obtained in last cleaning-step}
\item{wi.last}{constrained mle (inverse covariance matrix) obtained in last cleaning-step}
}
\description{
Differential Network
}
\details{
Remark:

* If include.mean=FALSE, then x1 and x2 have mean zero and DiffNet tests
  the hypothesis H0: Omega_1=Omega_2. You might need to center x1 and x2.
* If include.mean=TRUE, then DiffNet tests the hypothesis
  H0: mu_1=mu_2 & Omega_1=Omega_2
* However, we recommend to set include.mean=FALSE and to test equality of the means
  separately.
* You might also want to scale x1 and x2, if you are only interested in
  differences due to (partial) correlations.
}
\examples{

############################################################
##This example illustrates the use of Differential Network##
############################################################


##set seed
set.seed(1)

##sample size and number of nodes
n <- 40
p <- 10

##specifiy sparse inverse covariance matrices
gen.net <- generate_2networks(p,graph='random',n.nz=rep(p,2),
                              n.nz.common=ceiling(p*0.8))
invcov1 <- gen.net[[1]]
invcov2 <- gen.net[[2]]
plot_2networks(invcov1,invcov2,label.pos=0,label.cex=0.7)

##get corresponding correlation matrices
cor1 <- cov2cor(solve(invcov1))
cor2 <- cov2cor(solve(invcov2))

##generate data under null hypothesis (both datasets have the same underlying
## network)
library('mvtnorm')
x1 <- rmvnorm(n,mean = rep(0,p), sigma = cor1)
x2 <- rmvnorm(n,mean = rep(0,p), sigma = cor1)

##run diffnet (under null hypothesis)
dn.null <- diffnet_multisplit(x1,x2,b.splits=1,verbose=FALSE)
dn.null$ss.pval#single-split p-value

##generate data under alternative hypothesis (datasets have different networks)
x1 <- rmvnorm(n,mean = rep(0,p), sigma = cor1)
x2 <- rmvnorm(n,mean = rep(0,p), sigma = cor2)

##run diffnet (under alternative hypothesis)
dn.altn <- diffnet_multisplit(x1,x2,b.splits=1,verbose=FALSE)
dn.altn$ss.pval#single-split p-value
dn.altn$medagg.pval#median aggregated p-value

##typically we would choose a larger number of splits
# dn.altn <- diffnet_multisplit(x1,x2,b.splits=10,verbose=FALSE)
# dn.altn$ms.pval#multi-split p-values
# dn.altn$medagg.pval#median aggregated p-value
# plot(dn.altn)#histogram of single-split p-values
}
\author{
n.stadler
}

