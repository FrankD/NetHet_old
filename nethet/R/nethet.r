#' NetHet-package
#'
#' A bioconductor package for high-dimensional exploration of biological network heterogeneity 
#' 
#'
#' Includes:
#'   *Network-based clustering (MixGLasso)
#'   *Differential network (DiffNet)
#'   *Differential regression (DiffRegr)
#'   *Gene-set analysis based on graphical models (GGMGSA)
#'   *Plotting functions for exploring network heterogeneity
#'
#' @references St\"adler, N. and Mukherjee, S. (2013). Two-Sample Testing in High-Dimensional Models.
#' Preprint \url{http://arxiv.org/abs/1210.4584}.
#' @import glasso mvtnorm parcor GeneNet huge CompQuadForm ggm mclust multicore GSA limma multtest ICSNP glmnet 
#' @docType package
#' @name NetHet-package
#' @useDynLib nethet
NULL
