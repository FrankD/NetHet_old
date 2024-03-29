% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/mixglasso.R
\name{sparse_conc}
\alias{sparse_conc}
\title{Generates sparse inverse covariance matrices}
\usage{
sparse_conc(p, K, s, s.common, magn.nz = 0.5, scale.parcor = TRUE)
}
\arguments{
\item{p}{Dimensionality of inverse covariance matrix}

\item{K}{Number of inverse covariance matrices}

\item{s}{Number of non-zero entries per inverse covariance matrix}

\item{s.common}{Number of non-zero entries shared across different inverse covariance matrices}

\item{magn.nz}{Magnitude of non-zero elements}

\item{scale.parcor}{Should SigInv be scaled to have diagonal equal one, siginv=parcor ?}
}
\value{
SigInv: list of inverse covariance matrices
}
\description{
Generates sparse inverse covariance matrices
}
\author{
n.stadler
}
\keyword{internal}

