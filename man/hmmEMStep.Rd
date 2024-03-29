% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MVNSeq.r
\name{hmmEMStep}
\alias{hmmEMStep}
\title{HMM EM Update}
\usage{
hmmEMStep(K, y, pars)
}
\arguments{
\item{K}{the number of classes}

\item{y}{the sequence of observations}

\item{pars}{current sequence parameters}
}
\value{
the updated set of parameters
}
\description{
The HMM update step of the EM algorithm.
}
\details{
Perform a single step of the EM algorithm to fit the HMM component
of a multivariate Normal HMM model.
}
