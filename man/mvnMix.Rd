% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MVNSeq.r
\name{mvnMix}
\alias{mvnMix}
\alias{gmvnMix}
\alias{grmvnMix}
\title{Multivariate Normal Mixture Models}
\usage{
mvnMix(
  y,
  cl,
  min.iters = 10,
  max.iters = 50,
  tol = 0.001,
  verbose = interactive()
)

gmvnMix(
  y,
  cl,
  gr,
  common.fractions = FALSE,
  min.iters = 10,
  max.iters = 50,
  tol = 0.001,
  verbose = interactive()
)

grmvnMix(
  y,
  cl,
  gr,
  common.fractions = FALSE,
  min.iters = 10,
  max.iters = 100,
  tol = 0.001,
  verbose = interactive()
)
}
\arguments{
\item{y}{the sequence of observations}

\item{cl}{an integer vector allocating observations to classes}

\item{min.iters}{minimum number of EM iterations}

\item{max.iters}{maximum number of EM iterations}

\item{tol}{tolerance for the log likelihood}

\item{verbose}{should the log likelihood be reported.}

\item{gr}{an integer vector allocating observations to groups}

\item{common.fractions}{should the mixing fractions be common
across groups}
}
\value{
the fitted model
}
\description{
Fit multivariate Normal mixture models by EM.
}
\details{
These functions fit K component multivariate Normal mixtures to
sequences of observations. \code{mvnMix} fits a mixture model to a
single sequence of observations. \code{gmvnMix} and
\code{grmvnMix} fits separate mixture models to several groups
(sequences) of observations. \code{gmvnMix} fits a mixture to each
group so that each component has a different mean across groups,
but a common covariance. \code{grmvnMix} constrains the means of
the components to be Normally distributed across groups.
}
