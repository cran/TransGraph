% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trans_mean.R
\name{trans_mean}
\alias{trans_mean}
\title{Transfer learning for mean estimation.}
\usage{
trans_mean(t.mean.m, A.mean, n, clambda=1)
}
\arguments{
\item{t.mean.m}{The estimated target p-dimensional mean vector, where p is mean dimension.}

\item{A.mean}{A K*p matrix with the k-th row being the estimated p-dimensional mean vector of the k-th auxiliary domain.}

\item{n}{The target sample size.}

\item{clambda}{The coefficients set in tuning parameters used in transfer learning for mean eatimation, and the default setting is clambda.m * sqrt( log(p) / n ).}
}
\value{
t.mean.m.hat: The transfer learning estimation of the target p-dimensional mean vector.
}
\description{
Transfer learning for mean estimation.
}
\references{
Ren, M. and Wang J. (2023). Local transfer learning of Gaussian graphical mixture models.
}
\author{
Mingyang Ren \href{mailto:renmingyang17@mails.ucas.ac.cn}{renmingyang17@mails.ucas.ac.cn}.
}
