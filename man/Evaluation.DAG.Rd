% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Evaluation.DAG.R
\name{Evaluation.DAG}
\alias{Evaluation.DAG}
\title{Evaluation function for the estimated DAG.}
\usage{
Evaluation.DAG(estimated.adjace, true.adjace, type.adj=2)
}
\arguments{
\item{estimated.adjace}{The target data, a n * p matrix, where n is the sample size and p is data dimension.}

\item{true.adjace}{The auxiliary data in K auxiliary domains, a list with K elements, each of which is a nk * p matrix, where nk is the sample size of the k-th auxiliary domain.}

\item{type.adj}{The type of adjacency matrix. 1: the entries of matrix contains just two value, 0 and 1, which indicate the existence of edges; 2 (default): the matrix also measures connection strength, and 0 means no edge.}
}
\value{
A result list including Recall, FDR, F1score, MCC, Hamming Distance,and estimated error of adjacency matrix on F-norm.
}
\description{
Evaluation function for the estimated DAG.
}
\references{
Zhao, R., He X., and Wang J. (2022). Learning linear non-Gaussian directed acyclic graph with diverging number of nodes. Journal of Machine Learning Research.
}
\author{
Ruixaun Zhao \href{mailto:ruixuanzhao2-c@my.cityu.edu.hk}{ruixuanzhao2-c@my.cityu.edu.hk}.
}
