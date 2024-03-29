% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov_optim_methods.R
\name{obj_fun}
\alias{obj_fun}
\title{objective function for BFGS implementation using covariance matrices}
\usage{
obj_fun(A, B, tau)
}
\arguments{
\item{A}{Target Covariance Matrix}

\item{B}{list of background covariance(s).}

\item{tau}{lowerbound for root finding}
}
\value{
value of objective function
}
\description{
objective function for BFGS implementation using covariance matrices
}
