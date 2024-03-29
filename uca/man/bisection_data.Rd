% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_bisection_methods.R
\name{bisection_data}
\alias{bisection_data}
\title{Use bisection to find the optimal lagrange multiplier for a single background
without constructing covariance matrices}
\usage{
bisection_data(A, B, limit = 20L, maxit = 100000L, tol = 1e-06)
}
\arguments{
\item{A}{a n1 x p Target data matrix}

\item{B}{a n2 x p Background data matrix}

\item{limit}{upperbound for the lagrange multiplier.}

\item{maxit}{maximum iterations}

\item{tol}{tolerance for convergence criteria}
}
\value{
list of two elements:
\itemize{
 \item tau: the lagrange multiplier
 \item score: eigenvalue associated with tau
 }
}
\description{
Use bisection method to find the optimal Lagrange multiplier using data
matrices and a single background matrix. This is useful when constructing
the covariance matrix is computationally intensive.
}
