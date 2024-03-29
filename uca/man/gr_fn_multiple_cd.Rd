% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_optim_methods.R
\name{gr_fn_multiple_cd}
\alias{gr_fn_multiple_cd}
\title{gradient for coordinate descent BFGS implementation in multi-background
 setting without covariance matrix}
\usage{
gr_fn_multiple_cd(lambda, B, t_A, t_B, lambda_B, right, svd_right, j)
}
\arguments{
\item{lambda}{the lagrange multiplier}

\item{B}{background interested in calculating lagrange multiplier for}

\item{t_A}{transpose of A, a p x n1 data matrx}

\item{t_B}{list of transpose of B, a p x n2 data matrix}

\item{lambda_B}{list of background lagrange multiplier}

\item{right}{rbind(A, B), a (n1 + n2) x p data matrix}

\item{svd_right}{svd of right data matrix object}

\item{j}{which background matrix are we focusing on}
}
\value{
value of the gradient fn
}
\description{
gradient for multiple background implementation of data method
using coordinate descent
}
