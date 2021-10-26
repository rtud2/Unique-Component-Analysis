#' @importFrom Rcpp evalCpp
#' @useDynLib uca, .registration = TRUE
NULL

#' center_f
#'
#' Convenient function to center the data, rather than typing `scale' or `sweep`
#'
#' @param X a matrix
#' @return centered data matrix X
#' @export
#'
center_f <- function(X) {
  column_means <- colMeans(X)
  return(sweep(X, 2, column_means, "-"))
}

#' SVD of a difference of covariance matrices without constructing cov matrix
#'
#' calculate SVD of a product of matrices by using svd and QR decompositions
#'
#' @param left left side of a product
#' @param right right side of a product
#' @param nv number of unique components
#' @return top nv eigenvalues and associated eigenvectors
#' 
product_svd_R <- function(left, right, nv) {
  svd_right <- svd(right,nv = 0)
  qr_left_U <- qr(left %*% svd_right$u)

  RS_svd <- arma_svd(t(t(qr.R(qr_left_U)) * svd_right$d))
  u <- qr.Q(qr_left_U) %*% RS_svd$u

  #calculates diag(crossprod(u, left) %*% (right %*%u))
  eigenvalues <- colSums(u * (left %*% (right %*% u)))
  top_eig_vals <- order(eigenvalues, decreasing = T)[1:nv]
  list(values = eigenvalues[top_eig_vals],
       vectors = u[, top_eig_vals])
}


