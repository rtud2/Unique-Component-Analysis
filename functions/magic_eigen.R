#' magic_eigen
#' 
#' calculate the eigendecomposition via a product of matrices using SVD and QR decomposition. Much faster for large matrices
#'
#' @param A *Centered* Target Data Matrix. n1 x p dimensions
#' @param B *Centered* Background Data Matrix. n2 x p dimensions
#' @param lambda contrastive parameter
#' @param k number of eigenvalues and vectors to compute
#' @return list of the largest k eigenvalues and  assiciated eigenvectors
#' @importFrom methods as
#' @importFrom RSpectra eigs_sym
#' @export

magic_eigen <- function(A, B, lambda, k){
  
  A_divided = A/sqrt(nrow(A) - 1)
  B_divided = B/sqrt(nrow(B) - 1)  
  
  left <- cbind(t(A_divided), - lambda * t(B_divided))
  right <- rbind(A_divided, B_divided)
  
  svd_right <- svd(right, nv = 0)
  qr_left_U <- qr(left %*% svd_right$u)
  RS_svd <- svd( sweep(qr.R(qr_left_U), 2, svd_right$d, FUN = "*") , nv = 0)
  
  u = qr.Q(qr_left_U) %*% RS_svd$u
  eigenvalues <- diag(crossprod(u, left) %*% (right %*% u))
  top_eig_vals <- order(eigenvalues, decreasing = T)[1:k]
  
  list(values = eigenvalues[top_eig_vals],
       vectors = u[ ,top_eig_vals])
}

# t_happy2 <- scale(t_happy, scale = F)
# t_neutral2 <- scale(t_neutral, scale = F)

###### for 100*136
# > microbenchmark::microbenchmark(magic_eigen(t_happy2, t_neutral2, 2, 5), eigs_sym(cov_happy - 2*cov_neutral, 5, "LA"), times = 10)
# Unit: seconds
# expr      min       lq     mean   median        uq       max neval cld
# magic_eigen(t_happy2, t_neutral2, 2, 5) 3.647757 3.667107 3.749784 3.696459  3.805642  4.017985    10  a 
# eigs_sym(cov_happy - 2 * cov_neutral, 5, "LA") 9.804725 9.869240 9.947251 9.920517 10.016461 10.162978    10   b


###### apparently cbinding transposes is easier than   rbinding then transpose.
# > microbenchmark::microbenchmark(t(rbind(A_divided, - lambda * B_divided)), cbind(t(A_divided), - lambda*t(B_divided)))
# Unit: milliseconds
# expr     min      lq     mean  median       uq      max neval cld
# t(rbind(A_divided, -lambda * B_divided)) 26.5603 26.8648 28.52919 27.0243 27.65665 142.5401   100   b
# cbind(t(A_divided), -lambda * t(B_divided)) 14.5453 16.4540 19.09077 16.5184 16.69970 140.5627   100  a 


#tmp_eigs <- magic_eigen(t_happy2, t_neutral2, 2, 5)
#tmp_eigs2 <- eigs_sym(cov_happy - 2*cov_neutral, 5, "LA")

#tmp_eigs$values - tmp_eigs2$values

#sign flips on eigenvectors
#sum(abs(tmp_eigs$vectors) - abs(tmp_eigs2$vectors))
