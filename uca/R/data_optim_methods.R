#' BFGS implementation for single background using the data method
#' (without using covariance matrix)
#'
#' Compute the UCA dimensions for single background using SVD and QR
#' using optim. Faster for high-dimensional data because we avoid
#' constructing the covariance matrix
#' @param A a n1 x p data matrix
#' @param B a n2 x p data matrix
#' @param maxit maximum iterations
#' @importFrom stats optim
#' @return list of two elements:
#' \itemize{
#'  \item values: eigenvalues
#'  \item tau: the contrastive parameter
#'  }
optim_data_cd <- function(A, B, maxit = 5E2L) {
  right <-  rbind(A, B)
  svd_right <- arma_svd(right)
  t_A <- t(A)
  t_B <- t(B)
  optim_with_grad <- optim(par = 2,
                          fn = obj_fun_cpp,
                          gr = gr_fun_cpp,
                          t_A = t_A,
                          t_B = t_B,
                          B = B,
                          right = right,
                          right_u = svd_right$u,
                          right_d = svd_right$d,
                          method = "L-BFGS-B",
                          lower = 0,
                          control = list(maxit = maxit))
  if (optim_with_grad$convergence != 0) {
      warning("optim did not converge\n")
  }

  return(list(values = optim_with_grad$value - optim_with_grad$par,
              tau = optim_with_grad$par))
}

#' objective function for coordinate descent BFGS implementation in
#' multi-background setting without covariance matrix
#'
#' objective function for multiple background implementation of data method
#' using coordinate descent
#' @param B background interested in calculating lagrange multiplier for
#' @param t_A transpose of A, a p x n1 data matrx
#' @param t_B list of transpose of B, a p x n2 data matrix
#' @param right rbind(A, B), a (n1 + n2) x p data matrix
#' @param svd_right svd of right data matrix object
#' @param lambda the lagrange multiplier
#' @param lambda_B list of background lagrange multiplier
#' @param j which background matrix are we focusing on
#' @return value of the objective function
#'
obj_fn_multiple_cd <- function(lambda, B, t_A, t_B, lambda_B, right,
                               svd_right, j) {
  lambda_B[[j]] <-  -lambda * t_B[[j]]
  left <- do.call(cbind, c(list(t_A), lambda_B))
  objective_fn <- multiple_score_calc_cpp(left,
                                          right,
                                          svd_right$u,
                                          svd_right$d,
                                          B,
                                          lambda)
  return(objective_fn$values + lambda)
}

#' gradient for coordinate descent BFGS implementation in multi-background
#'  setting without covariance matrix
#'
#' gradient for multiple background implementation of data method
#' using coordinate descent
#' @param B background interested in calculating lagrange multiplier for
#' @param t_A transpose of A, a p x n1 data matrx
#' @param t_B list of transpose of B, a p x n2 data matrix
#' @param right rbind(A, B), a (n1 + n2) x p data matrix
#' @param svd_right svd of right data matrix object
#' @param lambda the lagrange multiplier
#' @param lambda_B list of background lagrange multiplier
#' @param j which background matrix are we focusing on
#' @return value of the gradient fn

gr_fn_multiple_cd <- function(lambda, B, t_A, t_B, lambda_B, right,
                              svd_right, j) {
  lambda_B[[j]] <-  -lambda * t_B[[j]]
  left <- do.call(cbind, c(list(t_A), lambda_B))
  objective_gr <- multiple_score_calc_cpp(left,
                                          right,
                                          svd_right$u,
                                          svd_right$d,
                                          B,
                                          lambda)
  return(objective_gr$score)
}

#' coordinate descent BFGS implementation in multi-background setting
#' without covariance matrix
#'
#' Compute the UCA lagrange multiplier for a single background in the multiple
#' backgrounds using SVD and QR in a coordinate descent manner using optim.
#' Faster for high-dimensional data since we avoid creating a covariance matrix
#' @param B_focus background interested in calculating lagrange multiplier for
#' @param t_A transpose of A, a p x n1 data matrx
#' @param t_B list of transpose of B, a p x n2 data matrix
#' @param right rbind(A, B), a (n1 + n2) x p data matrix
#' @param svd_right svd of right data matrix object
#' @param lambda the lagrange multiplier
#' @param j index of background matrix focused on
#' @param maxit maximum number of iterations for optim run. Default 5E2
#' @importFrom stats optim
#' @return list of two elements:
#' \itemize{
#'  \item values: optimal eigenvalue
#'  \item tau: the contrastive parameter (lagrange multiplier) for single background
#'  in the multi. background setting
#'  }
optim_cd_multiple <- function(B_focus, t_A, t_B, right, svd_right, lambda,
                              j, maxit = 5E2) {

  #constants that don't really change if focused on j-th background
  old_right <- t(do.call(cbind, c(list(t_A),  t_B[-j])))
  lambda_B <- Map("*", -lambda, t_B)
  old_left <- do.call(cbind, c(list(t_A), lambda_B[-j]))
  svd_right_check <- arma_svd(old_right)

  #checking bounds
  f_val <- multiple_score_calc_cpp(left = old_left,
                                    right = old_right,
                                    right_u = svd_right_check$u,
                                    right_d = svd_right_check$d,
                                    tau = 0,
                                    B = B_focus)
  if (f_val$score >= 0) {
    warning("Redundant Constraint: Lagrange Multiplier is negative. Setting lambda to 0 \n")
    return(f_val);
  }else{
  optim_res <- optim(par = 3,
        fn = obj_fn_multiple_cd,
        gr = gr_fn_multiple_cd,
        t_A = t_A,
        t_B = t_B,
        lambda_B = lambda_B,
        right = right,
        svd_right = svd_right,
        B = B_focus,
        j = j,
        method = "L-BFGS-B",
        lower = 0,
        control = list(maxit = maxit))
  }
  if (optim_res$convergence != 0) {
      warning("optim did not converge\n")
  }

return(list(values = optim_res$value - optim_res$par,
            tau = optim_res$par))
}

