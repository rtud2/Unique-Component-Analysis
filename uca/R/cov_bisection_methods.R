#' score_calc
#'
#' Calculate the derivative of Lagrangian.
#' 
#' @param A Target Covariance Matrix
#' @param B list of background covariance(s).
#' @param tau lowerbound for root finding
#' @return list of tau (the lagrangian), eigenvector associated with tau, and derivative value of the lagrange multiplier
#' @importFrom methods as
#' @importFrom RSpectra eigs_sym
#' 
score_calc <- function(A, B, tau) {
  eigen_calc <- eigs_sym(A - tau * B, 1L, "LA")
  return(list(
    score = 1 - as(crossprod(eigen_calc$vectors,
                             B %*% eigen_calc$vectors),
                   "numeric"),
    values = eigen_calc$values,
    tau = tau))
  }


#' bisection_cov
#'
#' Use bisection method to find the optimal Lagrangian for solving unique component analysis (uca)
#' 
#' @param A Target Covariance Matrix
#' @param B Background Covariance Matrix. 
#' @param nv number of eigenvectors to use
#' @param limit upperbound of lagrange multiplier
#' @param maxit maxium number of iterations for the algorithm to run
#' @param tol tolerance for when to stop the algorithm
#' @return the final list of tau (the optimal lagrange multiplier), vector (eigenvector)
#' associated with tau, value (eigenvalue), and score (derivative value of the lagrangian)

bisection_cov <- function(A, B, limit = 50, maxit = 1E5L, nv = 1, tol = 1E-6) {
  f_val <- vector(mode = "list", length = 2L)
  f_val[[1]] <- score_calc(A, B, 0)
  og_upper_lim <- f_val[[2]]$tau <- limit

  if (f_val[[1]]$score >= 0) {
    warning("Redundant Constraint: Lagrange Multiplier is negative. Setting lambda to 0 \n");
    return(f_val[[1]]);
  }else{
    for (iter in 1L:maxit) {
      limit <- c(f_val[[1]]$tau, f_val[[2]]$tau)
      if (limit[2] - limit[1] < tol * limit[1]) {
          break
      }
      if (iter == maxit) {
          warning("maximum iteration reached: solution may not be optimal \n")
      }
      tau_score <- score_calc(A, B, sum(limit) / 2)

      if (tau_score$score < 0) {
        f_val[[1]] <- tau_score
      } else {
        f_val[[2]] <- tau_score
      }
    }
    if (round(tau_score$tau) == og_upper_lim) {
      warning(paste("Lagrange Multiplier is near upperbound.
                    Consider increasing the upperbound. current limit is",
                    og_upper_lim, "\n"))
    }
    return(f_val[[which.min(abs(c(f_val[[1]]$score, f_val[[2]]$score)))]])
  }
}

#' cov_multiple
#'
#' UCA for multiple background datasets, with covariance matrices.
#' Use different methods to find the optimal Lagrangian for solving UCA of each background dataset.
#'
#' @param A Target Covariance Matrix
#' @param B list of background covariance(s).
#' @param lambda initial guess on what the Lagrange Multiplier should be. Default NULL
#' @param nv number of uca components to estimate
#' @param max_iter maximum iterations for coordinate descent, if tolerance is not reached. default 1E5
#' @param tol tolerance for stopping criteria of coordinate descent. default 1E-6
#' @param algo which algorithm to use. default algo == "bisection". If algo = "cd", L-BFGS-S optimization is used instead for coordinate descent. algo== "gd" than L-BFGS-S optimizes all lambda's simultaenously (gradient descent). 
#' @param ... other parameters to pass in to bisection(...)
#' @return for each background covariance matrix in B, return list of tau (the optimal lagrange multiplier), vectors(eigenvector)
#'  associated with largest value (eigenvalue)
#' @importFrom RSpectra eigs_sym

cov_multiple <- function(A,
                         B,
                         lambda = NULL,
                         nv = 2L,
                         max_iter = 1E5L,
                         tol = 1E-6,
                         algo = "bisection",
                         ...) {
  #initialize starting point if one isn't supplied
  if (length(lambda) == 0) {
    lambda <- c(0, sapply(seq_along(B)[-1], function(zz) {
                              bisection_cov(A,B[[zz]])$tau
                    }))
    }
  score <- Inf;
  if (algo == "bisection") {
    for (i in 1L:max_iter) {
      old_score <- score
      for (j in seq_along(B)) {
        bisection_j <- bisection_cov(A - Reduce("+", Map("*", lambda[-j], B[-j])),
                                     B[[j]],
                                     ...)
        lambda[j] <- bisection_j$tau
      }
      score <- sum(bisection_j$values, lambda)
      if (abs(old_score - score) < tol * abs(old_score)) {
          break
      }
   }
  }else if (algo == "cd") {
    for (i in 1L:max_iter) {
      old_score <- score
      for (j in seq_along(B)) {
        bisection_j <- optim_cov_cd(A - Reduce("+", Map("*", lambda[-j], B[-j])), B[[j]], ...) 
        lambda[j] <- bisection_j$tau
      }
      score <- sum(bisection_j$values, lambda)
      if (abs(old_score - score) < tol * abs(old_score)) {
          break
      }
    }
  } else if (algo == "gd") {
        optim_gd_tmp <- optim_cov_bfgs_gd(A = A, B = B, maxit = max_iter)
        lambda <- optim_gd_tmp$tau
  } else {
    stop(paste("algo", algo, " not recognized"))
  }


  dca <- eigs_sym(A - Reduce("+", Map("*", lambda, B)), nv, "LA")
 return(list(values = dca$values, vectors = dca$vectors, tau = lambda))
}
