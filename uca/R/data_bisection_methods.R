#' @importFrom Rcpp evalCpp
#' @useDynLib uca, .registration = TRUE
NULL

#' Use bisection to find the optimal lagrange multiplier for a single background
#' without constructing covariance matrices
#'
#' Use bisection method to find the optimal Lagrange multiplier using data
#' matrices and a single background matrix. This is useful when constructing
#' the covariance matrix is computationally intensive.
#' @param A a n1 x p Target data matrix
#' @param B a n2 x p Background data matrix
#' @param limit upperbound for the lagrange multiplier.
#' @param maxit maximum iterations
#' @param tol tolerance for convergence criteria
#' @return list of two elements:
#' \itemize{
#'  \item tau: the lagrange multiplier
#'  \item score: eigenvalue associated with tau
#'  }
bisection_data <- function(A, B, limit = 20L, maxit = 1E5L, tol = 1E-6) {
  right <- rbind(A, B)
  svd_right <- arma_svd(right)
  t_A <- t(A)
  t_B <- t(B)

  svd_right_check <- arma_svd(A)

  f_val <- vector(mode = "list", length = 2L)
  f_val[[1]] <- multiple_score_calc_cpp(left = t_A,
                                    right = A,
                                    right_u = svd_right_check$u,
                                    right_d = svd_right_check$d,
                                    tau = 0,
                                    B = B)
  og_upper_lim <- f_val[[2]]$tau <- limit

  if (f_val[[1]]$score >= 0) {
    warning("Redundant Constraint: Lagrange Multiplier is negative.
            Setting lambda to 0 \n")
    return(f_val[[1]]);
  } else {
    for (iter in 1L:maxit) {
      limit <- c(f_val[[1]]$tau, f_val[[2]]$tau)
      if (limit[2] - limit[1] < tol * limit[1]) {
         break
      }
      if (iter == maxit) {
          warning("maximum iteration reached: solution may not be optimal \n")
      }
      lambda <- 0.5 * sum(limit)
      tau_score <- multiple_score_calc_cpp(left = cbind(t_A, -lambda * t_B),
                                      right = right,
                                      right_u = svd_right$u,
                                      right_d = svd_right$d,
                                      tau = lambda,
                                      B = B)

      if (tau_score$score < 0) {
        f_val[[1]] <- tau_score
      } else {
        f_val[[2]] <- tau_score
      }
    }

    if (round(tau_score$tau) == og_upper_lim) {
      warning("Lagrange Multiplier is near upperbound.
              Consider increasing the upperbound.(default is 20)\n")
    }
    return(f_val[[which.min(abs(c(f_val[[1]]$score, f_val[[2]]$score)))]])
  }
}

#' Use bisection to find the optimal lagrange multiplier for multiple background
#' without constructing covariance matrices
#'
#' Use bisection method to find the optimal Lagrange multiplier using data
#' matrices and multiple background matrices. This is useful when constructing
#' the covariance matrix is computationally intensive.
#' @param B_focus n x p data matrix of background matrix to solve lagrangian for
#' @param t_A precomputed A transpose
#' @param t_B precomputed B transpose
#' @param right precomputed right long matrix: rbind(A, B)
#' @param svd_right precomputed svd of right matrix
#' @param lambda j dimensional vector with candidate lagrange multipliers for
#' each background data matrix
#' @param j index of specific background you're solving for
#' @param limit upperbound for the lagrange multiplier
#' @param maxit maximum iterations
#' @param tol tolerance for convergence criteria
#' @return list of two elements:
#' \itemize{
#'  \item tau: the lagrange multiplier
#'  \item score: eigenvalue associated with tau
#'  }
bisection_data_multiple <- function(B_focus, t_A, t_B, right, svd_right, lambda,
                                    j, limit = 20L, maxit = 1E5, tol = 1E-6) {
  #constants that don't really change if focused on j-th background
  old_right <- t(do.call(cbind, c(list(t_A),  t_B[-j])))
  lambda_B <- Map("*", -lambda, t_B)
  old_left <- do.call(cbind, c(list(t_A), lambda_B[-j]))
  svd_right_check <- arma_svd(old_right)

  #checking bounds
  f_val <- vector(mode = "list", length = 2L)
  f_val[[1]] <- multiple_score_calc_cpp(left = old_left,
                                    right = old_right,
                                    right_u = svd_right_check$u,
                                    right_d = svd_right_check$d,
                                    tau = 0,
                                    B = B_focus)
  og_upper_lim <- f_val[[2]]$tau <- limit

  if (f_val[[1]]$score >= 0) {
    warning("Redundant Constraint: Lagrange Multiplier is negative.
            Setting lambda to 0 \n");
    return(f_val[[1]])
  } else {
    for (iter in 1:maxit) {
      limit <- c(f_val[[1]]$tau, f_val[[2]]$tau)
      if (limit[2] - limit[1] < tol * limit[1]) {
          break
      }
      if (iter == maxit) {
        warning("maximum iteration reached: solution may not be optimal \n")
      }

      lambda_B[[j]] <-  -0.5 * sum(limit) * t_B[[j]]

      tau_score <- multiple_score_calc_cpp(
                        left = do.call(cbind, c(list(t_A), lambda_B)),
                        right = right,
                        right_u = svd_right$u,
                        right_d = svd_right$d,
                        tau = 0.5 * sum(limit),
                        B = B_focus)

      if (tau_score$score < 0) {
       f_val[[1]] <- tau_score
      } else {
       f_val[[2]] <- tau_score
      }
  }
      if (round(tau_score$tau) == og_upper_lim) {
        warning("Lagrange Multiplier is near upperbound.
                Consider increasing the upperbound.(default is 20) \n")
      }
      return(f_val[[which.min(abs(c(f_val[[1]]$score, f_val[[2]]$score)))]])
    }
}

#' UCA for multiple backgrounds using data matrices
#'
#' Wrapper function for UCA with multiple background datasets with data
#' matrices. We use a sequence of SVD and QR decomposition (Product-SVD).
#' Use different algorithms to find the optimal Lagrangian
#' @param A  Target Data Matrix. n1 x p dimensions
#' @param B  List of k Background Data Matrix. n_k x p dimensions
#' @param lambda initial guess for Lagrange Multiplier. Default NULL.
#' algo = "bisection" only
#' @param nv number of uca components to estimate
#' @param algo which algorithm to use. default algo == "bisection". If
#' algo = "cd", L-BFGS-B optimization is used instead for coordinate descent.
#' @param max_iter maximum iterations for all algorithms if tolerance is not
#' reached. default 1E5.
#' @param tol tolerance for stopping criteria for all algorithms. default 1E-6
#' @param ... additional parameters to pass in to `bisection_data_multiple()`
#' and `optim_data_cd()`
#' @return list of three elements:
#' \itemize{
#'  \item tau: the optimal lagrange multiplier(s)
#'  \item values: eigenvalue associated with tau(s)
#'  \item vectors: top `nv` eigenvectors associated with tau(s)
#'  }
data_multiple <- function(A, B, lambda = NULL, nv = 2L, max_iter = 1E5L,
                          tol = 1E-6, algo = "bisection", ...) {
  #initialize starting point if one isn't supplied. greedy start
  if (length(lambda) == 0) {
    # use A and B here instead of A/B_divided b/c already divided in
    # `bisection_multiple`
    # do not pre-initialize first iteration since it's overwritten in step 1 of
    # coordinate descent.
    lambda <- sapply(seq_along(B),
                     function(zz) {
                         optim_data_cd(A, B[[zz]], ...)$tau
                     })
  }

  t_A <- t(A)
  t_B <- lapply(B, t)
  right <- t(do.call(cbind, c(list(t_A), t_B)))
  svd_right <- arma_svd(right)
  score <- Inf

  #coordinate descent
  if (algo == "bisection") {
    for (i in 1L:max_iter) {
      old_score <- score
      for (j in seq_along(B)) {
        #calculate the optimal lagrange multiplier for each background j
        bisection_j <- bisection_data_multiple(B[[j]],
                                               t_A,
                                               t_B,
                                               right,
                                               svd_right,
                                               lambda,
                                               j,
                                               ...)
        lambda[j] <- bisection_j$tau
      }
      score <- sum(bisection_j$values, lambda)
      if (abs(old_score - score) < tol * abs(old_score)) {
         break
      }
    }
  } else if (algo == "cd") {
    for (i in 1L:max_iter) {
      old_score <- score
      for (j in seq_along(B)) {
        #calculate the optimal lagrange multiplier for each background j
        bisection_j <- optim_cd_multiple(B[[j]],
                                         t_A,
                                         t_B,
                                         right,
                                         svd_right,
                                         lambda,
                                         j,
                                         ...)
        lambda[j] <- bisection_j$tau
      }
      score <- sum(bisection_j$values, lambda)
      if (abs(old_score - score) < tol * abs(old_score)) {
         break
      }
    }
  } else {
    stop(paste("algo", algo, " not recognized"))
  }
    left <- do.call(cbind, c(list(t_A), Map("*", -lambda, t_B)))
    final_res <- broken_svd_cpp(left, right, nv)

    return(list(values = final_res$values,
                vectors = final_res$vectors,
                tau = lambda))
}

