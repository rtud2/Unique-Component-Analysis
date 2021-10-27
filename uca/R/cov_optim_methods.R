#' objective function for BFGS implementation using covariance matrices
#'
#' @param A Target Covariance Matrix
#' @param B list of background covariance(s).
#' @param tau lowerbound for root finding
#' @return value of objective function
#' @importFrom methods as
#' @importFrom RSpectra eigs_sym
obj_fun <- function(A, B, tau) {
  eigen_calc <- eigs_sym(A - tau * B, 1L, "LA")
  eigen_calc$values + tau
}

#' gradient function for BFGS implementation using covariance matrices
#'
#' @param A Target Covariance Matrix
#' @param B list of background covariance(s).
#' @param tau lowerbound for root finding
#' @return value of gradient function
#' @importFrom methods as
#' @importFrom RSpectra eigs_sym
gr_fun <- function(A, B, tau) {
  eigen_calc <- eigs_sym( A - tau * B, 1L, "LA")
   1 - as(crossprod(eigen_calc$vectors, B %*% eigen_calc$vectors), "numeric")
}

#' BFGS Implementation using covariance matrices for single background
#'
#' Compute the UCA dimensions for single background using optim-L-BFGS-B to
#' find the optimal Lagrange multiplier
#' @param A Target Covariance Matrix
#' @param B Background Covariance Matrix.
#' @param nv number of eigenvectors to use
#' @param maxit maxium number of iterations for the algorithm to run
#' @return list of two elements:
#' \itemize{
#'  \item values: optimal eigenvalue
#'  \item tau: the contrastive parameter (lagrange multiplier)
#'  }
optim_cov_cd <- function(A, B, maxit = 5E2L, nv = 1) {

optim_with_grad <- optim(par = 2,
                        fn = obj_fun,
                        gr = gr_fun,
                        A = A,
                        B = B,
                        method = "L-BFGS-B",
                        lower = 0,
                        control = list(maxit = maxit))
if (optim_with_grad$convergence != 0) {
    warning("optim did not converge \n")
}
tau <- optim_with_grad$par

return(list(values = optim_with_grad$value - tau, tau = tau))
}

#' objective function for gradient descent implementation in multi-background
#' with covariance matrices
#'
#' @param A Target Covariance Matrix
#' @param B_unlist B unlist of background covariance(s).
#' @param tau lowerbound for root finding
#' @return value of objective function
#' @importFrom methods as
#' @importFrom RSpectra eigs_sym
obj_fun_multiple_gd <- function(A, B_unlist, tau) {
    cov_length <- ncol(A)^2
    B_obj_global <- lapply(seq_along(tau),
                           function(iters) {
        matrix(B_unlist[((iters - 1) * cov_length) + (1:cov_length)],
               ncol = ncol(A))
      })
    contrastive_mat <- A - Reduce("+", Map("*", tau, B_obj_global))
    eigen_calc_global <- eigs_sym(contrastive_mat, 1L, "LA")
    eigen_calc_global$values + sum(tau)
}


#' Gradient Descent with BFGS implementation for single OR multi-background
#' with covariance matrices
#'
#' Use L-BFGS-B to find the optimal Lagrange Multipliers for the optimal
#' Lagrangian for a single, or multiple backgrounds, simulatneously, using
#' covariance matrices.
#' @param A Target Covariance Matrix
#' @param B List of background Covariance Matrices
#' @param nv number of eigenvectors to use
#' @param maxit maxium number of iterations for the algorithm to run
#' @return list of two elements:
#' \itemize{
#'  \item values: optimal eigenvalue
#'  \item tau: the contrastive parameter(s) (lagrange multiplier(s)) for each
#'  background in the multi. background setting
#'  }
optim_cov_bfgs_gd <- function(A, B, maxit = 5E2L, nv = 1) {
    n_bg <- length(B)
    gr_dsc_mthd <- optim(par = rep(2, n_bg),
                     fn = obj_fun_multiple_gd,
                     A = A,
                     B_unlist = unlist(B),
                     method = "L-BFGS-B",
                     lower = rep(0, n_bg),
                     control = list(maxit = maxit))

if (gr_dsc_mthd$convergence != 0) {
    warning("optim did not converge \n")
}
tau <- gr_dsc_mthd$par

return(list(values = gr_dsc_mthd$value - tau,
            tau = tau))
}
