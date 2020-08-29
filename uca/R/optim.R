#' objective function for optim implementation, cov method
#' 
#' @param A Target Covariance Matrix
#' @param B list of background covariance(s).
#' @param tau lowerbound for root finding
#' @return value of objective function 
#' @importFrom methods as
#' @importFrom RSpectra eigs_sym
obj_fun <- function(A, B, tau){
  eigen_calc <- eigs_sym ( A - tau * B, 1L, "LA")
  eigen_calc$values + tau
}

#' gradient function for optim implementation, cov method
#' 
#' @param A Target Covariance Matrix
#' @param B list of background covariance(s).
#' @param tau lowerbound for root finding
#' @return value of gradient function 
#' @importFrom methods as
#' @importFrom RSpectra eigs_sym
gr_fun = function(A, B, tau){
  eigen_calc <- eigs_sym ( A - tau * B, 1L, "LA")
   1 - as( crossprod(eigen_calc$vectors, B %*% eigen_calc$vectors), "numeric")
}

#' optim_bfgs implementaiton, cov method
#' 
#' Use optim-L-BFGS-S to find the optimal Lagrangian for solving unique component analysis (uca) 
#' @param A Target Covariance Matrix
#' @param B Background Covariance Matrix. 
#' @param nv number of eigenvectors to use
#' @param maxit maxium number of iterations for the algorithm to run
#' @return the final list of tau (the optimal lagrange multiplier), vector (eigenvector)
#'  associated with tau, value (eigenvalue), and score (derivative value of the lagrangian)

optim_bfgs = function(A, B, maxit = 5E2L, nv = 1){
  
optim_with_grad = optim(par = 2,
                        fn = obj_fun,
                        gr = gr_fun,
                        A = A,
                        B = B,
                        method = "L-BFGS-B",
                        lower = 0,
                        control = list(maxit = maxit))
if(optim_with_grad$convergence != 0){warning("optim did not converge \n")}
tau = optim_with_grad$par

return(list(values = optim_with_grad$value - tau, tau=tau))
}

#' optim_bfgs2: optim implementation for the data method
#' 
#' compute the UCA for single background using SVD and QR. good for bigger data where loading the covariance matrix is difficult using optim then bisection
#' @param A a n1xp data matrix
#' @param B a n2xp data matrix
#' @param maxit maximum iterations
#' @return list of tau, largest eigenvalue, and score
#' @importFrom Rfast transpose 

optim_bfgs2 = function(A, B, maxit = 5E2L, nv = 1){
  right =  rbind(A, B)
  svd_right <- arma_svd(right)
  t_A = Rfast::transpose(A)
  t_B = Rfast::transpose(B)
  optim_with_grad = optim(par = 2,
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
  if(optim_with_grad$convergence !=0){warning("optim did not converge\n")}
  
  return(list(values = optim_with_grad$value - optim_with_grad$par, tau=optim_with_grad$par ))
}

#' obj_fn_multiple
#' 
#' objective function multiple backgroundsimplementation, data method
#' @param B background interested in calculating lagrange multiplier for
#' @param t_A transpose of A, a pxn1 data matrx
#' @param t_B list of transpose of B, a pxn2 data matrix
#' @param right rbind(A, B), a (n1+n2)xp data matrix
#' @param svd_right svd of right data matrix object
#' @param lambda the lagrange multiplier
#' @param j which background matrix are we focusing on
#' @return value of the objective function
#' @importFrom Rfast transpose 
#' 
obj_fn_multiple = function(lambda, B, t_A, t_B, lambda_B, right, svd_right, j){
  lambda_B[[j]] =  -lambda*t_B[[j]]
  left <- do.call(cbind, c(list(t_A), lambda_B))
  objective_fn <- multiple_score_calc_cpp(left, right, svd_right$u, svd_right$d, B, lambda)
  return( objective_fn$values+lambda)
}

#' gr_fn_multiple
#' 
#' gradient function multiple backgroundsimplementation, data method
#' @param B background interested in calculating lagrange multiplier for
#' @param t_A transpose of A, a pxn1 data matrx
#' @param t_B list of transpose of B, a pxn2 data matrix
#' @param right rbind(A, B), a (n1+n2)xp data matrix
#' @param svd_right svd of right data matrix object
#' @param lambda the lagrange multiplier
#' @param j which background matrix are we focusing on
#' @return value of the gradient fn
#' @importFrom Rfast transpose 

gr_fn_multiple = function(lambda, B, t_A, t_B, lambda_B, right, svd_right, j){
  lambda_B[[j]] =  -lambda*t_B[[j]]
  left <- do.call(cbind, c(list(t_A), lambda_B))
  objective_gr <- multiple_score_calc_cpp(left, right, svd_right$u, svd_right$d, B, lambda)
  return( objective_gr$score)
}

#' optim_bfgs_multiple
#' 
#' objective function multiple backgroundsimplementation, data method
#' @param B_focus single background interested in calculating lagrange multiplier for
#' @param t_A transpose of A, a pxn1 data matrx
#' @param t_B list of transpose of B, a pxn2 data matrix
#' @param right rbind(A, B), a (n1+n2)xp data matrix
#' @param svd_right svd of right data matrix object
#' @param lambda the lagrange multiplier
#' @return optimum eigenvalue and lagrange multiplier for a single background in the multi. background setting
#' @importFrom Rfast transpose
optim_bfgs_multiple = function(B_focus, t_A, t_B, right, svd_right, lambda, j, maxit = 5E2){
  
  #constants that don't really change if focused on j-th background
  old_right <- Rfast::transpose(do.call(cbind, c(list(t_A),  t_B[-j]))) #for some reason, faster than rbind due to memory allocation.
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
  if(f_val$score >= 0){
    warning("Redundant Constraint: Lagrange Multiplier is negative. Setting lambda to 0 \n");
    return(f_val);
  }else{
  optim_res <- optim(par = 3,
        fn = obj_fn_multiple,
        gr = gr_fn_multiple,
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
  if(optim_res$convergence !=0){warning("optim did not converge\n")}

return(list(values = optim_res$value - optim_res$par, tau=optim_res$par ))
}




