#' score_calc: calculating the derivative for specific value of lagrangian
#'
#' Calculate the derivative of Lagrangian. depends on irlba
#' 
#' @param A Target Covariance Matrix
#' @param B list of background covariance(s).
#' @param nv number of eigenvectors to use
#' @param tau lowerbound for root finding
#' @return list of tau (the lagrangian), eigenvector associated with tau, and derivative value of the lagrange multiplier


score_calc = function(A, B, nv, tau){
  eigen_calc <- partial_eigen ( A - tau * B, n = nv, symmetric = T)
  v <- matrix(eigen_calc$vectors, ncol = nv)
  score = 1 - sum(diag(crossprod(v, B)%*% v)) #derivative of  sum(diag(crossprod(v, (A - tau * B) %*% v))) + tau
  return(list(vector = v, score = score, tau = tau, nv = nv))
  }

#' bisection
#'
#' Use bisection method to find the optimal Lagrangian for solving discriminant PCA. depends on irlba
#' 
#' @param A Target Covariance Matrix
#' @param B Background Covariance Matrix. 
#' @param nv number of eigenvectors to use
#' @param limit a vector of c(lower, upper) bounds
#' @param maxit maxium number of iterations for the algorithm to run
#' @param tol tolerance for when to stop the algorithm
#' @return the final list of tau (the optimal lagrangian), eigenvector associated with tau, and derivative value of the lagrange multiplier
 
bisection = function(A, B, nv = 1, limit = c(0,100), maxit = 1E5, tol = 1E-4, checks = F){
  
  if(checks == T){
    if(sign(score_calc(A, B, nv = nv, tau = upper)$score) + sign(score_calc(A, B, nv = nv, tau = lower)$score) !=0){
      stop("bisection will not work with the initial lower and upper values chosen.")
    }
    if(upper < lower){
      stop("upper value must be larger than lower value")
    }
  }

  for(iter in 1:maxit){
    if( max(limit) - min(limit) < tol) break;
    
    tau_score = score_calc(A, B, nv = nv, tau = sum(limit)/2)
    if(tau_score$score < 0){
      limit[1] = tau_score$tau
    }else{
      limit[2] = tau_score$tau
    }
  }
  return(tau_score)
}


#calculate A*, let j denote which bg data you're working on
# B list of matrix
# lambda <- 1:length(B)

#' bisection.multiple: bisection method but for multiple background datasets
#'
#' Use bisection method to find the optimal Lagrangian for solving discriminant PCA of each background dataset. depends on irlba
#' 
#' @param A Target Covariance Matrix
#' @param B list of background covariance(s).
#' @param ... other parameters to pass in to bisection(...)
#' @return for each background covariance matrix in B, return list of lambda (the optimal lagrangian), eigenvector associated with largest eigenvalue, and derivative value of the lagrange multiplier
#' 
 
bisection.multiple = function(A, B, lambda, max_iter = 1000, tol = 1E-5, ...){
 if(!is.list(B) & !is.matrix(B)){
   stop("B is not a list of matrices or a matrix")
 }else if(is.matrix(B)){
   warning("B is being coerced into a list")
   B <- list(B)
 }
  score = 0; #initialize
    for(i in 1:max_iter){
      old.score <- score
      for (j in 1:length(B)){
        A_star = A - Reduce("+", mapply("*", lambda[-j], B[[-j]], SIMPLIFY = F))
        lambda[j] = bisection(A_star, B[[j]], ...)$tau
      }
      score <- partial_eigen(A - Reduce( "+", mapply( "*", lambda, B, SIMPLIFY = F)), n = 1, symmetric = T)$values + sum(lambda) 
      if(score - old.score < tol) break;
    }
 return(lambda) 
}
  

#' dca: discriminant component analysis
#'
#' Run discriminant component analysis
#' 
#' @param A Target Covariance Matrix
#' @param B list of background covariance(s). 
#' @param n number of eigenvectors to keep from dca
#' @param ... other parameters to pass in to bisection(...)
#' @return eigenvalues and eigenvectors associated with discriminant component analysis


dca_f = function(A, B, n = 2, ...){
    lagrangian_res <- bisection.multiple(A, B, ...)  
  }
  if(is.list(B)){
    # mapply multiplies the lagrangian found by bisection to the background matrices. Reduce adds them together
    sigma = A -  Reduce("+", mapply("*", lapply(lagrangian_res, "[[",3), B, SIMPLIFY = FALSE))
  } else{
    sigma = A -  sapply(lagrangian_res, "[[",3) * B # when B is a matrix
  }
  
  partial_eigen(sigma, n = n, symmetric = T)
}

