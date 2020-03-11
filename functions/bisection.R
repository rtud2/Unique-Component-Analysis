#' score_calc: calculating the derivative for specific value of lagrangian
#'
#' Calculate the derivative of Lagrangian. depends on irlba
#' 
#' @param A Target Covariance Matrix
#' @param B list of background covariance(s).
#' @param tau lowerbound for root finding
#' @return list of tau (the lagrangian), eigenvector associated with tau, and derivative value of the lagrange multiplier


score_calc = function(A, B, tau){
  eigen_calc <- partial_eigen ( A - tau * B, n = 1, symmetric = T)
  v <- matrix(eigen_calc$vectors, ncol = 1)
  score = 1 - crossprod(v, B %*% v) #derivative of  sum(diag(crossprod(v, (A - tau * B) %*% v))) + tau
  return(list(score = score, values = eigen_calc$values, vectors = eigen_calc$vectors, tau = tau))
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
#' @return the final list of tau (the optimal lagrange multiplier), vector (eigenvector)
#'  associated with tau, value (eigenvalue), and score (derivative value of the lagrangian)
 
bisection = function(A, B, limit = c(0,100), maxit = 1E5, tol = 1E-6, checks = F){
  
  if(checks == T){
    if(sign(score_calc(A, B, tau = upper)$score) + sign(score_calc(A, B, nv = nv, tau = lower)$score) !=0){
      stop("bisection will not work with the initial lower and upper values chosen.")
    }
    if(upper < lower){
      stop("upper value must be larger than lower value")
    }
  }

  for(iter in 1:maxit){
    if( max(limit) - min(limit) < tol * min(limit)) break;
    
    tau_score = score_calc(A, B, tau = sum(limit)/2)
    if(tau_score$score < 0){
      limit[1] = tau_score$tau
    }else{
      limit[2] = tau_score$tau
    }
  }
  return(dca = tau_score)
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
#' @return for each background covariance matrix in B, return list of tau (the optimal lagrange multiplier), vectors(eigenvector)
#'  associated with largest value (eigenvalue)
 
bisection.multiple = function(A, B, lambda=NULL, nv = 2, max_iter = 1000, tol = 1E-6, ...){
  #initialize starting point if one isn't supplied
  if(length(lambda) == 0){lambda = rep(10, length(B))}
  score = Inf; 
  
    for(i in 1:max_iter){
      old.score <- score
      for (j in 1:length(B)){
        A_star = A - Reduce("+", Map("*", lambda[-j], B[-j]))
        bisection_j <- bisection(A_star, B[[j]], ...) 
        lambda[j] = bisection_j$tau
      }
      score <- bisection_j$values + sum(lambda)
      if(abs(old.score - score) < tol * old.score) break;
    }
  dca <- partial_eigen(A - Reduce("+", Map("*", lambda, B)), n = nv, symmetric = T)
 return(list(values = dca$values, vectors = dca$vectors, tau = lambda))
}
  

#' dca: discriminant component analysis
#'
#' Run discriminant component analysis
#' 
#' @param A Target Covariance Matrix
#' @param B list of background covariance(s). 
#' @param ... other parameters to pass in to bisection(...) and bisection.multiple(...)
#' @return values(eigenvalues), vectors (eigenvectors), tau (Lagrange Multiplier) associated with discriminant component analysis


dca = function(A, B, ...){
  if(!is.list(B) & !is.matrix(B)){
    stop("B is not a list of matrices or a matrix")
  }
  
  if(is.list(B) & length(B) > 1){
    bisection.multiple(A, B, ... )
  }else{
    bisection(A, B, ...)
  }
}
