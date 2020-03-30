#' score_calc: calculating the derivative for specific value of lagrangian
#'
#' Calculate the derivative of Lagrangian. depends on RSpectra
#' 
#' @param A Target Covariance Matrix
#' @param B list of background covariance(s).
#' @param tau lowerbound for root finding
#' @return list of tau (the lagrangian), eigenvector associated with tau, and derivative value of the lagrange multiplier


score_calc = function(A, B, tau){
  eigen_calc <- eigs_sym ( A - tau * B, 1, "LA")
  v <- matrix(eigen_calc$vectors, ncol = 1)
  score = as.numeric(1 - crossprod(v, B %*% v)) #derivative of  sum(diag(crossprod(v, (A - tau * B) %*% v))) + tau
  return(list(score = score, values = eigen_calc$values, tau = tau))
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

bisection = function(A, B, limit = c(0,20), maxit = 1E5, nv = 1, tol = 1E-6){
  
  f_val = lapply(limit, function(z) score_calc(A, B, z))
  maxit_flag = F;
  initial <- sapply(f_val, "[[", 1)
  
  if(sum(sign(initial)) == 2){
    warning("initial f(a) and f(b) are positive. Setting lambda to 0 \n")
    return(f_val[[1]])
  }else{
    
    for(iter in 1:maxit){
      limit <- sapply(f_val, "[[", 3)
      if( max(limit) - min(limit) < tol * min(limit)) break;
      if(iter == maxit) maxit_flag == T;
      
      tau_score = score_calc(A, B, tau = sum(limit)/2)
      
      if(tau_score$score < 0){
        f_val[[1]] = tau_score
      }else{
        f_val[[2]] = tau_score
      }
    }  
    
    if(maxit_flag == T){
      warning("maximum iteration reached: solution may not be optimal \n")
    }
    min_fval <- which.min(abs(sapply(f_val, "[[", 1)))
    return(f_val[[min_fval]]) 
  }
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
 
bisection.multiple = function(A, B, lambda=NULL, nv = 2, max_iter = 1E5, tol = 1E-6, ...){
  
  #initialize starting point if one isn't supplied
  if(length(lambda) == 0){lambda = sapply(1:length(B), function(zz){bisection(A,B[[zz]])$tau})}
  score = Inf; 
  
    for(i in 1:max_iter){
      old.score <- score
      for (j in 1:length(B)){
        A_star = A - Reduce("+", Map("*", lambda[-j], B[-j]))
        bisection_j <- bisection(A_star, B[[j]], ...) 
        lambda[j] = bisection_j$tau
      }
      score <- bisection_j$values + sum(lambda)
      if(abs(old.score - score) < tol * abs(old.score)) break;
    }
  dca <- eigs_sym(A - Reduce("+", Map("*", lambda, B)), nv, "LA")
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


dca = function(A, B, nv = 2, ...){
  if(!is.list(B) & !is.matrix(B)){
    stop("B is not a list of matrices or a matrix")
  }
  
  if(is.list(B) & length(B) > 1){
    bisection.multiple(A=A, B=B, nv=nv, ... )
  }else{
    if(is.list(B)) B = B[[1]]
    tmp_res <- bisection(A=A, B=B, nv=nv, ...)
    res <- eigs_sym(A - tmp_res$tau*B, nv, "LA")
    return(list(values = res$values, vectors = res$vectors, tau = tmp_res$tau))
  }
}
