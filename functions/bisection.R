#' score_calc: calculating the derivative for specific value of lagrangian
#'
#' Calculate the derivative of Lagrangian. depends on irlba
#' 
#' @param A dataset - dataset of interest
#' @param B Background dataset(s). Input multiple datasets as a list. 
#' @param nv number of eigenvectors to use
#' @param tau lowerbound for root finding
#' @return list of tau and value of derivative, the optimal lagrangian


score_calc = function(A, B, nv, tau, return_all = F){
  eigen_calc <- partial_eigen ( A - tau * B, n = nv, symmetric = T)
  v <- matrix(eigen_calc$vectors, ncol = nv)
  score = 1 - sum(diag(crossprod(v, B)%*% v)) #derivative of  sum(diag(crossprod(v, (A - tau * B) %*% v))) + tau
  if(return_all){
    list(vector = v, score = score, tau = tau, nv = nv)
  }else{
    return(score)
    
  }
  }

#' bisection
#'
#' Use bisection method to find the optimal Lagrangian for solving discriminant PCA. depends on irlba
#' 
#' @param A dataset - dataset of interest
#' @param B Background dataset(s). 
#' @param nv number of eigenvectors to use
#' @param limit a vector of c(lower, upper) bounds
#' @param tol tolerance for when to stop the algorithm
#' @return tau, the optimal lagrangian

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
    
    tau = sum(limit)/2
    if(score_calc(A, B, nv = nv, tau = tau) < 0){
      limit[1] = tau
    }else{
      limit[2] = tau
    }
  }
  final = score_calc(A, B, nv = nv, tau = tau, return_all = T) 
   
  return(final)
}

#' bisection.multiple: bisection method but for multiple background datasets
#'
#' Use bisection method to find the optimal Lagrangian for solving discriminant PCA of each background dataset. depends on irlba
#' 
#' @param A dataset - dataset of interest
#' @param B list of background dataset(s). 
#' @param ... other parameters to pass in to bisection(...)
#' @return tau, the optimal lagrangian


bisection.multiple = function(A, B, ...){
 if(!is.list(B)){
   stop("B is not a list of matrices")
 }else if(!is.list(A)){
   warning(paste0("replicating A ", length(B), " times to match length(B)"))
   A = replicate(n = length(B),expr =  A, simplify = F)
 }
  
  if( length(list(...)) > 0){
    mapply(FUN = bisection, A, B, list(...), SIMPLIFY = F)
  }else{
    mapply(FUN = bisection, A, B, SIMPLIFY = F)
  }
}

