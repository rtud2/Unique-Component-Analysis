#' bisection: bisection's method
#'
#' Use bisection method to find the optimal Lagrangian for solving discriminant PCA. depends on irlba
#' 
#' @param A dataset - dataset of interest
#' @param B Background dataset(s). Input multiple datasets as a list. 
#' @param nv number of eigenvectors to use
#' @param lower lowerbound for root finding
#' @param upper upperbound for root finding
#' @param tol tolerance for when to stop the algorithm
#' @return tau, the optimal lagrangian


score_calc = function(A, B, nv, tau){
  eigen_calc <- partial_eigen ( A - tau * B, n = nv, symmetric = T)
  v <- matrix(eigen_calc$vectors, ncol = nv)
  score = 1 - sum(diag(crossprod(v, B %*% v))) #derivative of  sum(diag(crossprod(v, (A - tau * B) %*% v))) + tau
  return(
    list( score = score, tau = tau)) #v = v))
}


bisection = function(A, B, nv = 1, lower = 0, upper = 100, tol = 1E-5, lower_score = NULL, upper_score = NULL){
  
  if(length(lower_score)==0 & length(upper_score)==0 ){
    upper_score = score_calc(A, B, nv = nv, tau = upper)
    lower_score = score_calc(A, B, nv = nv, tau = lower)
    
    if(sign(upper_score$score) + sign(lower_score$score) !=0){
      stop("bisection will not work with the lower and upper values chosen")
    }
  }
  
  prev_score <- c(lower_score$score, upper_score$score)
  
  tau = (lower_score$tau + upper_score$tau)/2
  proposed_score = score_calc(A, B, nv = nv, tau = tau)
  
  to_replace = which(sign(prev_score) == sign(proposed_score$score))
  score_diff = as.numeric(proposed_score$score - prev_score[-to_replace])
  
  if(abs(score_diff) > tol & to_replace == 1){  
    bisection(A, B, lower_score = proposed_score, upper_score = upper_score)
  }else if(abs(score_diff) > tol & to_replace == 2){
    bisection(A, B, lower_score = lower_score, upper_score = proposed_score)
  }else{
      return(tau)
  }
}

